subroutine shock_tube_ion
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,col,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,flag_ir,&
       x,y,z,dx,dy,dz,n_fraction,debug_direction,debug_option,T0,n0,Nexcite,&
       f_p_ini,f_p_p_ini,n0fac,Gm_rec_ref,colrat,gm_ion,gm_rec,expinttab,&
       flag_rad,radrat,rad_temp,gm_rec_rad,gm_ion_rad
  use scheme_rot,only:pv2cq_mhd,pv2cq_hd,get_Te_MHD
  use model_rot, only:set_coordinate,setcq
  use PIP_rot, only:get_col_ion_coeff,expintread,get_radrat_fixed
  implicit none
  double precision :: ro_h(1:ix,1:jx,1:kx),ro_m(1:ix,1:jx,1:kx)
  double precision :: vx_h(1:ix,1:jx,1:kx),vx_m(1:ix,1:jx,1:kx)
  double precision :: vy_h(1:ix,1:jx,1:kx),vy_m(1:ix,1:jx,1:kx)
  double precision :: vz_h(1:ix,1:jx,1:kx),vz_m(1:ix,1:jx,1:kx) 
  double precision :: P_h (1:ix,1:jx,1:kx),P_m (1:ix,1:jx,1:kx)
  double precision :: B_x (1:ix,1:jx,1:kx)
  double precision :: B_y (1:ix,1:jx,1:kx)
  double precision :: B_z (1:ix,1:jx,1:kx)
  double precision ::mask(ix,jx,kx)
  double precision ::v_para(ix,jx,kx),v_perp(ix,jx,kx)
  double precision ::b_para(ix,jx,kx),b_perp(ix,jx,kx)
  double precision ::v_b_para(ix,jx,kx),b_b_para(ix,jx,kx)
  double precision ::v_theta(ix,jx,kx),v_phi(ix,jx,kx)
  double precision ::b_theta(ix,jx,kx),b_phi(ix,jx,kx)
  double precision ::Te_p(ix,jx,kx),Gm_ion0(ix,jx,kx),Gm_rec0(ix,jx,kx)
  double precision::Eion(6) !Energy to ionise
  double precision f_n,f_p,f_p_n,f_p_p,start(3),end(3),B0
  double precision theta_p,phi_p,tmp,v_L(8),v_R(8),wtr,ioneq,Te_0,tfac,n1(6),tstep
  double precision,parameter::kbhat=1.38064852,mehat=9.10938356,hhat=6.62607004
  double precision,parameter::kboltz=1.38064852e-23 !Boltzmann Constant [m^2 kg s^-2 K^-1]
  integer i,j,k,ti

  !Find the equilibrium neutral fraction
  if(flag_IR .eq. 4) then
      print*,'Calculating LTE excitation state'
      allocate(Nexcite(ix,jx,kx,6)) !Allocate the fractional array
      Eion=[13.6,3.4,1.51,0.85,0.54,0.0] !in eV
      Eion=Eion/13.6*2.18e-18 !Convert to joules (to be dimensionally correct)
      Nexcite(:,:,:,6)=n0
      Nexcite(:,:,:,1)=(2.0/n0/2.d0*(2.0*pi*mehat*kbhat*T0/hhat/hhat*1.0e14)**(3.0/2.0)*exp(-Eion(1)/kboltz/T0))
      Nexcite(:,:,:,2)=(2.0/n0/8.d0*(2.0*pi*mehat*kbhat*T0/hhat/hhat*1.0e14)**(3.0/2.0)*exp(-Eion(2)/kboltz/T0))
      Nexcite(:,:,:,3)=(2.0/n0/18.d0*(2.0*pi*mehat*kbhat*T0/hhat/hhat*1.0e14)**(3.0/2.0)*exp(-Eion(3)/kboltz/T0))
      Nexcite(:,:,:,4)=(2.0/n0/32.d0*(2.0*pi*mehat*kbhat*T0/hhat/hhat*1.0e14)**(3.0/2.0)*exp(-Eion(4)/kboltz/T0))
      Nexcite(:,:,:,5)=(2.0/n0/50.d0*(2.0*pi*mehat*kbhat*T0/hhat/hhat*1.0e14)**(3.0/2.0)*exp(-Eion(5)/kboltz/T0))
      Nexcite(:,:,:,1:5)=n0/Nexcite(:,:,:,1:5)
      Nexcite=Nexcite/n0
      Nexcite=Nexcite/sum(Nexcite(1,1,1,:))
print*,'finding equilibrium'
      !print*,Nexcite(1,1,1,:)
allocate(expinttab(4,10000)) !table for the exponential integral table
allocate(Colrat(ix,jx,kx,6,6))
if (flag_rad .eq. 1) allocate(radrat(ix,jx,kx,6,6))
call expintread
      call get_col_ion_coeff(T0+0.d0*U_m(:,:,:,1),n0+0.d0*U_m(:,:,:,1),Gm_ion,Gm_rec) 
!      call get_col_ion_coeff_aprox(spread(spread(spread(T0,1,ix),2,jx),3,kx),&
!            spread(spread(spread(n0,1,ix),2,jx),3,kx),Gm_ion,Gm_rec) 
!print*,T0,n0
!print*,colrat(1,1,1,:,:)
!stop
!print*,'Colrat'
!print*,colrat(1,1,1,:,:)
      tstep=0.1*minval(1.0/colrat(1,1,1,:,:))
      if (flag_rad .eq. 1) then
        call get_radrat_fixed(rad_temp,T0+0.d0*U_m(:,:,:,1),n0+0.d0*U_m(:,:,:,1),Gm_ion_rad,Gm_rec_rad)
        tstep=min(tstep,0.1*minval(1.0/radrat(1,1,1,:,:)))
!print*,'Radrat'
!print*,radrat(1,1,1,:,:)
      endif
      do ti=0,50000 !Should be a while loop
        !print*,nexcite(1,1,1,1:6)/sum(nexcite(1,1,1,1:6))!,colrat
        n1=nexcite(1,1,1,:)
        do i=1,6
        do j=1,6
            nexcite(1,1,1,i)=nexcite(1,1,1,i)+tstep*(n1(j)*colrat(1,1,1,j,i) - n1(i)*colrat(1,1,1,i,j)) !colisional step
            if (flag_rad .eq. 1) then
                nexcite(1,1,1,i)=nexcite(1,1,1,i)+tstep*(n1(j)*radrat(1,1,1,j,i) - n1(i)*radrat(1,1,1,i,j)) !radiative step
            endif
        enddo
            nexcite(:,:,:,i)=max(nexcite(1,1,1,i),0.0) !prevent negative masses
        enddo
    enddo
      f_n=sum(Nexcite(1,1,1,1:5))/sum(Nexcite(1,1,1,1:6))
      f_p=1.d0-f_n!Nexcite(1,1,1,6)
      f_p_n=f_n/(f_n+2.0d0*f_p)
      f_p_p=2.0d0*f_p/(f_n+2.0d0*f_p)
      f_p_ini=f_p
      f_p_p_ini=f_p_p
      !n0fac=Nexcite(1,1,1,6)
      n0fac=f_p
!print*,Nexcite(1,1,1,:)
!stop
      print*,f_n,f_p
      if (f_n*f_p .le. 0) then
            print*,'Something wrong with initial conditions'
            stop
      endif     !stop
  else 
      Te_0=T0/1.1604e4
      ioneq=(2.6e-19/dsqrt(Te_0))/(2.91e-14/(0.232+13.6/Te_0)*(13.6/Te_0)**0.39*dexp(-13.6/Te_0))
      f_n=ioneq/(ioneq+1.0d0)
      f_p=1.0d0-f_n
      f_p_n=f_n/(f_n+2.0d0*f_p)
      f_p_p=2.0d0*f_p/(f_n+2.0d0*f_p)
  endif

  if (my_rank.eq.0) print*,'Neutral fraction = ',f_n
!stop
!  !set ionization fraction-----------------
  if(flag_pip.eq.0) then
     f_n=1.0d0
     f_p=1.0d0
     f_p_n=1.0d0
     f_p_p=1.0d0
  endif
  !----------------------------------------

  !Set coordinate (uniform grid)--------------------------
  !!set lower and upper coordinate
  start(1)=0.0d0 ;end(1)=1500.0d0          !400.0d0/f_p
  start(2)=-1.0d0 ;end(2)=1.0d0
  start(3)=-1.0d0 ;end(3)=1.0d0
  call set_coordinate(start,end)
  !---------------------------------------
  
  !!default boundary condition----------------------
  if (flag_bnd(1) .eq.-1) flag_bnd(1)=3
  if (flag_bnd(2) .eq.-1) flag_bnd(2)=10
  if (flag_bnd(3) .eq.-1) flag_bnd(3)=10
  if (flag_bnd(4) .eq.-1) flag_bnd(4)=10
  if (flag_bnd(5) .eq.-1) flag_bnd(5)=10
  if (flag_bnd(6) .eq.-1) flag_bnd(6)=10
  !-------------------------------------------------

  !!!========================================================
  !parameters---------------------
!  B0=sqrt(2.0*gm/beta)
  B0=1.0d0
!  B0=0.0d0
!/sqrt(2.0)
  !-------------------------
  if(debug_direction.eq.1) then
     mask=spread(spread(pi*x,2,jx),3,kx)
     phi_p=0.0d0
     theta_p=pi/2.0d0
     tmp=0.0
  else if(debug_direction.eq.2) then
     mask=spread(spread(pi*y,1,ix),3,kx)     
     phi_p=pi/2.0d0
     theta_p=pi/2.0d0     
     tmp=0.0
  else if(debug_direction.eq.3) then
     mask=spread(spread(pi*z,1,jx),1,ix)
     phi_p=0.0d0
     theta_p=0.0d0
     tmp=0.0
  else if(debug_direction.eq.4) then
     do k=1,kx;do j=1,jx; do i=1,ix
        if(x(i)+y(j).gt.0)then
           mask(i,j,k)=1
        else
           mask(i,j,k)=-1
        endif
     enddo;enddo;enddo
     phi_p=pi/4.0d0
     theta_p=pi/2.0d0
     tmp=pi/2.0d0
  else if(debug_direction.eq.5) then
     do k=1,kx;do j=1,jx; do i=1,ix
        if(x(i)+z(k).gt.0)then
           mask(i,j,k)=1
        else
           mask(i,j,k)=-1
        endif
     enddo;enddo;enddo
     phi_p=0.0d0
     theta_p=pi/4.0d0
     tmp=0.0
  else if(debug_direction.eq.6) then
     do k=1,kx;do j=1,jx; do i=1,ix
        if(y(j)+z(k).gt.0)then
           mask(i,j,k)=1
        else
           mask(i,j,k)=-1
        endif
     enddo;enddo;enddo
     phi_p=pi/2.0d0
     theta_p=pi/4.0d0
     tmp=0.0
  else if(debug_direction.eq.7) then
     do k=1,kx;do j=1,jx; do i=1,ix
        if(x(i)+y(j)+z(k).gt.0)then
           mask(i,j,k)=1
        else
           mask(i,j,k)=-1
        endif
     enddo;enddo;enddo
     phi_p=pi/4.0d0
     theta_p=pi/4.0d0
     tmp=0.0
  endif  
  select case(debug_option)
  case (2)  !Switch-off shock (Falle et al.1998)
     v_l=(/1.368d0,1.769d0,0.269d0,1.0d0,0.0d0,1.0d0,0.0d0,0.0d0/)
     v_r=(/1.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,1.0d0,0.0d0/)     
  case default
!     v_l=(/1.0d0,1.0d0,0.0d0,0.0d0,0.0d0,B0*0.75d0,B0,0.0d0/)
!     v_r=(/0.125d0,0.1d0,0.0d0,0.0d0,0.0d0,B0*0.75d0,-B0,0.0d0/)
     v_l=(/1.0d0,beta*B0**2/2.d0,0.0d0,0.0d0,0.0d0,B0*0.3d0,B0,0.0d0/)
!     v_r=(/1.0d0,beta*B0**2/2.d0,0.0d0,0.0d0,0.0d0,B0*0.3d0,B0,0.0d0/)
     v_r=(/1.0d0,beta*B0**2/2.d0,0.0d0,0.0d0,0.0d0,B0*0.3d0,-B0,0.0d0/)
!     Density, pressure, velocity * 3, magnetic field *3


  end select

  where(mask<0)
     ro_h=f_n*v_l(1)
     p_h=f_p_n*v_l(2)
     ro_m=f_p*v_l(1)
     p_m=f_p_p*v_l(2)
     vx_h=v_l(3)
     vy_h=v_l(4)
     vz_h=v_l(5)
     vx_m=v_l(3)
     vy_m=v_l(4)
     vz_m=v_l(5)
     b_para=v_l(6)
     b_perp=v_l(7)
  elsewhere
     ro_h=f_n*v_r(1)
     p_h=f_p_n*v_r(2)
     ro_m=f_p*v_r(1)
     p_m=f_p_p*v_r(2)
     vx_h=v_r(3)
     vy_h=v_r(4)
     vz_h=v_r(5)
     vx_m=v_r(3)
     vy_m=v_r(4)
     vz_m=v_r(5)
     b_para=v_r(6)
     b_perp=v_r(7)
  end where
!  wtr=0.001d0
  wtr=dx(1)*10.d0
!0.001d0
  do k=1,kx;do j=1,jx;do i=1,ix
     ro_h(i,j,k)=f_n*(v_l(1)+(v_r(1)-v_l(1))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
     ro_m(i,j,k)=f_p*(v_l(1)+(v_r(1)-v_l(1))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
     p_h(i,j,k)=f_p_n*(v_l(2)+(v_r(2)-v_l(2))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
     p_m(i,j,k)=f_p_p*(v_l(2)+(v_r(2)-v_l(2))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
     vx_h(i,j,k)=(v_l(3)+(v_r(3)-v_l(3))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
     vx_m(i,j,k)=(v_l(3)+(v_r(3)-v_l(3))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
     vy_h(i,j,k)=(v_l(4)+(v_r(4)-v_l(4))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
     vy_m(i,j,k)=(v_l(4)+(v_r(4)-v_l(4))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
     vz_h(i,j,k)=(v_l(5)+(v_r(5)-v_l(5))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
     vz_m(i,j,k)=(v_l(5)+(v_r(5)-v_l(5))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
     b_para(i,j,k)=(v_l(6)+(v_r(6)-v_l(6))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
     b_perp(i,j,k)=(v_l(7)+(v_r(7)-v_l(7))*(1.0d0+tanh(mask(i,j,k)/wtr))*0.5d0)
  enddo;enddo;enddo


  b_theta=b_perp*sin(tmp)
  b_phi  =b_perp*cos(tmp)
  b_x=b_para*sin(theta_p)*cos(phi_p)+ &
       b_theta*cos(theta_p)*cos(phi_p)-b_phi*sin(phi_p)
  b_y=b_para*sin(theta_p)*sin(phi_p)+ &
       b_theta*cos(theta_p)*sin(phi_p)+b_phi*cos(phi_p)
  b_z=b_para*cos(theta_p)-b_theta*sin(theta_p) 
  !!!========================================================

  !convert PV2cq and set that value to global variable 'U_h' and/or 'U_m'
  call setcq(ro_m,vx_m,vy_m,vz_m,p_m,B_x,B_y,B_z, &
       ro_h,vx_h,vy_h,vz_h,p_h)
  !---------------------------------------------------------------------

  !Set the reference recombination rate
!  if (flag_IR .eq. 4) then
!        call get_Te_MHD(U_m,Te_p)
!        tfac=beta/2.0d0*f_p_p_ini*5.0d0/6.0d0/f_p_ini
!print*,'Te_p',Te_p
!print*,'T0',T0
!print*,'tfac',tfac
!print*,U_m(:,:,:,1)
!print*,n0
!print*,n0fac
!        call get_col_ion_coeff(Te_p*T0/tfac,U_m(:,:,:,1)*n0/n0fac,Gm_ion0,Gm_rec0)
!        Gm_rec_ref=Gm_rec0(1,1,1) 
!  endif

  !set default output period and time duration--------------------------
  if(tend.lt.0.0) then
     tend=0.4d0
     dtout=tend/10.0
     if(flag_mpi.eq.0 .or.my_rank.eq.0)      print *,"TEND",dtout,tend
  endif
  !---------------------------------------------------------------------

end subroutine shock_tube_ion
