subroutine tearing
  use parameters,only:pi
  !Tearing mode test (Lundi et al.2008)
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,&
       x,y,z,dx,dy,dz,n_fraction,ndim,debug_direction,margin,xi_n,flag_amb,&
       flag_ir
  use scheme_rot,only:pv2cq_mhd,pv2cq_hd
  use model_rot, only:set_coordinate,setcq
  use PIP_rot,only:set_xin,get_initial_xin
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
  double precision ::b_para(ix,jx,kx),b_perp(ix,jx,kx),b_tang(ix,jx,kx)
  double precision ::v_b_para(ix,jx,kx),b_b_para(ix,jx,kx)
  double precision ::v_theta(ix,jx,kx),v_phi(ix,jx,kx)
  double precision ::b_theta(ix,jx,kx),b_phi(ix,jx,kx)
  double precision ::P_tot(ix,jx,kx),Te_tot(ix,jx,kx)
  double precision ::N_tot(ix,jx,kx),xi_n0(ix,jx,kx)
  double precision ::f_n(ix,jx,kx),f_p(ix,jx,kx)
  double precision ::f_p_n(ix,jx,kx),f_p_p(ix,jx,kx)
  double precision start(3),end(3)
  double precision thick,b0,Bz0
  double precision theta_p,phi_p,tmp,dTe,r_x,k_x,k_z,phase,epsilon,delta
  integer i,j,k,orix,oriy,oriz,nk,ni,ii,kk,L_x,L_z


  !----------------------------------------
  !Set coordinate (uniform grid)--------------------------
  !!set lower and upper coordinate
  L_x=70.0d0*pi
  L_z=1.0d0*pi
  start(1)=-0.5d0*L_x ;end(1)=0.5d0*L_x
!  start(1)=0.0d0 ;end(1)=4.0d0
  start(2)=0.00d0 ;end(2)=pi/2.0d0
  start(3)=-1.0d0 ;end(3)=1.0d0
  call set_coordinate(start,end)
  !---------------------------------------  
  if(flag_amb.eq.1.or.flag_pip.eq.1) then
     call set_xin(U_h,U_m)
  endif



  !set ionization fraction-----------------
  if(flag_pip.eq.0) then
     f_n=1.0d0
     f_p=1.0d0
     f_p_n=1.0d0
     f_p_p=1.0d0
  else
     f_n=xi_n
     f_p=1.0d0-xi_n
     f_p_n=f_n/(f_n+2.0d0*f_p)
     f_p_p=2.0d0*f_p/(f_n+2.0d0*f_p)
  endif

  
  !!default boundary condition----------------------
  if (flag_bnd(1) .eq.-1) flag_bnd(1)=11
  if (flag_bnd(2) .eq.-1) flag_bnd(2)=11
  if (flag_bnd(3) .eq.-1) flag_bnd(3)=11
  if (flag_bnd(4) .eq.-1) flag_bnd(4)=11
  if (flag_bnd(5) .eq.-1) flag_bnd(5)=11
  if (flag_bnd(6) .eq.-1) flag_bnd(6)=11
  !-------------------------------------------------

  !!!========================================================
  !set parameter
  thick=0.5d0
  B0=sqrt(2.0d0/(beta*gm))
  Bz0=1.0d0*B0*0
  tmp=0.0d0
  !set debug direction----------------------------------------------
  if(debug_direction.eq.1) then
     mask=spread(spread(x,2,jx),3,kx)
     phi_p=0.0d0
     theta_p=pi/2.0d0
     flag_bnd(1)=3
     flag_bnd(2)=11
!     flag_bnd(2)=2
     flag_bnd(3)=4
  else if(debug_direction.eq.2) then
     mask=spread(spread(y,1,ix),3,kx)     
     phi_p=pi/2.0d0
     theta_p=pi/2.0d0     
     flag_bnd(1)=11
     flag_bnd(2)=11
     flag_bnd(3)=4
     flag_bnd(4)=3
  else if(debug_direction.eq.3) then
     mask=spread(spread(z,1,jx),1,ix)
     phi_p=0.0d0
     theta_p=0.0d0
     flag_bnd(5)=2
     flag_bnd(6)=2
  else if(debug_direction.eq.4) then
     do k=1,kx;do j=1,jx; do i=1,ix
        mask(i,j,k)=pi*(x(i)+y(j))
     enddo;enddo;enddo
     phi_p=pi/4.0d0
     theta_p=pi/2.0d0
     tmp=pi/2.0
  else if(debug_direction.eq.5) then
     do k=1,kx;do j=1,jx; do i=1,ix
        mask(i,j,k)=(x(i)+z(k))
     enddo;enddo;enddo
     phi_p=0.0d0
     theta_p=pi/4.0d0
  else if(debug_direction.eq.6) then
     do k=1,kx;do j=1,jx; do i=1,ix
        mask(i,j,k)=(y(j)+z(k))
     enddo;enddo;enddo
     phi_p=pi/2.0d0
     theta_p=pi/4.0d0
  else if(debug_direction.eq.7) then
     do k=1,kx;do j=1,jx; do i=1,ix
        mask(i,j,k)=(x(i)+y(j)+z(k))
     enddo;enddo;enddo
     phi_p=pi/4.0d0
     theta_p=pi/4.0d0
  endif  
  !----------------------------------------------------------------
  vx_h=0.0d0
  vy_h=0.0d0
  vz_h=0.0d0
  vx_m=0.0d0
  vy_m=0.0d0
  vz_m=0.0d0  
  b_perp=B0*tanh(mask/thick)
  b_tang=Bz0  
  p_tot=(1.0d0/gm+0.5d0*(b0*b0-b_perp**2))

  if(flag_ir.eq.2) then
     te_tot=1.0d0
     r_x=10.0d0
     dTe=0.3d0
     do k=1,kx;do j=1,jx;do i=1,ix
        te_tot(i,j,k)=te_tot(i,j,k)-dTe+dTe*((tanh((x(i)-r_x)/2.0d0)+1.0d0)*0.5d0 &
             +(1.0d0-tanh((x(i)+r_x)/2.0d0))*0.5d0)
     enddo;enddo;enddo     
     call get_initial_xin(P_tot,te_tot,n_tot,xi_n0)
     f_n=xi_n0
     f_p=1.0d0-xi_n0
     f_p_n=f_n/(f_n+2.0d0*f_p)
     f_p_p=2.0d0*f_p/(f_n+2.0d0*f_p)
  else
     te_tot=1.0d0
     n_tot=1.0d0
  endif
  
  p_h=f_p_n*p_tot
  p_m=f_p_p*p_tot
  ro_h=f_n*n_tot
  ro_m=f_p*n_tot
          

  
  b_theta=b_perp*sin(tmp)+b_tang*cos(tmp)
  b_phi  =b_perp*cos(tmp)-b_tang*sin(tmp)
  b_x=b_para*sin(theta_p)*cos(phi_p)+ &
       b_theta*cos(theta_p)*cos(phi_p)-b_phi*sin(phi_p)
  b_y=b_para*sin(theta_p)*sin(phi_p)+ &
       b_theta*cos(theta_p)*sin(phi_p)+b_phi*cos(phi_p)
  b_z=b_para*cos(theta_p)-b_theta*sin(theta_p) 

  !set perturvation
  nk=1
  ni=ix/2
  epsilon=1.0d-7
  delta=0.1d0
  do kk=1,nk;do ii=1,ni
     call random_number(phase)
     phase=phase*2.0*pi
     k_z=2.0d0*pi*kk/L_z
     k_x=2.0d0*pi*ii/L_x
     do k=1,kx;do j=1,jx;do i=1,ix
        vy_m(i,j,k)=vy_m(i,j,k)+epsilon*tanh(y(j)/thick)/cosh(y(j)/thick)* &
             sin(k_x*x(i)+k_z*z(k)+phase)
     enddo;enddo;enddo
  enddo;enddo
  do k=1,kx;do j=1,jx;do i=1,ix
     vx_m(i,j,k)=0.1d0/cosh(y(j)/(thick*0.5d0))*tanh(x(i)/(0.25*L_x))
  enddo;enddo;enddo
    

  !!!========================================================

  !convert PV2cq and set that value to global variable 'U_h' and/or 'U_m'
  call setcq(ro_m,vx_m,vy_m,vz_m,p_m,B_x,B_y,B_z, &
       ro_h,vx_h,vy_h,vz_h,p_h)
  !---------------------------------------------------------------------

  !set default output period and time duration--------------------------
  if(tend.lt.0.0) then
     tend=20.0d0
     dtout=tend/20.0
     if(flag_mpi.eq.0 .or.my_rank.eq.0)      print *,"TEND",dtout,tend
  endif
  !---------------------------------------------------------------------

end subroutine tearing
