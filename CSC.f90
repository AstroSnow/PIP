subroutine CSC
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,&
       x,y,z,dx,dy,dz,dxc,dyc,dzc,col, &
       n_fraction,ndim,debug_direction,margin,xi_n,flag_amb,&
       flag_ir,mpi_siz,mpi_pos,s_order,debug_parameter
  use scheme_rot,only:pv2cq_mhd,pv2cq_hd
  use model_rot, only:set_coordinate,setcq
  use PIP_rot,only:set_xin,get_initial_xin,initialize_IR,&
       get_NT_from_PX
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
  double precision thick,b0,Bz0,inc_angle
  double precision theta_p,phi_p,tmp,dTe,r_x,wtr,dx0,dy0,factor
  integer i,j,k,orix,oriy,oriz,dummy



  !----------------------------------------
  !Set coordinate (uniform grid)--------------------------
  !!set lower and upper coordinate
  if(debug_direction.eq.2) then
     start(1)=-30.0d0 ;end(1)=30.0d0
     start(2)=0.0d0 ;end(2)=4.0d0
     start(3)=-1.0d0 ;end(3)=1.0d0
  else if(debug_direction.eq.1) then
     start(1)=0.0d0 ;end(1)=4.0d0
     start(2)=-1.0d0 ;end(2)=1.0d0
     start(3)=-1.0d0 ;end(3)=1.0d0     
  endif

  call set_coordinate(start,end)

!  call set_coordinate_NUG(mpi_siz(1)*(ix-2*margin(1))+2*margin(1),&
!       ix,x,dx,dxc,s_order,0.0d0,2.0d0,1.0d-2,1.015d0,1.0d3,mpi_pos(1),margin(1),0)
  

  !---------------------------------------  
!  if(flag_amb.eq.1.or.flag_pip.eq.1) then
!     call set_xin(U_h,U_m)
!  endif



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
  tmp=0.0d0
!  factor=1.015
  factor=1.015
  !set debug direction----------------------------------------------
  if(debug_direction.eq.1) then
     dx0=dx(1)
     mask=spread(spread(x,2,jx),3,kx)
     phi_p=0.0d0
     theta_p=pi/2.0d0
     flag_bnd(1)=3
     flag_bnd(2)=11
!     flag_bnd(2)=2
     flag_bnd(3)=4
     call set_coordinate_NUG(mpi_siz(1)*(ix-2*margin(1))+2*margin(1),&
          ix,x,dx,dxc,s_order,0.0d0,2.0d0,dx0,factor,1.0d3,mpi_pos(1),margin(1),0)
  else if(debug_direction.eq.2) then
     mask=spread(spread(y,1,ix),3,kx)     
     phi_p=pi/2.0d0
     theta_p=pi/2.0d0     
     flag_bnd(1)=11
     flag_bnd(2)=11
     flag_bnd(3)=4
     flag_bnd(4)=11
     dx0=dx(1)
     dy0=dy(1)
     call set_coordinate_NUG(mpi_siz(1)*(ix-2*margin(1))+2*margin(1),&
          ix,x,dx,dxc,s_order,0.0d0,2.0d1,dx0,factor, &
          2.0d-1,mpi_pos(1),margin(1),1)
     call set_coordinate_NUG(mpi_siz(2)*(jx-2*margin(2))+2*margin(2),&
          jx,y,dy,dyc,s_order,0.0d0,2.0d0,dy0,factor, &
          1.0d-1,mpi_pos(2),margin(2),0)

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
  thick=0.5d0

  B0=sqrt(2.0d0/(beta*gm))
  B0=1.0d0
  inc_angle=pi*debug_parameter
  Bz0=B0*tan(inc_angle)
  b_perp=B0*tanh(mask/thick)
  b_tang=Bz0  
!  p_tot=(1.0d0/gm+0.5d0*(b0*b0+b_tang**2-b_perp**2-b_tang**2))
  p_tot=0.5d0*B0*B0*(1.0d0+beta)-0.5d0*b_perp**2
  
  if(flag_pip.ge.1.or.flag_amb.ge.1) then
     te_tot=1.0d0
     r_x=10.0d0
!     dTe=0.35d0
     dTe=0.3d0
     wtr=2.0d0
     if(debug_direction.eq.2) then
        do k=1,kx;do j=1,jx;do i=1,ix
           te_tot(i,j,k)=te_tot(i,j,k)-dTe+dTe*((tanh((x(i)-r_x)/wtr)+1.0d0)*0.5d0 &
                +(1.0d0-tanh((x(i)+r_x)/wtr))*0.5d0)
        enddo;enddo;enddo     
     else
        te_tot=te_tot-dte
     endif
     dummy=1
!     n_tot=gm*p_tot/te_tot
     n_tot=1.0d0
     te_tot=gm*p_tot
     xi_n0(:,:,:)=n_fraction
!calculate total number density and neutral fraction from P_tot, T_tot     
!     call get_initial_xin(P_tot,te_tot,n_tot,xi_n0)     
!calculate total number density and temperature form P_tot,xi_n0
!     call get_NT_from_PX(P_tot,xi_n0,N_tot,Te_tot)
          
     f_n=xi_n0
     f_p=1.0d0-xi_n0
     f_p_n=f_n/(f_n+2.0d0*f_p)
     f_p_p=2.0d0*f_p/(f_n+2.0d0*f_p)
     xi_n=xi_n0
     if(flag_amb.ne.0) then
        f_p_p=f_p_n
        f_p=f_n
     endif
  else
     te_tot=1.0d0
     n_tot=gm*p_tot/te_tot
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


  !!!========================================================

  !convert PV2cq and set that value to global variable 'U_h' and/or 'U_m'
  call setcq(ro_m,vx_m,vy_m,vz_m,p_m,B_x,B_y,B_z, &
       ro_h,vx_h,vy_h,vz_h,p_h)
  !---------------------------------------------------------------------

  !set default output period and time duration--------------------------
  if(tend.lt.0.0) then
!     tend=2.0d0*col*n_fraction*(1.0d0-n_fraction)
     tend=2.0d0*col*minval(ro_h)*minval(ro_m)
     dtout=tend/40.0
     if(flag_mpi.eq.0 .or.my_rank.eq.0)      print *,"TEND",dtout,tend
  endif
  !---------------------------------------------------------------------

end subroutine CSC
