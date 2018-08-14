subroutine Orszag_Tang
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,&
       x,y,z,dx,dy,dz,n_fraction,gra,flag_grav,scl_height,margin
  use scheme_rot,only:pv2cq_mhd,pv2cq_hd
  use model_rot, only:set_coordinate,setcq
  use matrix_rot,only:inverse_tridiagonal
  implicit none
  double precision :: ro_h(1:ix,1:jx,1:kx),ro_m(1:ix,1:jx,1:kx)
  double precision :: vx_h(1:ix,1:jx,1:kx),vx_m(1:ix,1:jx,1:kx)
  double precision :: vy_h(1:ix,1:jx,1:kx),vy_m(1:ix,1:jx,1:kx)
  double precision :: vz_h(1:ix,1:jx,1:kx),vz_m(1:ix,1:jx,1:kx) 
  double precision :: P_h (1:ix,1:jx,1:kx),P_m (1:ix,1:jx,1:kx)
  double precision :: B_x (1:ix,1:jx,1:kx)
  double precision :: B_y (1:ix,1:jx,1:kx)
  double precision :: B_z (1:ix,1:jx,1:kx)
  double precision f_n,f_p,f_p_n,f_p_p,start(3),end(3)
  double precision Atwood,ro_l,ro_u,vx_l,vx_u,w_lay,v0,b0
  double precision A(jx,3),b(jx),P_y(jx)
  integer j
  !set ionization fraction-----------------
  if(flag_pip.eq.0) then
     f_n=1.0d0
     f_p=1.0d0
     f_p_n=1.0d0
     f_p_p=1.0d0
  else
     f_n=n_fraction
     f_p=1.0d0-n_fraction     
     f_p_n=f_n/(f_n+2.0d0*f_p)
     f_p_p=2.0d0*f_p/(f_n+2.0d0*f_p)
  endif
  !----------------------------------------

  !Set coordinate (uniform grid)--------------------------
  !!set lower and upper coordinate
  start(1)=0.0d0 ;end(1)=1.0d0
  start(2)=0.0d0 ;end(2)=1.0d0
  start(3)=0.0d0 ;end(3)=1.0d0
  call set_coordinate(start,end)
  !---------------------------------------
  
  !!default boundary condition----------------------
  if (flag_bnd(1) .eq.-1) flag_bnd(1)=1
  if (flag_bnd(2) .eq.-1) flag_bnd(2)=1
  if (flag_bnd(3) .eq.-1) flag_bnd(3)=1
  if (flag_bnd(4) .eq.-1) flag_bnd(4)=1
  if (flag_bnd(5) .eq.-1) flag_bnd(5)=1
  if (flag_bnd(6) .eq.-1) flag_bnd(6)=1
  !-------------------------------------------------

  !!!========================================================
  !density of lower fluid is unity
  
  flag_grav=0
  v0=1.0d0
  b0=1.0/sqrt(4.0*pi)  
!  b0=1.0d0

  ro_h(:,:,:)=25.0d0/(36.0*pi)*f_n
  ro_m(:,:,:)=25.0d0/(36.0*pi)*f_p
  P_h(:,:,:)=5.0/(12.0*pi)*f_p_n
  P_m(:,:,:)=5.0/(12.0*pi)*f_p_p

  vx_h=v0*spread(spread(-sin(2.0*pi*y),1,ix),3,kx)
  vx_m=v0*spread(spread(-sin(2.0*pi*y),1,ix),3,kx)
  vy_h=v0*spread(spread(sin(2.0*pi*x),2,jx),3,kx)
  vy_m=v0*spread(spread(sin(2.0*pi*x),2,jx),3,kx)
  
  b_x=b0*spread(spread(-sin(2.0*pi*y),1,ix),3,kx)
  b_y=b0*spread(spread(sin(4.0*pi*x),2,jx),3,kx)
  b_z=0.0d0
  vz_h=0.0
  vz_m=0.0
 



  !!!========================================================

  !convert PV2cq and set that value to global variable 'U_h' and/or 'U_m'
  call setcq(ro_m,vx_m,vy_m,vz_m,p_m,B_x,B_y,B_z, &
       ro_h,vx_h,vy_h,vz_h,p_h)
  !---------------------------------------------------------------------

  !set default output period and time duration--------------------------
  if(tend.lt.0.0) then
     tend=2.5d0
     dtout=tend/20.0d0
     if(flag_mpi.eq.0 .or.my_rank.eq.0)      print *,"TEND",dtout,tend
  endif
  !---------------------------------------------------------------------

end subroutine Orszag_Tang
