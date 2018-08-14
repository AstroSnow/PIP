subroutine Flare
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,&
       x,y,z,dx,dy,dz,n_fraction,ndim,debug_direction,margin,debug_parameter
  use scheme_rot,only:pv2cq_mhd,pv2cq_hd
  use model_rot, only:set_coordinate,setcq
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
  double precision B0,angle,Te_max,ro_chr,ro_cor,h_chr,w_tr

  integer i,j,k

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
  start(1)=0.0d0 ;end(1)=10.0d0
  start(2)=0.0d0 ;end(2)=debug_parameter
  start(3)=-1.0d0 ;end(3)=1.0d0
  call set_coordinate(start,end)
  !---------------------------------------
  
  !!default boundary condition----------------------
  if (flag_bnd(1) .eq.-1) flag_bnd(1)=3
  if (flag_bnd(2) .eq.-1) flag_bnd(2)=2
  if (flag_bnd(3) .eq.-1) flag_bnd(3)=11
  if (flag_bnd(4) .eq.-1) flag_bnd(4)=4
  if (flag_bnd(5) .eq.-1) flag_bnd(5)=1
  if (flag_bnd(6) .eq.-1) flag_bnd(6)=1
  !-------------------------------------------------

  !!!========================================================
  !write some code to set physical variables
  B0=sqrt(2.0d0/(gm*beta))
  angle=pi/4.0d0
  Te_max=10.0d0
  ro_chr=1.0d5
  ro_cor=1.0d0
  h_chr=1.0d0
  w_tr=0.2d0

!  B_x=B0
!  B_y=B0*sin(angle)
!  B_z=0.0d0
  do k=1,kx;do j=1,jx; do i=1,ix
     p_h(i,j,k)=1.0d0/gm*f_p_n
     p_m(i,j,k)=1.0d0/gm*f_p_p
     ro_h(i,j,k)=f_n*(ro_chr+(ro_cor-ro_chr)*(tanh((y(j)-h_chr)/w_tr)+1.0d0)*0.5d0)
     ro_m(i,j,k)=f_p*(ro_chr+(ro_cor-ro_chr)*(tanh((y(j)-h_chr)/w_tr)+1.0d0)*0.5d0)
     B_y(i,j,k)=B0*tanh(x(i)/0.5d0)
     B_z(i,j,k)=B0/cosh(x(i)/0.5d0)
  enddo;enddo;enddo
  vx_h=0.0d0;vx_m=0.0d0
  vy_h=0.0d0;vy_m=0.0d0
  vz_h=0.0d0;vz_m=0.0d0
  B_x=0.0d0
  

  !!!========================================================

  !convert PV2cq and set that value to global variable 'U_h' and/or 'U_m'
  call setcq(ro_m,vx_m,vy_m,vz_m,p_m,B_x,B_y,B_z, &
       ro_h,vx_h,vy_h,vz_h,p_h)
  !---------------------------------------------------------------------

  !set default output period and time duration--------------------------
  if(tend.lt.0.0) then
     tend=0.05d0
     dtout=tend/5.0
     if(flag_mpi.eq.0 .or.my_rank.eq.0)      print *,"TEND",dtout,tend
  endif
  !---------------------------------------------------------------------

end subroutine Flare
