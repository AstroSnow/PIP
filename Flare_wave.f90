!! Non uniform mesh used
!! Gravitationally stratified atmosphere
!! A simple force free field
subroutine Flare_wave
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,&
       x,y,z,dx,dy,dz,n_fraction,ndim,debug_direction,margin,debug_parameter, &
       gra,flag_grav
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
  integer ierr,flag_err_NUG


  if(flag_grav.eq.0) then
     print*,'=== flag_grav should be 1 ==='
     call MPI_ABORT(mpi_comm_world,ierr)
  endif
  
  
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
  ! !!set lower and upper coordinate
  ! start(1)=0.0d0 ;end(1)=10.0d0
  ! start(2)=0.0d0 ;end(2)=debug_parameter
  ! start(3)=-1.0d0 ;end(3)=1.0d0
  ! call set_coordinate(start,end)
  !---------------------------------------

  !Set coordinate (Non-uniform grid)--------------------------
  !! X-direction
  x_uni0 = 0.d0
  x_uni1 = 4.d0
  dx_uni = 2.5d-2
  x_fact = 1.015d0
  dx_max = 5.d0*dx_uni
  call set_coordinate_NUG(mpi_siz(1)*(ix-2*margin(1))+2*margin(1), &
       ix,x,dx,dxc,s_order,x_uni0,x_uni1,dx_uni,x_fact,dx_max,mpi_pos(1),margin(1), &
       0,flag_err_NUG)
  if(flag_err_NUG.ne.0) then
     print*,'Grid # in x-direction is too small. STOP'
     call MPI_ABORT(mpi_comm_world,ierr)
  endif
  
  !! Y-direction
  y_uni0 = 0.d0
  y_uni1 = 24.d0
  dy_uni = 2.d-2
  y_fact = 1.015d0
  dy_max = 5.d0*dy_uni
  call set_coordinate_NUG(mpi_siz(2)*(jx-2*margin(2))+2*margin(2), &
       jx,y,dy,dyc,s_order,y_uni0,y_uni1,dy_uni,y_fact,dy_max,mpi_pos(2),margin(2), &
       0,flag_err_NUG)
  if(flag_err_NUG.ne.0) then
     print*,'Grid # in y-direction is too small. STOP'
     call MPI_ABORT(mpi_comm_world,ierr)
  endif
  
  !!default boundary condition----------------------
  if (flag_bnd(1) .eq.-1) flag_bnd(1)=3  ! sym
  if (flag_bnd(2) .eq.-1) flag_bnd(2)=2  ! sym
  if (flag_bnd(3) .eq.-1) flag_bnd(3)=11 ! free 
  if (flag_bnd(4) .eq.-1) flag_bnd(4)=10 ! free
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

  !! gravity
  flag_grav = 1
  gra(:,:,:,:) = 0.d0
  gra(:,:,:,2) = -1.d0/gm
  if(minval(y).lt.0.d0) then
     gra(:,1:margin(2),:,2) = 1.d0/gm
  endif

  !! gravitatioinally stratified atmosphere
  !! The density at the reconnection height is set to unity
  
  do k=1,kx;do j=1,jx; do i=1,ix
     ! p_h(i,j,k)=1.0d0/gm*f_p_n
     ! p_m(i,j,k)=1.0d0/gm*f_p_p
     ! ro_h(i,j,k)=f_n*(ro_chr+(ro_cor-ro_chr)*(tanh((y(j)-h_chr)/w_tr)+1.0d0)*0.5d0)
     ! ro_m(i,j,k)=f_p*(ro_chr+(ro_cor-ro_chr)*(tanh((y(j)-h_chr)/w_tr)+1.0d0)*0.5d0)
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

end subroutine Flare_Wave
