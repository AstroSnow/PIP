subroutine HC_tests
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,&
       x,y,z,dx,dy,dz,n_fraction,gra,flag_grav,scl_height,margin,&
       debug_direction
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
  double precision T0(1:ix)
  double precision, allocatable:: p_tot(:),ro_tot(:),temp(:),b0(:)
  integer i,j,k
  integer, allocatable:: new(:), old(:)
  integer size, seed(2), gseed(2), hiseed(2), zseed(2)


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

if(debug_direction .eq. 1) then
	print*,'1D test of thermal conduction from Rempel 2017'
	!Set coordinate (uniform grid)--------------------------
	!!set lower and upper coordinate
	start(1)=0.0d0 ;end(1)=1.0d0
	start(2)=0.0d0 ;end(2)=1.0d0
	start(3)=0.0d0 ;end(3)=1.0d0
	call set_coordinate(start,end)
	!---------------------------------------

	!!default boundary condition----------------------
	if (flag_bnd(1) .eq.-1) flag_bnd(1)=10
	if (flag_bnd(2) .eq.-1) flag_bnd(2)=10
	if (flag_bnd(3) .eq.-1) flag_bnd(3)=1
	if (flag_bnd(4) .eq.-1) flag_bnd(4)=1
	if (flag_bnd(5) .eq.-1) flag_bnd(5)=1
	if (flag_bnd(6) .eq.-1) flag_bnd(6)=1
	!-------------------------------------------------

	!Initial temperature profile, see EQN 12 in https://iopscience.iop.org/article/10.3847/1538-4357/834/1/10/pdf
	T0(1:ix)=0.1+0.9*x**5
	
	!Set MHD properties
	ro_m=1.0
	do j=1,jx; do k=1,kx
		P_m(1:ix,j,k)=ro_m(1:ix,j,k)*T0/gm
	enddo;enddo
	vx_m=0.0
	vy_m=0.0
	vz_m=0.0
	B_x=1.0
	B_y=0.0
	B_z=0.0


	!Set HD values to zero but should be used anyway!	
	ro_h=0.0
	P_h=0.0
	vx_h=0.0
	vy_h=0.0
	vz_h=0.0

endif

  !!!========================================================
  !convert PV2cq and set that value to global variable 'U_h' and/or 'U_m'
  call setcq(ro_m,vx_m,vy_m,vz_m,p_m,B_x,B_y,B_z, &
       ro_h,vx_h,vy_h,vz_h,p_h)
  !---------------------------------------------------------------------

  !set default output period and time duration--------------------------
  if(tend.lt.0.0) then
     tend=4.0d0
     dtout=tend/4.0/5.d0
     if(flag_mpi.eq.0 .or.my_rank.eq.0)      print *,"TEND",dtout,tend
  endif
  !---------------------------------------------------------------------

end subroutine HC_tests
