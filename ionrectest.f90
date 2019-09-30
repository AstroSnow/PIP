subroutine ionrectest
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,col,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,&
       x,y,z,dx,dy,dz,n_fraction,debug_direction,debug_option,gra,flag_grav,total_prc,&
       margin,dxc,s_order,mpi_siz,mpi_pos, flag_IR, T0
  use scheme_rot,only:pv2cq_mhd,pv2cq_hd
  use model_rot, only:set_coordinate,setcq
  use procedures
  implicit none
  double precision :: ro_h(1:ix,1:jx,1:kx),ro_m(1:ix,1:jx,1:kx)
  double precision :: vx_h(1:ix,1:jx,1:kx),vx_m(1:ix,1:jx,1:kx)
  double precision :: vy_h(1:ix,1:jx,1:kx),vy_m(1:ix,1:jx,1:kx)
  double precision :: vz_h(1:ix,1:jx,1:kx),vz_m(1:ix,1:jx,1:kx) 
  double precision :: P_h (1:ix,1:jx,1:kx),P_m (1:ix,1:jx,1:kx)
  double precision :: B_x (1:ix,1:jx,1:kx)
  double precision :: B_y (1:ix,1:jx,1:kx)
  double precision :: B_z (1:ix,1:jx,1:kx)
  double precision :: temperature(1:ix+1),energy(1:ix+1),rho(1:ix+1)
  double precision ::mask(ix,jx,kx)
  double precision ::v_para(ix,jx,kx),v_perp(ix,jx,kx)
  double precision ::b_para(ix,jx,kx),b_perp(ix,jx,kx)
  double precision ::v_b_para(ix,jx,kx),b_b_para(ix,jx,kx)
  double precision ::v_theta(ix,jx,kx),v_phi(ix,jx,kx)
  double precision ::b_theta(ix,jx,kx),b_phi(ix,jx,kx)
  double precision f_n,f_p,f_p_n,f_p_p,start(3),end(3),B0,P_tot,P_top,P_base
  integer i,j,k


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(debug_direction.eq.1) then
! 0D problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!set ionization fraction-----------------
if ((flag_pip.eq.0) .OR. (flag_IR.eq.0)) then
	print*,'=== flag_PIP and flag_IR should be 1 ==='
	stop
else
	f_n=n_fraction
	f_p=1.0d0-n_fraction     
	f_p_n=f_n/(f_n+2.0d0*f_p)
	f_p_p=2.0d0*f_p/(f_n+2.0d0*f_p)
endif
!----------------------------------------

!Set coordinate (uniform grid)--------------------------
!!set lower and upper coordinate
start(1)=0.0d0 ;end(1)=1.0d0          !400.0d0/f_p
start(2)=-1.0d0 ;end(2)=1.0d0
start(3)=-1.0d0 ;end(3)=1.0d0
call set_coordinate(start,end)
!---------------------------------------

!initialise arrays
do k=1,kx;do j=1,jx;do i=1,ix
	ro_h(i,j,k)=1.0d0*f_n
	ro_m(i,j,k)=1.0d0*f_p
	p_h(i,j,k)=2.0d0*f_p_n/gm
	p_m(i,j,k)=2.0d0*f_p_p/gm
	vx_h(i,j,k)=0.0d0
	vx_m(i,j,k)=0.0d0
	vy_h(i,j,k)=0.0d0
	vy_m(i,j,k)=0.0d0
	vz_h(i,j,k)=0.0d0
	vz_m(i,j,k)=0.0d0
	b_x(i,j,k)=0.0d0
	b_y(i,j,k)=0.0d0
	b_z(i,j,k)=0.0d0
enddo;enddo;enddo

!!set boundary condition----------------------
flag_bnd(1)=1
flag_bnd(2)=1
flag_bnd(3)=1
flag_bnd(4)=1
flag_bnd(5)=1
flag_bnd(6)=1
!-------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
elseif(debug_direction.eq.2) then
! Leake 2013 reconnection paper TO DO!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

endif


!convert PV2cq and set that value to global variable 'U_h' and/or 'U_m'
call setcq(ro_m,vx_m,vy_m,vz_m,p_m,B_x,B_y,B_z, &
ro_h,vx_h,vy_h,vz_h,p_h)
!---------------------------------------------------------------------

!set default output period and time duration--------------------------
if(tend.lt.0.0) then
	tend=0.4d0
	dtout=tend/10.0
	if(flag_mpi.eq.0 .or.my_rank.eq.0)      print *,"TEND",dtout,tend
endif
!---------------------------------------------------------------------

end subroutine ionrectest
