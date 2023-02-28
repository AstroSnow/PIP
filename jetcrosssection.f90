subroutine jetcrosssection
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,&
       x,y,z,dx,dy,dz,dxc,dyc,dzc, &
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
     start(1)=-1.0d0 ;end(1)=1.0d0
     start(2)=-1.0d0 ;end(2)=1.0d0
     start(3)=-1.0d0 ;end(3)=1.0d0

  call set_coordinate(start,end)

  
  !!default boundary condition----------------------
  if (flag_bnd(1) .eq.-1) flag_bnd(1)=10
  if (flag_bnd(2) .eq.-1) flag_bnd(2)=10
  if (flag_bnd(3) .eq.-1) flag_bnd(3)=10
  if (flag_bnd(4) .eq.-1) flag_bnd(4)=10
  if (flag_bnd(5) .eq.-1) flag_bnd(5)=10
  if (flag_bnd(6) .eq.-1) flag_bnd(6)=10
  !-------------------------------------------------

  !Initialise arrays
  vx_h=0.0d0
  vy_h=0.0d0
  vz_h=0.0d0
  ro_h=0.d0
  p_h=0.d0
  vx_m=0.0d0
  vy_m=0.0d0
  vz_m=0.0d0
  ro_m=0.d0
  p_m=0.d0  
  b_x=0.d0!I am assuming no magnetic field is wanted for these simulations?
  b_y=0.d0
  b_z=0.d0

  do k=1,kx;do j=1,jx;do i=1,ix
  
  	vx_m(i,j,k)=0.14d0
  	vy_m(i,j,k)=0.14d0
  	vz_m(i,j,k)=0.9d0!You want this localised to the jet so replace with a formula
  	ro_m(i,j,k)=1.d0
  	p_m(i,j,k)=1.0/gm !Density and pressure choice are such that the sound speed is 1 in the normalised system
  
  
  enddo;enddo;enddo
          

  !!!========================================================

  !convert PV2cq and set that value to global variable 'U_h' and/or 'U_m'
  call setcq(ro_m,vx_m,vy_m,vz_m,p_m,B_x,B_y,B_z, &
       ro_h,vx_h,vy_h,vz_h,p_h)
  !---------------------------------------------------------------------

  !set default output period and time duration--------------------------
  if(tend.lt.0.0) then
     tend=10.0d0
     dtout=tend/20.0
     if(flag_mpi.eq.0 .or.my_rank.eq.0)      print *,"TEND",dtout,tend
  endif
  !---------------------------------------------------------------------

end subroutine jetcrosssection
