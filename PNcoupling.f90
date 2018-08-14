subroutine PNcoupling
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,&
       x,y,z,dx,dy,dz,n_fraction,ndim,debug_direction,margin
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
  double precision ::mask(ix,jx,kx)
  double precision ::v_para(ix,jx,kx),v_perp(ix,jx,kx)
  double precision ::b_para(ix,jx,kx),b_perp(ix,jx,kx),b_tang(ix,jx,kx)
  double precision ::v_b_para(ix,jx,kx),b_b_para(ix,jx,kx)
  double precision ::v_theta(ix,jx,kx),v_phi(ix,jx,kx)
  double precision ::b_theta(ix,jx,kx),b_phi(ix,jx,kx)
  double precision f_n,f_p,f_p_n,f_p_p,start(3),end(3)
  double precision thick,b0,P_tot,Bz0,v0
  double precision theta_p,phi_p,tmp,T_n,T_p
  integer i,j,k,orix,oriy,oriz

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
  start(1)=-0.0d0 ;end(1)=2.0d0
  start(2)=0.0d0 ;end(2)=10.0d0
  start(3)=-1.0d0 ;end(3)=1.0d0
  call set_coordinate(start,end)
  !---------------------------------------
  
  !!default boundary condition----------------------
  if (flag_bnd(1) .eq.-1) flag_bnd(1)=10
  if (flag_bnd(2) .eq.-1) flag_bnd(2)=10
  if (flag_bnd(3) .eq.-1) flag_bnd(3)=10
  if (flag_bnd(4) .eq.-1) flag_bnd(4)=10
  if (flag_bnd(5) .eq.-1) flag_bnd(5)=10
  if (flag_bnd(6) .eq.-1) flag_bnd(6)=10
  !-------------------------------------------------

  !!!========================================================
  tmp=0.0d0
  !set debug direction----------------------------------------------
  if(debug_direction.eq.1) then
     mask=spread(spread(x,2,jx),3,kx)
     phi_p=0.0d0
     theta_p=pi/2.0d0
     flag_bnd(1)=1
     flag_bnd(2)=1     
     flag_bnd(1:2)=1
     flag_bnd(3)=3
  else if(debug_direction.eq.2) then
     mask=spread(spread(y,1,ix),3,kx)     
     phi_p=pi/2.0d0
     theta_p=pi/2.0d0     
     flag_bnd(1)=3
     flag_bnd(3)=3
     flag_bnd(4)=2
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
  v0=0.0d0
  
  ro_h=1.0d0*f_n
  ro_m=1.0d0*f_p
  vx_h=v0*0.5d0
  vy_h=0.0d0
  vz_h=0.0d0
  vx_m=-v0*0.5d0
  vy_m=0.0d0
  vz_m=0.0d0
  
!  p_h=1.0d0*f_p_n/gm
!  p_m=1.0d0*f_p_p/gm
  T_n=2.0d0
  T_p=0.5d0

  P_h=ro_h*T_n/gm
  P_m=2.0d0*ro_m*T_p/gm
  B_x=0.0d0
  B_y=0.0d0
  B_z=0.0d0

  !!!========================================================

  !convert PV2cq and set that value to global variable 'U_h' and/or 'U_m'
  call setcq(ro_m,vx_m,vy_m,vz_m,p_m,B_x,B_y,B_z, &
       ro_h,vx_h,vy_h,vz_h,p_h)
  !---------------------------------------------------------------------

  !set default output period and time duration--------------------------
  if(tend.lt.0.0) then
     tend=5.0d0
     dtout=tend/20.0
     if(flag_mpi.eq.0 .or.my_rank.eq.0)      print *,"TEND",dtout,tend
  endif
  !---------------------------------------------------------------------

end subroutine PNcoupling
