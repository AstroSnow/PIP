subroutine Prova_plasmoide
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
  integer i,j,k
  double precision wave_number,eps,B0,temp

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
  start(1)=-1.0d0 ;end(1)=1.0d0
  start(2)=0.0d0 ;end(2)=4.0d0
  start(3)=-1.0d0 ;end(3)=1.0d0
  call set_coordinate(start,end)
  !---------------------------------------
  
  !!default boundary condition----------------------
  if (flag_bnd(1) .eq.-1) flag_bnd(1)=1
  if (flag_bnd(2) .eq.-1) flag_bnd(2)=1
  if (flag_bnd(3) .eq.-1) flag_bnd(3)=4
  if (flag_bnd(4) .eq.-1) flag_bnd(4)=11
  if (flag_bnd(5) .eq.-1) flag_bnd(5)=1
  if (flag_bnd(6) .eq.-1) flag_bnd(6)=1
  !-------------------------------------------------
 
  !!!========================================================
  !write some code to set physical variables
  B0=sqrt(2.0d0/(gm*beta))
  wave_number=pi*3
  eps=debug_parameter
  ro_h=1.0d0*f_n
  ro_m=1.0d0*f_p
  vy_h=0.0d0
  vz_h=0.0d0
  vy_m=0.0d0
  vz_m=0.0d0 
  
  b_z=0.0d0
  do k=1,kx;do j=1,jx; do i=1,ix
     temp=cosh(wave_number*y(j))+eps*cos(wave_number*x(i))
     b_x(i,j,k)=-B0*sinh(wave_number*y(j))/temp          
     b_y(i,j,k)=-eps*B0*sin(wave_number*x(i))/temp
     p_h(i,j,k)=B0**2/2.0d0*(1.0d0-eps**2)/(temp**2)
     p_m(i,j,k)=B0**2/2.0d0*(1.0d0-eps**2)/(temp**2)
     vx_m(i,j,k)=-0.05*sin(wave_number*x(i)/2)*exp(-y(j)**2)*0
  enddo;enddo;enddo
  p_h=f_p_n*p_h
  p_m=f_p_p*p_m  

  !!!========================================================

  !convert PV2cq and set that value to global variable 'U_h' and/or 'U_m'
  call setcq(ro_m,vx_m,vy_m,vz_m,p_m,B_x,B_y,B_z, &
       ro_h,vx_h,vy_h,vz_h,p_h)
  !---------------------------------------------------------------------

  !set default output period and time duration--------------------------
  if(tend.lt.0.0) then
     tend=2.0d0
     dtout=tend/10.0
     if(flag_mpi.eq.0 .or.my_rank.eq.0)      print *,"TEND",dtout,tend
  endif
  !---------------------------------------------------------------------

end subroutine Prova_plasmoide
