subroutine field_compression
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,&
       x,y,z,dx,dy,dz,n_fraction,ndim,debug_direction,margin,debug_parameter,T0
  use scheme_rot,only:pv2cq_mhd,pv2cq_hd
  use model_rot, only:set_coordinate,setcq
  implicit none
  double precision :: ro_h(1:ix,1:jx,1:kx),ro_m(1:ix,1:jx,1:kx)
  double precision :: vx_h(1:ix,1:jx,1:kx),vx_m(1:ix,1:jx,1:kx)
  double precision :: vy_h(1:ix,1:jx,1:kx),vy_m(1:ix,1:jx,1:kx)
  double precision :: vz_h(1:ix,1:jx,1:kx),vz_m(1:ix,1:jx,1:kx) 
  double precision :: P_h (1:ix,1:jx,1:kx),P_m (1:ix,1:jx,1:kx)
  double precision :: B_x (1:ix,1:jx,1:kx)
  double precision :: B_y (1:ix,1:jx,1:kx),mask(1:ix,1:jx,1:kx)
  double precision :: B_z (1:ix,1:jx,1:kx),p_tot(1:ix,1:jx,1:kx)
  double precision f_n,f_p,f_p_n,f_p_p,start(3),end(3),wtr,b0,ioneq,Te_0,v0


  !Find the equilibrium neutral fraction
  Te_0=T0/1.1604e4
  ioneq=(2.6e-19/dsqrt(Te_0))/(2.91e-14/(0.232+13.6/Te_0)*(13.6/Te_0)**0.39*dexp(-13.6/Te_0))
  f_n=ioneq/(ioneq+1.0d0)
  f_p=1.0d0-f_n
  f_p_n=f_n/(f_n+2.0d0*f_p)
  f_p_p=2.0d0*f_p/(f_n+2.0d0*f_p)

  if (my_rank.eq.0) print*,'Neutral fraction = ',f_n

!  !set ionization fraction-----------------
  if(flag_pip.eq.0) then
     f_n=1.0d0
     f_p=1.0d0
     f_p_n=1.0d0
     f_p_p=1.0d0
  endif
  
  !set ionization fraction-----------------
  !if(flag_pip.eq.0) then
  !   f_n=1.0d0
  !   f_p=1.0d0
  !   f_p_n=1.0d0
  !   f_p_p=1.0d0
  !else
  !   f_n=n_fraction
  !   f_p=1.0d0-n_fraction     
  !   f_p_n=f_n/(f_n+2.0d0*f_p)
  !   f_p_p=2.0d0*f_p/(f_n+2.0d0*f_p)
  !endif
  !----------------------------------------

  !Set coordinate (uniform grid)--------------------------
  !!set lower and upper coordinate
  start(1)=-10.0d0 ;end(1)=10.0d0
  start(2)=-1.0d0 ;end(2)=1.0d0
  start(3)=-1.0d0 ;end(3)=1.0d0
  call set_coordinate(start,end)
  !---------------------------------------
  
  !!default boundary condition----------------------
  if (flag_bnd(1) .eq.-1) flag_bnd(1)=11
  if (flag_bnd(2) .eq.-1) flag_bnd(2)=11 !in precedenza c'era boundary 2
  if (flag_bnd(3) .eq.-1) flag_bnd(3)=1
  if (flag_bnd(4) .eq.-1) flag_bnd(4)=1
  if (flag_bnd(5) .eq.-1) flag_bnd(5)=1
  if (flag_bnd(6) .eq.-1) flag_bnd(6)=1
  !-------------------------------------------------

  !!!========================================================
  !write some code to set physical variables
  wtr=0.1d0
  mask=spread(spread(x,2,jx),3,kx)
  ro_h=f_n*1.0d0
  ro_m=f_p*1.0d0
  
  B0=sqrt(2.0d0/gm/beta)
  v0=0.05
  
  b_x=0.0d0
  b_y=B0*(tanh(mask/wtr))
  b_z=0.0d0
  vx_h=v0*(tanh(mask/wtr));vy_h=0.0d0;vz_h=0.0d0
  vx_m=v0*(tanh(mask/wtr));vy_m=0.0d0;vz_m=0.0d0

  p_tot=1.0d0/gm
  p_h=f_p_n*p_tot
  p_m=f_p_p*p_tot
    
 

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

end subroutine field_compression
