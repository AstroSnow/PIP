subroutine Alfven_damping
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,col,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,&
       x,y,z,dx,dy,dz,n_fraction,debug_option,debug_direction,flag_amb, &
       flag_resi,eta_0,debug_parameter

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
  double precision ::b_para(ix,jx,kx),b_perp(ix,jx,kx)
  double precision ::v_b_para(ix,jx,kx),b_b_para(ix,jx,kx)
  double precision ::v_theta(ix,jx,kx),v_phi(ix,jx,kx)
  double precision ::b_theta(ix,jx,kx),b_phi(ix,jx,kx),ro_t(ix,jx,kx)
  double precision f_n,f_p,f_p_n,f_p_p,start(3),end(3),omega_r,omega_i
  double precision amp,v_A,v_F,v_S,v_ph,B0,theta,theta_b,theta_p,phi_p,tmp,kn
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
  start(1)=-1.0d0 ;end(1)=9.0d0
  start(2)=-1.0d0 ;end(2)=1.0d0
  start(3)=-1.0d0 ;end(3)=1.0d0
  call set_coordinate(start,end)
  !---------------------------------------
  
  !!default boundary condition----------------------
  if (flag_bnd(1) .eq.-1) flag_bnd(1)=11
  if (flag_bnd(2) .eq.-1) flag_bnd(2)=11
  if (flag_bnd(3) .eq.-1) flag_bnd(3)=1
  if (flag_bnd(4) .eq.-1) flag_bnd(4)=1
  if (flag_bnd(5) .eq.-1) flag_bnd(5)=1
  if (flag_bnd(6) .eq.-1) flag_bnd(6)=1
  !-------------------------------------------------

  !!!========================================================
  !write some code to set physical variables
  !set parameters----
  amp=1.0d-8
  amp=1.0d-1
  kn=pi*(1.0d0+debug_parameter)
  B0=sqrt(2.0d0/(gm*beta))
  tmp=0.0
!  mask=spread(spread(kn*x,2,jx),3,kx)
  do i=1,ix
     if(x(i).le.0.0) then
        mask(i,:,:)=kn*(x(i)+0.5d0)
     else
        mask(i,:,:)=pi/2.0d0
     endif
  enddo
  phi_p=0.0d0
  theta_p=pi/2.0d0
  !------------------
  !Alfven mode (if flag_amb.eq.1 .or. flag_pip.eq.1 then
  !             Damping Alfven wave)
  v_A=B0
  if(flag_pip.eq.1.or.flag_amb.eq.1.or.flag_resi.eq.1) then 
     if (flag_resi.eq.1) then
        omega_i=0.5d0*kn**2*eta_0
     else
        omega_i=0.5d0*kn**2*(B0**2/(col*(1.0-n_fraction)/n_fraction))
     endif
     if (omega_i.le. v_A*kn) then 
        omega_r=sqrt((v_A*kn)**2-omega_i**2)
     else
        omega_r=0.0d0
        omega_i=omega_i+sqrt(omega_i**2-(v_A*kn)**2)
     endif
  else
     omega_r=v_A*kn
     omega_i=0.0d0
  endif
     
  v_ph=v_A

  do i=1,ix
     if(x(i).le.1.0d0) then
        ro_t(i,:,:)=1.0d0
     else if(x(i).le.8.0d0) then
        ro_t(i,:,:)=1.0d0-0.99*(x(i)-1.0d0)/7.0d0
     else
        ro_t(i,:,:)=0.01d0
     endif
  enddo

  ro_h=f_n*ro_t
  ro_m=f_p*ro_t
  p_h =f_p_n*(1.0d0/gm)
  p_m =f_p_p*(1.0d0/gm)
  vx_h =0.0
  vx_m =0.0
  vy_h =0.0d0
  vy_m =0.0d0
  
  B_x  =B0
  B_y  =0.0d0
  vz_h =amp*cos(mask)
  vz_m =amp*cos(mask)
  print *,"OMEGA Real and Imaginary : ",omega_r,omega_i
  B_z  =-b0/(v_A**2*kn)*amp*(omega_r*cos(mask)+omega_i*sin(mask))   
  if(flag_mpi.eq.0 .or.my_rank.eq.0)print *,"phase velocity=",v_ph
  !!!========================================================

  !convert PV2cq and set that value to global variable 'U_h' and/or 'U_m'
  call setcq(ro_m,vx_m,vy_m,vz_m,p_m,B_x,B_y,B_z, &
       ro_h,vx_h,vy_h,vz_h,p_h)
  !---------------------------------------------------------------------

  !set default output period and time duration--------------------------
  if(tend.lt.0.0) then
     tend=2.0
     if(debug_option.eq.3.and.omega_R.eq.0.0)tend=3.0d0/omega_I
     dtout=tend/10.0
     if(flag_mpi.eq.0 .or.my_rank.eq.0)      print *,"TEND",dtout,tend
  endif
  !---------------------------------------------------------------------


end subroutine Alfven_damping
