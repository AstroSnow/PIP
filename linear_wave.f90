subroutine linear_wave
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
  double precision ::b_theta(ix,jx,kx),b_phi(ix,jx,kx)
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
  start(1)=-1.0d0 ;end(1)=1.0d0
  start(2)=-1.0d0 ;end(2)=1.0d0
  start(3)=-1.0d0 ;end(3)=1.0d0
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
  !write some code to set physical variables
  !set parameters----
  amp=1.0d-8
!  amp=1.0d0
  kn=pi*(1.0d0+debug_parameter)
  B0=sqrt(2.0d0/(gm*beta))
  tmp=0.0
  if(debug_direction.eq.1) then
     mask=spread(spread(kn*x,2,jx),3,kx)
     phi_p=0.0d0
     theta_p=pi/2.0d0
  else if(debug_direction.eq.2) then
     mask=spread(spread(pi*y,1,ix),3,kx)     
     phi_p=pi/2.0d0
     theta_p=pi/2.0d0     
  else if(debug_direction.eq.3) then
     mask=spread(spread(pi*z,1,jx),1,ix)
     phi_p=0.0d0
     theta_p=0.0d0
  else if(debug_direction.eq.4) then
     do k=1,kx;do j=1,jx; do i=1,ix
        mask(i,j,k)=pi*(x(i)+y(j))
     enddo;enddo;enddo
     phi_p=pi/4.0d0
     theta_p=pi/2.0d0
     tmp=pi/2.0
  else if(debug_direction.eq.5) then
     do k=1,kx;do j=1,jx; do i=1,ix
        mask(i,j,k)=pi*(x(i)+z(k))
     enddo;enddo;enddo
     phi_p=0.0d0
     theta_p=pi/4.0d0
  else if(debug_direction.eq.6) then
     do k=1,kx;do j=1,jx; do i=1,ix
        mask(i,j,k)=pi*(y(j)+z(k))
     enddo;enddo;enddo
     phi_p=pi/2.0d0
     theta_p=pi/4.0d0
  else if(debug_direction.eq.7) then
     do k=1,kx;do j=1,jx; do i=1,ix
        mask(i,j,k)=pi*(x(i)+y(j)+z(k))
     enddo;enddo;enddo
     phi_p=pi/4.0d0
     theta_p=pi/4.0d0
  endif  
  !------------------
  if(flag_mhd.eq.0.or.debug_option.eq.0) then
     v_ph=1.0d0
     ro_h=f_n*(1.0d0+amp*sin(mask))
     ro_m=f_p*(1.0d0+amp*sin(mask))
     p_h =f_p_n*(1.0d0/gm+amp*sin(mask))
     p_m =f_p_p*(1.0d0/gm+amp*sin(mask))
     v_para=amp*sin(mask)
     v_perp=0.0d0
     vx_h=v_para*sin(theta_p)*cos(phi_p)
     vx_m=v_para*sin(theta_p)*cos(phi_p)
     vy_h=v_para*sin(theta_p)*sin(phi_p)
     vy_m=v_para*sin(theta_p)*sin(phi_p)
     vz_h=v_para*cos(theta_p)
     vz_m=v_para*cos(theta_p)
     B_x =0.0d0
     B_y =0.0d0
     B_z =0.0d0
  else if(debug_option.eq.1.or.debug_option.eq.2) then
     !Fast mode
     if(debug_option.eq.1) then
        !angle of magnetic field from propagation direction
        !wave speed
        theta=pi/2.0d0
        v_A=B0
        v_f=sqrt(0.5d0*((1.0d0+v_A**2)+ &
             sqrt((1.0d0+v_A**2)**2-4.0d0*v_A**2*cos(theta)**2)))
        v_s=sqrt(0.5d0*((1.0d0+v_A**2)- &
             sqrt((1.0d0+v_A**2)**2-4.0d0*v_A**2*cos(theta)**2)))     
        v_ph=v_f
        print*,'v_ph=v_f: ',v_ph
     else
        theta=0.0d0
        v_A=B0
        v_f=sqrt(0.5d0*((1.0d0+v_A**2)+ &
             sqrt((1.0d0+v_A**2)**2-4.0d0*v_A**2*cos(theta)**2)))
        v_s=sqrt(0.5d0*((1.0d0+v_A**2)- &
             sqrt((1.0d0+v_A**2)**2-4.0d0*v_A**2*cos(theta)**2)))     
        v_ph=v_s
        print*,'v_ph=v_s: ',v_ph        
     endif
     ro_h=f_n*(1.0d0+amp*sin(mask))
     ro_m=f_p*(1.0d0+amp*sin(mask))

     p_h =f_p_n*(1.0d0/gm+amp*sin(mask))
     p_m =f_p_p*(1.0d0/gm+amp*sin(mask))
     
!     v_b_para=1.0/v_ph*amp*sin(mask)*cos(theta)
!     v_perp=(v_ph**2-cos(theta)**2)/(v_ph*sin(theta))*amp*sin(mask)
!     b_b_para=b0*sin(theta)+b0*v_perp*sin(theta)/v_ph
!     b_perp  =b0*cos(theta)-b0*v_perp*cos(theta)/v_ph
!     v_para=v_b_para*cos(theta)+v_perp*sin(theta)
!     v_phi=-v_b_para*sin(theta)+v_perp*cos(theta)
!     v_theta=0.0
!     b_para=b_b_para*cos(theta)+b_perp*sin(theta)
!     b_phi=-b_b_para*sin(theta)+b_perp*cos(theta)
!     b_theta=0.0

     v_b_para=v_ph*amp*sin(mask)
     b_b_para=B0*cos(theta)
     v_perp=-(v_A*cos(theta))**2/(v_ph**2-(v_A*cos(theta))**2)*v_b_para
     b_perp=B0*(sin(theta)+v_ph/(v_ph**2-(v_A*cos(theta))**2)*v_b_para)   
     v_para=v_b_para
     b_para=b_b_para
!     select case(debug_direction)
!     case(5)        
!        v_theta=v_perp*cos(theta_p)
!        v_phi  =-v_perp*sin(theta_p)*sin(phi_p) 
!        b_theta=b_perp*cos(theta_p)
!        b_phi  =-b_perp*sin(theta_p)*sin(phi_p) 
!     case default
     v_theta=v_perp*sin(tmp)     
     v_phi  =v_perp*cos(tmp)     
     b_theta=b_perp*sin(tmp)     
     b_phi  =b_perp*cos(tmp)
!     end select

     vx_h=v_para*sin(theta_p)*cos(phi_p)+ &
          v_theta*cos(theta_p)*cos(phi_p)-v_phi*sin(phi_p)
     vx_m=v_para*sin(theta_p)*cos(phi_p)+ &
          v_theta*cos(theta_p)*cos(phi_p)-v_phi*sin(phi_p)
     vy_h=v_para*sin(theta_p)*sin(phi_p)+ &
          v_theta*cos(theta_p)*sin(phi_p)+v_phi*cos(phi_p)
     vy_m=v_para*sin(theta_p)*sin(phi_p)+ &
          v_theta*cos(theta_p)*sin(phi_p)+v_phi*cos(phi_p)
     vz_h=v_para*cos(theta_p)-v_theta*sin(theta_p)
     vz_m=v_para*cos(theta_p)-v_theta*sin(theta_p)

     b_x=b_para*sin(theta_p)*cos(phi_p)+ &
          b_theta*cos(theta_p)*cos(phi_p)-b_phi*sin(phi_p)
     b_y=b_para*sin(theta_p)*sin(phi_p)+ &
          b_theta*cos(theta_p)*sin(phi_p)+b_phi*cos(phi_p)
     b_z=b_para*cos(theta_p)-b_theta*sin(theta_p)

  else if(debug_option.eq.3) then
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
     ro_h=f_n
     ro_m=f_p
     p_h =f_p_n*(1.0d0/gm)
     p_m =f_p_p*(1.0d0/gm)
     print*,'v_ph=v_A: ',v_ph

     select case(debug_direction)
     case(1)
        B_x  =B0
        B_y  =0.0d0
        vx_h =0.0
        vx_m =0.0
        vy_h =0.0d0
        vy_m =0.0d0
        vz_h =0.0d0
        vz_m =amp*cos(mask)
        print *,"OMEGA Real and Imaginary : ",omega_r,omega_i
        B_z  =-b0/(v_A**2*kn)*amp*(omega_r*cos(mask)+omega_i*sin(mask))   
     case(2)
        B_y  =B0
        B_z  =0.0d0
        vy_h =0.0d0
        vy_m =0.0d0
        vz_h =0.0d0
        vz_m =0.0
        vx_h =0.0
        vx_m =amp*cos(mask)
        print *,"OMEGA Real and Imaginary : ",omega_r,omega_i
        B_x  =-b0/(v_A**2*kn)*amp*(omega_r*cos(mask)+omega_i*sin(mask))   
     case(3)
        B_z  =B0
        B_x  =0.0d0
        vz_h =0.0
        vz_m =0.0
        vx_h =0.0d0
        vx_m =0.0d0
        vy_h =0.0d0
        vy_m =amp*cos(mask)
        print *,"OMEGA Real and Imaginary : ",omega_r,omega_i
        B_y  =-b0/(v_A**2*kn)*amp*(omega_r*cos(mask)+omega_i*sin(mask))   
     end select
  endif
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


end subroutine linear_wave
