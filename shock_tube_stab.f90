subroutine shock_tube_stab
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,col,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,&
       x,y,z,dx,dy,dz,n_fraction,debug_direction,debug_option
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
  double precision f_n,f_p,f_p_n,f_p_p,start(3),end(3),B0
  double precision theta_p,phi_p,tmp,v_L(8),v_R(8),wtr,wtrf
  double precision mach,rcom,rpres,alf,ang, byrat,vyrat,vu,bxu,byu,rou,pru,vxu
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
  start(1)=0.0d0/f_p ;end(1)=1500.0d0          !400.0d0/f_p
  start(2)=-1000.0d0 ;end(2)=1000.0d0
  start(3)=-1.0d0 ;end(3)=1.0d0
  call set_coordinate(start,end)
  !---------------------------------------
  
  !!default boundary condition----------------------
  if (flag_bnd(1) .eq.-1) flag_bnd(1)=10
  if (flag_bnd(2) .eq.-1) flag_bnd(2)=10
  if (flag_bnd(3) .eq.-1) flag_bnd(3)=1
  if (flag_bnd(4) .eq.-1) flag_bnd(4)=1
  if (flag_bnd(5) .eq.-1) flag_bnd(5)=10
  if (flag_bnd(6) .eq.-1) flag_bnd(6)=10
  !-------------------------------------------------

  !!!========================================================
  !parameters---------------------
!  B0=sqrt(2.0*gm/beta)
  B0=1.0d0
!/sqrt(2.0)
  !-------------------------
  if(debug_direction.eq.1) then
     mask=spread(spread(pi*x,2,jx),3,kx)
     phi_p=0.0d0
     theta_p=pi/2.0d0
     tmp=0.0
  else if(debug_direction.eq.2) then
     mask=spread(spread(pi*y,1,ix),3,kx)     
     phi_p=pi/2.0d0
     theta_p=pi/2.0d0     
     tmp=0.0
  else if(debug_direction.eq.3) then
     mask=spread(spread(pi*z,1,jx),1,ix)
     phi_p=0.0d0
     theta_p=0.0d0
     tmp=0.0
  else if(debug_direction.eq.4) then
     do k=1,kx;do j=1,jx; do i=1,ix
        if(x(i)+y(j).gt.0)then
           mask(i,j,k)=1
        else
           mask(i,j,k)=-1
        endif
     enddo;enddo;enddo
     phi_p=pi/4.0d0
     theta_p=pi/2.0d0
     tmp=pi/2.0d0
  else if(debug_direction.eq.5) then
     do k=1,kx;do j=1,jx; do i=1,ix
        if(x(i)+z(k).gt.0)then
           mask(i,j,k)=1
        else
           mask(i,j,k)=-1
        endif
     enddo;enddo;enddo
     phi_p=0.0d0
     theta_p=pi/4.0d0
     tmp=0.0
  else if(debug_direction.eq.6) then
     do k=1,kx;do j=1,jx; do i=1,ix
        if(y(j)+z(k).gt.0)then
           mask(i,j,k)=1
        else
           mask(i,j,k)=-1
        endif
     enddo;enddo;enddo
     phi_p=pi/2.0d0
     theta_p=pi/4.0d0
     tmp=0.0
  else if(debug_direction.eq.7) then
     do k=1,kx;do j=1,jx; do i=1,ix
        if(x(i)+y(j)+z(k).gt.0)then
           mask(i,j,k)=1
        else
           mask(i,j,k)=-1
        endif
     enddo;enddo;enddo
     phi_p=pi/4.0d0
     theta_p=pi/4.0d0
     tmp=0.0
  endif  
  select case(debug_option)
  case (2)  !Switch-off shock (Falle et al.1998)
     v_l=(/1.368d0,1.769d0,0.269d0,1.0d0,0.0d0,1.0d0,0.0d0,0.0d0/)
     v_r=(/1.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,1.0d0,0.0d0/)  
  case (3)  !parallel shock
	print*,'Parallel shock'
	mach=2.d0
	rcom=(gm+1.d0)*mach**2/(2.d0+(gm-1.d0)*mach**2)
	rpres=1.d0+gm*mach**2*(1.d0-1.d0/rcom)
!     v_l=(/rcom,rpres,10.0d0*sqrt(gm)/rcom+10.0d0*sqrt(gm),0.0d0,0.0d0,2.0d0*sqrt(gm)*mach,0.0d0,0.0d0/)
!     v_r=(/1.0d0,1.0d0,0.0d0,0.0d0,0.0d0,2.0d0*sqrt(gm)*mach,0.0d0,0.0d0/)  
     v_l=(/rcom,rpres/gm,-mach/rcom+mach,0.0d0,0.0d0,sqrt(2.d0/gm/beta),0.0d0,0.0d0/)
     v_r=(/1.0d0,1.0d0/gm,0.0d0,0.0d0,0.0d0,sqrt(2.d0/gm/beta),0.0d0,0.0d0/)  
!print*,'r',rcom
  case (4)  !switch-off shock NOT WORKING YET
	print*,'Switch-off shock'
!	mach=10.d0
	bxu=0.1d0
	byu=-1.d0!dsqrt(1.d0-bxu**2)
	rou=1.d0
!	pru=beta*(bxu**2+byu**2)/2.d0
	pru=beta*(byu**2)/2.d0
	ang=dtan(byu/bxu)
	alf=dsqrt(bxu**2+byu**2)/dsqrt(rou)
	vxu=-alf*dcos(ang)
	mach=vxu/dsqrt(gm*pru/rou)
	rcom=(gm+1.d0)*mach**2/(2.d0+(gm-1.d0)*mach**2)
	rpres=1.d0+gm*mach**2*(rcom-1.d0)/rcom*(1.d0- &
 (rcom*alf**2*((rcom+1.d0)*vxu**2-2.d0*rcom*alf**2*dcos(ang)**2))/&
(2.d0*(vxu**2-rcom*alf**2*dcos(ang)**2)**2))
     v_l=(/rcom*rou,rpres*pru,vxu/rcom,0.0d0,0.0d0,bxu,0.0d0,0.0d0/)
     v_r=(/rou,pru,vxu,0.0d0,0.0d0,bxu,byu,0.0d0/)
print*,'beta=',beta,2.d0*pru/byu**2
!print*,rcom 
  case (5)  !oblique shock
	mach=10.d0
	alf=1.0d0
	ang=0.5d0*3.1415d0
	rcom=(gm+1.d0)*mach**2/(2.d0+(gm-1.d0)*mach**2)
	vu=mach*sqrt(gm)
	byrat=rcom*((vu**2-alf**2.d0*dcos(ang)**2)/(vu**2-rcom*alf**2.d0*dcos(ang)**2))
	vyrat=(vu**2-alf**2.d0*dcos(ang)**2)/(vu**2-rcom*alf**2.d0*dcos(ang)**2)
	rpres=1.d0+gm*vu**2*(rcom-1.d0)/gm/rcom*(1.d0- &
 (rcom*alf**2*((rcom+1.d0)*rcom**2-2.d0*rcom*alf**2*dcos(ang)**2))/&
(2.d0*(vu**2-rcom*alf**2*dcos(ang)**2)))
     v_l=(/rcom,rpres,vu/rcom-mach,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0/)
     v_r=(/1.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0/)
  case (6)  !switch-off shock from simulation NOT WORKING YET
     v_l=(/1.9052d0,0.4256d0,0.0d0,-0.9551d0,0.0d0,0.3d0,0.0d0,0.0d0/)
     v_r=(/0.8372d0,0.3730d0,-1.8383d0,-0.0565d0,0.0d0,0.3d0,-0.822d0,0.0d0/)     
  case default
!     v_l=(/1.0d0,1.0d0,0.0d0,0.0d0,0.0d0,B0*0.75d0,B0,0.0d0/)
!     v_r=(/0.125d0,0.1d0,0.0d0,0.0d0,0.0d0,B0*0.75d0,-B0,0.0d0/)
     v_l=(/1.0d0,beta*B0**2/2.d0,0.0d0,0.0d0,0.0d0,B0*0.3d0,B0,0.0d0/)
     v_r=(/1.0d0,beta*B0**2/2.d0,0.0d0,0.0d0,0.0d0,B0*0.3d0,-B0,0.0d0/)
!     Density, pressure, velocity * 3, magnetic field *3


  end select

  where(mask<0)
     ro_h=f_n*v_l(1)
     p_h=f_p_n*v_l(2)
     ro_m=f_p*v_l(1)
     p_m=f_p_p*v_l(2)
     vx_h=v_l(3)
     vy_h=v_l(4)
     vz_h=v_l(5)
     vx_m=v_l(3)
     vy_m=v_l(4)
     vz_m=v_l(5)
     b_para=v_l(6)
     b_perp=v_l(7)
  elsewhere
     ro_h=f_n*v_r(1)
     p_h=f_p_n*v_r(2)
     ro_m=f_p*v_r(1)
     p_m=f_p_p*v_r(2)
     vx_h=v_r(3)
     vy_h=v_r(4)
     vz_h=v_r(5)
     vx_m=v_r(3)
     vy_m=v_r(4)
     vz_m=v_r(5)
     b_para=v_r(6)
     b_perp=v_r(7)
  end where
!  wtr=0.001d0
  wtr=1.0d0!dx(1)*1.d0
  wtrf=dx(1)*4.d0
!0.001d0
  do k=1,kx;do j=1,jx;do i=1,ix
     ro_h(i,j,k)=f_n*(v_l(1)+(v_r(1)-v_l(1))*(1.0d0+tanh(mask(i,j,k)/wtr-wtrf))*0.5d0)
     ro_m(i,j,k)=f_p*(v_l(1)+(v_r(1)-v_l(1))*(1.0d0+tanh(mask(i,j,k)/wtr-wtrf))*0.5d0)
     p_h(i,j,k)=f_p_n*(v_l(2)+(v_r(2)-v_l(2))*(1.0d0+tanh(mask(i,j,k)/wtr-wtrf))*0.5d0)
     p_m(i,j,k)=f_p_p*(v_l(2)+(v_r(2)-v_l(2))*(1.0d0+tanh(mask(i,j,k)/wtr-wtrf))*0.5d0)
     vx_h(i,j,k)=(v_l(3)+(v_r(3)-v_l(3))*(1.0d0+tanh(mask(i,j,k)/wtr-wtrf))*0.5d0)
     vx_m(i,j,k)=(v_l(3)+(v_r(3)-v_l(3))*(1.0d0+tanh(mask(i,j,k)/wtr-wtrf))*0.5d0)
     vy_h(i,j,k)=(v_l(4)+(v_r(4)-v_l(4))*(1.0d0+tanh(mask(i,j,k)/wtr-wtrf))*0.5d0)
     vy_m(i,j,k)=(v_l(4)+(v_r(4)-v_l(4))*(1.0d0+tanh(mask(i,j,k)/wtr-wtrf))*0.5d0)
     vz_h(i,j,k)=(v_l(5)+(v_r(5)-v_l(5))*(1.0d0+tanh(mask(i,j,k)/wtr-wtrf))*0.5d0)
     vz_m(i,j,k)=(v_l(5)+(v_r(5)-v_l(5))*(1.0d0+tanh(mask(i,j,k)/wtr-wtrf))*0.5d0)
     b_para(i,j,k)=(v_l(6)+(v_r(6)-v_l(6))*(1.0d0+tanh(mask(i,j,k)/wtr-wtrf))*0.5d0)
     b_perp(i,j,k)=(v_l(7)+(v_r(7)-v_l(7))*(1.0d0+tanh(mask(i,j,k)/wtr-wtrf))*0.5d0)

	if (x(i) .GE. 4000.d0) then	
!		ro_m(i,j,k)=ro_m(i,j,k)+0.2d0*dsin((x(i)-4000.d0)/100.0/3.14)*dcos((y(j))/20.d0/3.14d0)
	endif
  enddo;enddo;enddo


  b_theta=b_perp*sin(tmp)
  b_phi  =b_perp*cos(tmp)
  b_x=b_para*sin(theta_p)*cos(phi_p)+ &
       b_theta*cos(theta_p)*cos(phi_p)-b_phi*sin(phi_p)
  b_y=b_para*sin(theta_p)*sin(phi_p)+ &
       b_theta*cos(theta_p)*sin(phi_p)+b_phi*cos(phi_p)
  b_z=b_para*cos(theta_p)-b_theta*sin(theta_p) 
  !!!========================================================

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

end subroutine shock_tube_stab
