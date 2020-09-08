subroutine shock_tube_stab3
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
  double precision rorand1,rorand2
  integer i,j,k,dpl

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
  start(1)=-100.0d0 ;end(1)=100.d0!1200.0d0          !400.0d0/f_p
  start(2)=0.0d0 ;end(2)=10.0d0
  start(3)=-1.0d0 ;end(3)=1.0d0
  call set_coordinate(start,end)
  !---------------------------------------
  
  !!default boundary condition----------------------
  if (flag_bnd(1) .eq.-1) flag_bnd(1)=20
  if (flag_bnd(2) .eq.-1) flag_bnd(2)=20
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
     mask=spread(spread(pi*x,2,jx),3,kx)
     phi_p=0.0d0
     theta_p=pi/2.0d0
     tmp=0.0

	print*,'Parallel shock'
	mach=2.d0
	rcom=(gm+1.d0)*mach**2/(2.d0+(gm-1.d0)*mach**2)
	rpres=1.d0+gm*mach**2*(1.d0-1.d0/rcom)
!     v_l=(/rcom,rpres,10.0d0*sqrt(gm)/rcom+10.0d0*sqrt(gm),0.0d0,0.0d0,2.0d0*sqrt(gm)*mach,0.0d0,0.0d0/)
!     v_r=(/1.0d0,1.0d0,0.0d0,0.0d0,0.0d0,2.0d0*sqrt(gm)*mach,0.0d0,0.0d0/)  
     v_l=(/rcom,rpres/gm,-mach/rcom,0.0d0,0.0d0,dsqrt(2.d0/gm/beta),0.0d0,0.0d0/) !Use this one!
     v_r=(/1.0d0,1.0d0/gm,-mach,0.0d0,0.0d0,dsqrt(2.d0/gm/beta),0.0d0,0.0d0/)  
!v_l=v_r
!print*,'r',rpres

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
  wtr=0.0d0!0.001d0!dx(1)*1.d0
  wtrf=0.0d0!dx(1)*4.d0
!0.001d0
  do k=1,kx;do j=1,jx;do i=1,ix
     ro_h(i,j,k)=f_n*(v_l(1)+(v_r(1)-v_l(1))*(1.0d0+dtanh(mask(i,j,k)/wtr-wtrf))*0.5d0)
     ro_m(i,j,k)=f_p*(v_l(1)+(v_r(1)-v_l(1))*(1.0d0+dtanh(mask(i,j,k)/wtr-wtrf))*0.5d0)
     p_h(i,j,k)=f_p_n*(v_l(2)+(v_r(2)-v_l(2))*(1.0d0+dtanh(mask(i,j,k)/wtr-wtrf))*0.5d0)
     p_m(i,j,k)=f_p_p*(v_l(2)+(v_r(2)-v_l(2))*(1.0d0+dtanh(mask(i,j,k)/wtr-wtrf))*0.5d0)
     vx_h(i,j,k)=(v_l(3)+(v_r(3)-v_l(3))*(1.0d0+dtanh(mask(i,j,k)/wtr-wtrf))*0.5d0)
     vx_m(i,j,k)=(v_l(3)+(v_r(3)-v_l(3))*(1.0d0+dtanh(mask(i,j,k)/wtr-wtrf))*0.5d0)
     vy_h(i,j,k)=(v_l(4)+(v_r(4)-v_l(4))*(1.0d0+dtanh(mask(i,j,k)/wtr-wtrf))*0.5d0)
     vy_m(i,j,k)=(v_l(4)+(v_r(4)-v_l(4))*(1.0d0+dtanh(mask(i,j,k)/wtr-wtrf))*0.5d0)
     vz_h(i,j,k)=(v_l(5)+(v_r(5)-v_l(5))*(1.0d0+dtanh(mask(i,j,k)/wtr-wtrf))*0.5d0)
     vz_m(i,j,k)=(v_l(5)+(v_r(5)-v_l(5))*(1.0d0+dtanh(mask(i,j,k)/wtr-wtrf))*0.5d0)
     b_para(i,j,k)=(v_l(6)+(v_r(6)-v_l(6))*(1.0d0+dtanh(mask(i,j,k)/wtr-wtrf))*0.5d0)
     b_perp(i,j,k)=(v_l(7)+(v_r(7)-v_l(7))*(1.0d0+dtanh(mask(i,j,k)/wtr-wtrf))*0.5d0)

call srand(81728)
wtr=1.d0

	if ((x(i) .GE. 10.0d0) .AND. (x(i) .LE. 20.0d0)) then	
		do dpl=1,10
		rorand1=rand()
		rorand2=rand()
		ro_m(i,j,k)=ro_m(i,j,k)+f_p*0.1d0*rorand1*&
 dsin((x(i)-10.0d0)*3.14d0/10.0d0)* dcos((y(j)-rorand2*10.d0)*3.14d0/10.0d0*2.d0*wtr)
		ro_h(i,j,k)=ro_h(i,j,k)+f_n*0.1d0*rorand1*&
 dsin((x(i)-10.0d0)*3.14d0/10.0d0)* dcos((y(j)-rorand2*10.d0)*3.14d0/10.0d0*2.d0*wtr)
		wtr=wtr+1.d0
		enddo
	endif
  enddo;enddo;enddo
!print*,f_p_p,rpres/gm,p_m(1,1,1)
!print*,f_p,rcom,ro_m(1,1,1)

  b_theta=b_perp*dsin(tmp)
  b_phi  =b_perp*dcos(tmp)
  b_x=b_para*dsin(theta_p)*dcos(phi_p)+ &
       b_theta*dcos(theta_p)*dcos(phi_p)-b_phi*dsin(phi_p)
  b_y=b_para*dsin(theta_p)*dsin(phi_p)+ &
       b_theta*dcos(theta_p)*dsin(phi_p)+b_phi*dcos(phi_p)
  b_z=b_para*dcos(theta_p)-b_theta*dsin(theta_p) 
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

end subroutine shock_tube_stab3
