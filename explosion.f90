subroutine explosion
  use parameters,only:pi
  use globalvar,only:ix,jx,kx,U_h,U_m,flag_bnd,beta,flag_b_stg,dtout,&
       flag_mhd,flag_mpi,my_rank,flag_pip,gm,beta,tend,&
       x,y,z,dx,dy,dz,n_fraction,ndim,debug_direction,margin,dxc,mpi_siz,s_order,mpi_pos
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
  double precision f_n,f_p,f_p_n,f_p_p,start(3),end(3)
  double precision radius,P_max,V_max,x0,y0,z0,theta_b,theta_p,phi_p,b0
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
  start(1)=-1.0d0 ;end(1)=1.0d0
  start(2)=-1.0d0 ;end(2)=1.0d0
  start(3)=-1.0d0 ;end(3)=1.0d0
  call set_coordinate(start,end)
!  call set_coordinate_NUG(mpi_siz(1)*(ix-2*margin(1))+2*margin(1),&
!       ix,x,dx,dxc,s_order,0.0d0,2.0d0,2.0d-2,1.0d0,1.0d3,mpi_pos(1),margin(1),1)

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
  !set parameter----------------------------------------  
  radius=0.3d0
  P_max=0.1d1
  V_max=0.5d0*0
  theta_b=pi*0.5*0
  !(x0,y0,z0) center of explosion
  orix=ix/2+1 ;oriy=jx/2+1;oriz=kx/2+1
  
  x0=x(orix)-0.5d0*dx(1)
  if(ix-margin(1)*2.eq.1)x(margin(1)+1)=x0  
  if (ndim .ge.2) then 
     y0=y(oriy)-0.5d0*dy(1)
  else
     y0=y(oriy)
  endif
  if (ndim.ge.3) then
     z0=z(oriz)-0.5d0*dz(1)
  else
     z0=z(oriz)
  endif
  !-------------------------------------------------

  do k=1,kx
     do j=1,jx
        do i=1,ix   
           P_h(i,j,k)=f_p_n*(1.0d0/gm+P_max* &
                exp(-((x(i))**2+(y(j)-y0)**2+(z(k)-z0)**2)/radius**2))
           P_m(i,j,k)=f_p_p*(1.0d0/gm+P_max* &
                exp(-((x(i))**2+(y(j)-y0)**2+(z(k)-z0)**2)/radius**2))
           vx_h(i,j,k)=V_max* &
                exp(-((x(i))**2+(y(j)-y0)**2+(z(k)-z0)**2)/radius**2)
           vx_m(i,j,k)=V_max* &
                exp(-((x(i))**2+(y(j)-y0)**2+(z(k)-z0)**2)/radius**2)
        enddo
     enddo
  enddo
  B0=sqrt(2.0/(gm*beta))
  ro_h=1.0d0*f_n
  ro_m=1.0d0*f_p
  vy_h=0.0d0
  vz_h=0.0d0
  vy_m=0.0d0
  vz_m=0.0d0
  if(debug_direction.eq.1) then
     phi_p=0.0d0
     theta_p=pi/2.0d0
  else if(debug_direction.eq.2) then
     phi_p=pi/2.0d0
     theta_p=pi/2.0d0     
  else if(debug_direction.eq.3) then
     phi_p=0.0d0
     theta_p=0.0d0
  else if(debug_direction.eq.4) then
     phi_p=pi/4.0d0
     theta_p=pi/2.0d0
  else if(debug_direction.eq.5) then
     phi_p=0.0d0
     theta_p=pi/4.0d0
  else if(debug_direction.eq.6) then
     phi_p=pi/2.0d0
     theta_p=pi/4.0d0
  else if(debug_direction.eq.7) then
     phi_p=pi/4.0d0
     theta_p=pi/4.0d0
  endif  
  b_para=B0
  b_theta=0
  b_phi=0
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

end subroutine explosion
