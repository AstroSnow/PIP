subroutine emerging_flux(n_fraction)
  use globalvar
  use scheme_routines
  implicit none
  double precision :: ro_h(1:ix,1:jx,1:kx),ro_m(1:ix,1:jx,1:kx)
  double precision :: vx_h(1:ix,1:jx,1:kx),vx_m(1:ix,1:jx,1:kx)
  double precision :: vy_h(1:ix,1:jx,1:kx),vy_m(1:ix,1:jx,1:kx)
  double precision :: vz_h(1:ix,1:jx,1:kx),vz_m(1:ix,1:jx,1:kx) 
  double precision :: P_h (1:ix,1:jx,1:kx),P_m (1:ix,1:jx,1:kx)
  double precision :: B_x (1:ix,1:jx,1:kx)
  double precision :: B_y (1:ix,1:jx,1:kx)
  double precision :: B_z (1:ix,1:jx,1:kx)
  double precision xs,xe,ys,ye,zs,ze,dx0,dy0,dz0,angle,ds,s,x0,y0,z0,radius
  double precision Pmax
  integer shock_cen,shock_hwid
  integer orix,oriy,oriz
  double precision,intent(inout)::n_fraction
  double precision f_n,f_p,f_p_n,f_p_p,B0,theta_B
  double precision p_tot(jx),n_tot(jx),xi_n0(jx),g0(jx),rbeta0(jx),te0(jx)
  double precision wtr,tcor,tadg,zcnv,ymin,af1,af2,zf1,zf2,ztr
  !Gravity limit around boundary
  double precision y_lower_bnd,y_upper_bnd,w_bnd,gym,wtr2
  double precision amp,xptb,yptb1,yptb2,wptb


  theta_B=0.5*pi

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


    
!Set coordinate
!!set lower and upper coordinate
  xs=0.0d0 ;  xe=40.0d0
  ys=-10.0d0 ;  ye=40.0d0
  zs=0.0d0 ;  ze=1.0d0
!!set grid space  
  dx0=(xe-xs)/(ix-2*margin(1)) ;dy0=(ye-ys)/(jx-2*margin(2));dz0=(ze-zs)/(kx-2*margin(3))

!!set origin
  orix=margin(1)+1
  oriy=int(-ys/dy0)+margin(2)+2
  oriz=margin(3)+1
  dx(:)=dx0
  x(orix)=0.5d0*dx0
  do i=orix+1,ix
     x(i)=x(i-1)+dx(i-1)
  enddo
  do i=orix-1,1,-1
     x(i)=x(i+1)-dx(i)
  enddo

  dy(:)=dy0
  y(oriy)=0.5d0*dy0
  do j=oriy+1,jx
     y(j)=y(j-1)+dy(j-1)
  enddo
  do j=oriy-1,1,-1
     y(j)=y(j+1)-dy(j)
  enddo

  dz(:)=dz0
  z(oriz)=0.5d0*dz(oriz)
  do k=oriz+1,kx
     z(k)=z(k-1)+dz(k-1)
  enddo
  do k=oriz-1,1,-1
     z(k)=z(k+1)-dz(k)
  enddo

!  do i=1,6
     if(flag_bnd(1).eq.0) flag_bnd(1)=3
     if(flag_bnd(2).eq.0) flag_bnd(2)=3
     if(flag_bnd(3).eq.0) flag_bnd(3)=2
     if(flag_bnd(4).eq.0) flag_bnd(4)=2
     if(flag_bnd(5).eq.0) flag_bnd(5)=1
     if(flag_bnd(6).eq.0) flag_bnd(6)=1
!  enddo     

  x0=x(orix)-0.5d0*dx0
  if (ndim .ge.2) then 
     y0=y(oriy)-0.5d0*dy0
  else
     y0=y(oriy)
  endif
  if (ndim.ge.3) then
     z0=z(oriz)-0.5d0*dz0
  else
     z0=z(oriz)
  endif

  wtr=0.5d0
  tcor=100.0
  ztr=10d0
  tadg=2.0d0
  zcnv=0.0d0

  y_lower_bnd=-8.0d0
  y_upper_bnd=35.0d0
  w_bnd=0.5d0



  gra(:,:,:,1)=0.0d0   
  gra(:,:,:,3)=0.0d0
  do j=1,jx
     af1=(tanh((y(j)-y_lower_bnd)/w_bnd)+1.0d0)/2.0d0
     af2=(-tanh((y(j)-y_upper_bnd)/w_bnd)+1.0d0)/2.0d0
     g0(j)=-af1*af2/gm
     gra(:,j,:,2)=g0(j)
  enddo

  do j=1,jx
     te0(j)=1.0d0 &
          -tadg*(gm-1.0d0)/gm*(y(j)-zcnv)*0.5*(1.0-tanh(y(j)-zcnv)) &
          +0.5d0*(tcor-1.0d0)*(tanh((y(j)-ztr)/wtr)+1.0d0)
  enddo
  
  zf1=-2.0
  zf2=0.0

  do j=1,jx
     af1=(tanh((y(j)-zf1)/0.5d0)+1.0d0)/2.0d0
     af2=(-tanh((y(j)-zf2)/0.5d0)+1.0d0)/2.0d0
     rbeta0(j)=1.0d0/beta*af1*af2     
  enddo

!Gravitational stratification
  n_tot(oriy)=1.0d0
  p_tot(oriy)=1.0d0/gm*te0(oriy)
  do j=oriy+1,jx 
     gym=0.5d0*(g0(j-1)+g0(j))
     n_tot(j)=n_tot(j-1) &
          *((1.0d0)*te0(j-1)+0.5d0*gm*gym*dy(j-1)) &
          /((1.0d0)*te0(j)-0.5d0*gm*gym*dy(j-1))
     p_tot(j)=p_tot(oriy)*(n_tot(j)/n_tot(oriy))*(te0(j)/te0(oriy))
  enddo
  
  do j=oriy-1,1,-1
     gym=0.5d0*(g0(j+1)+g0(j))
     n_tot(j)=n_tot(j+1) &
          *((1.0d0)*te0(j+1)-0.5d0*gm*gym*dy(j)) &
          /((1.0d0)*te0(j)    +0.5d0*gm*gym*dy(j))
     p_tot(j)=p_tot(oriy)*(n_tot(j)/n_tot(oriy))*(te0(j)/te0(oriy))
  enddo

  if(flag_pip.eq.1) then   
     call initial_atmos(te0,g0,xi_n0,p_tot,n_tot,jx,oriy,gm,dy)
     do k=1,kx
        do j=1,jx
           do i=1,ix   
              !           P_h(i,j,k)=f_n*P_tot(j)
              !           P_m(i,j,k)=f_p*P_tot(j)
              P_h(i,j,k)=xi_n0(j)/(2.0d0-xi_n0(j))*P_tot(j)           
              P_m(i,j,k)=(2.0d0-2.0d0*xi_n0(j))/(2.0d0-xi_n0(j))*P_tot(j)
              ro_h(i,j,k)=xi_n0(j)*n_tot(j)
              ro_m(i,j,k)=(1.0d0-xi_n0(j))*n_tot(j)           
              B_x(i,j,k)=sqrt(2.0d0*P_tot(j)*rbeta0(j))
           enddo
        enddo
     enddo
  else if(flag_mhd.eq.1) then
     print *,"EF_MHD"
     do k=1,kx
        do j=1,jx
           do i=1,ix   
              P_m(i,j,k)=P_tot(j)
              ! P_h(i,j,k)=f_n*P_tot(j)
              !           P_m(i,j,k)=f_p*P_tot(j)
!              P_h(i,j,k)=xi_n0(j)/(2.0d0-xi_n0(j))*P_tot(j)           
!              P_m(i,j,k)=(2.0d0-2.0d0*xi_n0(j))/(2.0d0-xi_n0(j))*P_tot(j)
!              ro_h(i,j,k)=xi_n0(j)*n_tot(j)
              ro_m(i,j,k)=n_tot(j)           
              B_x(i,j,k)=sqrt(2.0d0*P_tot(j)*rbeta0(j))
           enddo
        enddo
     enddo


     !perturbation 
     amp=0.05d0
     xptb=20.0d0
     yptb1=-2.0d0
     yptb2=0.0d0
     wptb=0.5d0

     do k=1,kx
        do j=1,jx
           do i=1,ix
              vy_m(i,j,k)=amp*cos(2*pi*x(i)/xptb) &
                   *0.5*(tanh((x(i)+0.75*xptb)/wptb)-tanh((x(i)-0.75*xptb)/wptb))&
                   *0.5*(tanh((y(j)-yptb1)/wptb)-tanh((y(j)-yptb2)/wptb))
           enddo
        enddo
     enddo
  endif

  

  vx_h=0.0d0
  vy_h=0.0d0
  vz_h=0.0d0
  vx_m=0.0d0
!  vy_m=0.0d0
  vz_m=0.0d0
!  B_x=sqrt(2.0d0/(beta*gm))*cos(theta_B)
  B_y=0.0d0
  B_z=0.0d0
  


  if(flag_divb.eq.2) then
     bx_cf(2:ix,:,:)=0.5d0*(B_x(1:ix-1,:,:)+B_x(1:ix-1,:,:))
     bx_cf(1,:,:)=bx_cf(2,:,:)
     bx_cf(ix+1,:,:)=bx_cf(ix,:,:)
     if(ndim.ge.2) then
        by_cf(:,2:jx,:)=0.5d0*(B_y(:,1:jx-1,:)+B_y(:,1:jx-1,:))
        by_cf(:,1,:)=by_cf(:,2,:)
        by_cf(:,jx+1,:)=by_cf(:,jx,:)     
        if(ndim.ge.3) then
           bz_cf(:,:,2:kx)=0.5d0*(B_z(:,:,1:kx-1)+B_z(:,:,1:kx-1))
           bz_cf(:,:,1)=bz_cf(:,:,2)
           bz_cf(:,:,kx+1)=bz_cf(:,:,kx)
        else
           bz_cf=B_z
        endif
     else
        by_cf=B_y
     endif     
  endif


  if(flag_mhd.eq.1) then 
     if(flag_divb.ne.2) then 
        call pv2cq_mhd(ro_m,vx_m,vy_m,vz_m,p_m,B_x,B_y,B_z)
     else
        call pv2cq_mhd(ro_m,vx_m,vy_m,vz_m,p_m,Bx_cf,By_cf,Bz_cf)
     endif
     if(flag_pip.eq.1) then
        call pv2cq_hd(ro_h,vx_h,vy_h,vz_h,p_h)
     endif
  else
     call pv2cq_hd(ro_h,vx_h,vy_h,vz_h,p_h)
  endif
  
end subroutine emerging_flux
