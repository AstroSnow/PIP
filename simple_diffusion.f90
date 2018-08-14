subroutine simple_diffusion(n_fraction)
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
  double precision f_n,f_p,f_p_n,f_p_p,B0,theta_B,wtr,P_tot
  double precision l_x,k_x,l_y,k_y,amp

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
  l_x=4.0d0;
  l_y=10.0d0;

  xs=0.0d0 ;  xe=l_x
  ys=0.0d0 ;  ye=l_y
  zs=0.0d0 ;  ze=1.0d0
!!set grid space  
  dx0=(xe-xs)/(ix-2*margin(1)) ;dy0=(ye-ys)/(jx-2*margin(2));dz0=(ze-zs)/(kx-2*margin(3))
!!set origin
  orix=margin(1)+1; oriy=margin(2)+1 ;oriz=margin(3)+1
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
     if(flag_bnd(1).eq.0) flag_bnd(1)=2
     if(flag_bnd(2).eq.0) flag_bnd(2)=2
!     if(flag_bnd(3).eq.0) flag_bnd(3)=3
!     if(flag_bnd(4).eq.0) flag_bnd(4)=10
     if(flag_bnd(3).eq.0) flag_bnd(3)=1
     if(flag_bnd(4).eq.0) flag_bnd(4)=1
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
  radius=0.1d0
  Pmax=1.0d1
  wtr=0.5d0
  B0=sqrt(2.0d0/(beta*gm))

  print *,'n_fraction',n_fraction
  do k=1,kx
     do j=1,jx
        do i=1,ix   
           B_y(i,j,k)=B0*exp(-(x(i)/wtr)**2)
!           B_z(i,j,k)=B0/cosh(x(i)/wtr)
           P_tot=1.0d0/gm+0.5d0*(B0**2-B_y(i,j,k)**2)
!           P_tot=1.0d0/gm           
           P_h(i,j,k)=f_n*P_tot
           P_m(i,j,k)=f_p*P_tot
        enddo
     enddo
  enddo


  ro_h=1.0d0*f_n
  ro_m=1.0d0*f_p
  vx_h=0.0d0
  vy_h=0.0d0
  vz_h=0.0d0
  vx_m=0.0d0
  vy_m=0.0d0
  vz_m=0.0d0
!  B_x=sqrt(2.0d0/(beta*gm))*cos(theta_B)
  B_x=0.0d0
!  B_z=0.0d0


  !perturbation
  do k=1,kx
     do j=1,jx
        do i=1,ix
           vx_h(i,j,k)=-amp*x(i)*exp(-(x(i)/wtr)**2)*cos(k_y*y(j))
           vx_m(i,j,k)=-amp*x(i)*exp(-(x(i)/wtr)**2)*cos(k_y*y(j))
           vy_h(i,j,k)=amp*exp(-(x(i)/wtr)**2)*sin(k_y*y(j))
           vy_m(i,j,k)=amp*exp(-(x(i)/wtr)**2)*sin(k_y*y(j))
        enddo
     enddo
  enddo

  


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
  
end subroutine simple_diffusion
