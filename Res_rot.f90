module Res_rot
  use globalvar,only:ix,jx,kx,eta,nvar_h,nvar_m,eta_0,x,y,z,debug_parameter,&
       ndim,s_order,dxc,dyc,dzc,my_rank,j_cri,vd_cri
  use Util_rot,only:get_rotation
  implicit none  
  integer,save::resi_type
contains
  subroutine initialize_resistivity(flag_resi)
    integer,intent(inout)::flag_resi    
    if (flag_resi.eq.0) return
    allocate(eta(ix,jx,kx))        
    resi_type=mod((flag_resi/10),10)
    flag_resi=mod(flag_resi,10)
    if(my_rank.eq.0) print*,'resi_type, flag_resi: ',resi_type,flag_resi
  end subroutine initialize_resistivity
  subroutine set_resistivity(U_h,U_m)
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m) 
    double precision,parameter::r_x=0.8d0,r_y=0.8d0,r_z=5.0d0
!    double precision,parameter::J_cri=10.0d0
    double precision cur(ix,jx,kx,3),c2(ix,jx,kx),vdri(ix,jx,kx)
    double precision tmp
    integer i,j,k    

!    if(my_rank.eq.0) print*,'resi_type: ',resi_type
    select case(resi_type)
    case(0)!uniform resistivity 
       if(my_rank.eq.0) print*,'Uniform resistivity'
       eta(:,:,:)=eta_0
    case(1)!artificially localized resistivity
       do k=1,kx
          do j=1,jx
             do i=1,ix
                eta(i,j,k)=eta_0*exp(-((x(i)/r_x)**2+((y(j)-debug_parameter)/r_y)**2))
             enddo
          enddo
       enddo              
    case(2)!anomalous resistivity (Current density dependent)
       eta = 0.d0 ! initialization
       call get_rotation(ix,jx,kx,ndim,s_order,dxc,dyc,dzc,U_m(:,:,:,6:8),cur)
       c2=sqrt(cur(:,:,:,1)**2+cur(:,:,:,2)**2+cur(:,:,:,3)**2)
       do k=1,kx
          do j=1,jx
             do i=1,ix
                if(c2(i,j,k).ge.j_cri) eta(i,j,k)=eta_0*(c2(i,j,k)/j_cri-1.0d0)
             enddo
          enddo
       enddo                     
    case(3)!anomalous resistivity (drift velocity dependent)
!       if(my_rank.eq.0) print*,'anom resistivity'
       eta = 0.d0 ! initialization
       call get_rotation(ix,jx,kx,ndim,s_order,dxc,dyc,dzc,U_m(:,:,:,6:8),cur)
       vdri=sqrt(cur(:,:,:,1)**2+cur(:,:,:,2)**2+cur(:,:,:,3)**2)/U_m(:,:,:,1)
       do k=1,kx
          do j=1,jx
             do i=1,ix
                if(vdri(i,j,k).ge.vd_cri) eta(i,j,k)=eta_0*(vdri(i,j,k)/vd_cri-1.0d0)**2
             enddo
          enddo
       enddo                     
    end select
  end subroutine set_resistivity
end module Res_rot
