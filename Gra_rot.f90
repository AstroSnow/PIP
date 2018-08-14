module Gra_rot
  use globalvar,only:ndim,ix,jx,kx,flag_mhd,flag_pip,gra,nvar_h,nvar_m,gm
  implicit none

  double precision,save::scl_height
contains
  subroutine initialize_gravity(flag_gra)
    integer,intent(inout)::flag_gra
    if (flag_gra.eq.0) return
    if(scl_height.eq.0.0) scl_height=1.0    
    allocate(gra(ix,jx,kx,3))        
  end subroutine initialize_gravity  
  
  subroutine set_gravity(U_h,U_m)
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)
!    gra(:,:,:,:)=0.0d0
!    gra(:,:,:,ndim)=-1.0/(gm*scl_height)
  end subroutine set_gravity

  subroutine source_gravity(S_h,S_m,U_h,U_m)
    double precision,intent(inout)::S_h(ix,jx,kx,nvar_h),S_m(ix,jx,kx,nvar_m)
    double precision,intent(in)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)
    integer n
    if(flag_mhd.eq.0.or.flag_pip.eq.1) then
       do n=2,4
          S_h(:,:,:,n)=S_h(:,:,:,n)+U_h(:,:,:,1)*gra(:,:,:,n-1)
       enddo
       S_h(:,:,:,5) = S_h(:,:,:,5) &
            + U_h(:,:,:,2)*gra(:,:,:,1) &
            + U_h(:,:,:,3)*gra(:,:,:,2) &
            + U_h(:,:,:,4)*gra(:,:,:,3)
    endif
    if(flag_mhd.eq.1) then
       do n=2,4
          S_m(:,:,:,n)=S_m(:,:,:,n)+U_m(:,:,:,1)*gra(:,:,:,n-1)
       enddo
       S_m(:,:,:,5) = S_m(:,:,:,5) &
            + U_m(:,:,:,2)*gra(:,:,:,1) &
            + U_m(:,:,:,3)*gra(:,:,:,2) &
            + U_m(:,:,:,4)*gra(:,:,:,3)
    endif
  end subroutine source_gravity

end module Gra_rot
