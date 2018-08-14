module HLL_rot

  !====================================================================
  ! This is the HLL-type routines module
  ! 
  ! First version - 2013/06/04 NN
  !====================================================================
  use Util_rot,only:minmod_func  
  use globalvar,only:ndim,ix,jx,kx,nvar_h,nvar_m,s_order,gm &
       ,ro_lim,pr_lim,tiny,flag_divb,cmax
  implicit none
  double precision,parameter::epsi=1.0d-10,acc=1.03d0
  !  use Util_rot,only:hllc_flux,hlld_flux,MinMod_func
contains
  subroutine get_Side(U,U_L,U_R,direction,nvar)
    integer,intent(in)::nvar
    double precision,intent(in)  :: U(ix,jx,kx,nvar)
    double precision,intent(out) :: U_L(ix,jx,kx,nvar)
    double precision,intent(out) :: U_R(ix,jx,kx,nvar)    
    integer,intent(in)::direction

    double precision dU(ix,jx,kx,nvar)
    integer ::dif(3)

    U_L = 0.d0 ; U_R = 0.d0

    dif(:)=0
    if(s_order.eq.1) then
       U_L=U
       dif(direction)=1
       U_R(1:ix-dif(1),1:jx-dif(2),1:kx-dif(3),:)= &
            U(1+dif(1):ix,1+dif(2):jx,1+dif(3):kx,:)
    else if(s_order.eq.2) then       
       if(direction.eq.1) then
          du(2:ix-1,:,:,:)= MinMod_func(U(1:ix-2,:,:,:),U(2:ix-1,:,:,:),U(3:ix,:,:,:), &
               ix-2,jx,kx,nvar)
          U_L(2:ix-1,:,:,:)=U(2:ix-1,:,:,:)+0.5*dU(2:ix-1,:,:,:)
          U_R(1:ix-2,:,:,:)=U(2:ix-1,:,:,:)-0.5*dU(2:ix-1,:,:,:)
       else if(direction.eq.2) then
          du(:,2:jx-1,:,:)= &
               MinMod_func(U(:,1:jx-2,:,:),U(:,2:jx-1,:,:),U(:,3:jx,:,:), &
               ix,jx-2,kx,nvar)
          U_L(:,2:jx-1,:,:)=U(:,2:jx-1,:,:)+0.5*dU(:,2:jx-1,:,:)
          U_R(:,1:jx-2,:,:)=U(:,2:jx-1,:,:)-0.5*dU(:,2:jx-1,:,:)
       else if(direction.eq.3) then
          du(:,:,2:kx-1,:)= &
               MinMod_func(U(:,:,1:kx-2,:),U(:,:,2:kx-1,:),U(:,:,3:kx,:), &
               ix,jx,kx-2,nvar)
          U_L(:,:,2:kx-1,:)=U(:,:,2:kx-1,:)+0.5*dU(:,:,2:kx-1,:)
          U_R(:,:,1:kx-2,:)=U(:,:,2:kx-1,:)-0.5*dU(:,:,2:kx-1,:)
       endif
    else
       U_L=U
       dif(direction)=1
       U_R(1:ix-dif(1),1:jx-dif(2),1:kx-dif(3),:)= &
            U(1+dif(1):ix,1+dif(2):jx,1+dif(3):kx,:)       
    endif    

    !! avoid ro=0
    U_L(:,:,:,1) = max(U_L(:,:,:,1),ro_lim)
    U_R(:,:,:,1) = max(U_R(:,:,:,1),ro_lim)    
    
  end subroutine get_Side

  subroutine hll_fluxes_ideal_hd(F_h,U_h)
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h)
    double precision,intent(out)::F_h(ix,jx,kx,nvar_h,3)
    integer direction0
    double precision V_h(ix,jx,kx,nvar_h)
    double precision V_L(ix,jx,kx,nvar_h),V_R(ix,jx,kx,nvar_h)    

    do direction0=1,ndim
       call u2v_hd(V_h,U_h)
       call get_Side(V_h,V_L,V_R,direction0,nvar_h)
       call set_HLL_fluxes_hd(F_h,V_L,V_R,direction0)       
    enddo
  end subroutine hll_fluxes_ideal_hd

  subroutine hll_fluxes_ideal_mhd(F_m,U_m)
    double precision,intent(inout)::U_m(ix,jx,kx,nvar_m)
    double precision,intent(out)::F_m(ix,jx,kx,nvar_m,3)
    integer direction0
    double precision V_m(ix,jx,kx,nvar_m)
    double precision V_L(ix,jx,kx,nvar_m),V_R(ix,jx,kx,nvar_m)    

    do direction0=1,ndim
       call u2v_mhd(V_m,U_m)
       call get_Side(V_m,V_L,V_R,direction0,nvar_m)
       call set_HLL_fluxes_mhd(F_m,V_L,V_R,direction0)       
    enddo
  end subroutine hll_fluxes_ideal_mhd

  
  subroutine hllc_fluxes_ideal(F_h,U_h)
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h)
    double precision,intent(out)::F_h(ix,jx,kx,nvar_h,3)
    double precision V_h(ix,jx,kx,nvar_h)
    double precision V_L(ix,jx,kx,nvar_h),V_R(ix,jx,kx,nvar_h)    
    integer direction0

    do direction0=1,ndim
       call u2v_hd(V_h,U_h)
       call get_Side(V_h,V_L,V_R,direction0,nvar_h)
       call set_HLLC_fluxes(F_h,V_L,V_R,direction0)       
    enddo    

  end subroutine hllc_fluxes_ideal

  subroutine hlld_fluxes_ideal(F_m,U_m)
    double precision,intent(inout)::U_m(ix,jx,kx,nvar_m)
    double precision,intent(out)::F_m(ix,jx,kx,nvar_m,3)
    integer direction0
    double precision V_m(ix,jx,kx,nvar_m)
    double precision V_L(ix,jx,kx,nvar_m),V_R(ix,jx,kx,nvar_m)  

    do direction0=1,ndim
       call u2v_mhd(V_m,U_m)
       call get_Side(V_m,V_L,V_R,direction0,nvar_m)
       call set_HLLD_fluxes(F_m,V_L,V_R,direction0) 
    enddo    
  end subroutine hlld_fluxes_ideal

  

  subroutine set_HLL_fluxes_hd(F_h,V_L,V_R,direction)
    integer,intent(in)::direction
    !    double precision,intent(out)::F_h(ix,jx,kx,nvar_h,3)
    double precision,intent(inout)::F_h(ix,jx,kx,nvar_h,3)    
    double precision,intent(in)::V_L(ix,jx,kx,nvar_h)
    double precision,intent(in)::V_R(ix,jx,kx,nvar_h)
    double precision flx_tmp(nvar_h)
    double precision v_l1(nvar_h)
    double precision v_r1(nvar_h)
    integer i,j,k

    if (direction.eq.1) then 
       do k=1,kx; do j=1,jx ;do i=1,ix-1
          v_l1=V_L(i,j,k,:)
          v_r1=V_R(i,j,k,:)
          call hll_flux_hd(v_l1,v_r1,flx_tmp,nvar_h,gm)
          F_h(i,j,k,:,1)=flx_tmp
       enddo;enddo;enddo
    else if(direction.eq.2)then
       do k=1,kx; do j=1,jx-1 ;do i=1,ix
          v_l1(1)=V_L(i,j,k,1);v_r1(1)=V_R(i,j,k,1)
          v_l1(2)=V_L(i,j,k,3);v_r1(2)=V_R(i,j,k,3)
          v_l1(3)=V_L(i,j,k,2);v_r1(3)=V_R(i,j,k,2)
          v_l1(4)=V_L(i,j,k,4);v_r1(4)=V_R(i,j,k,4)
          v_l1(5)=V_L(i,j,k,5);v_r1(5)=V_R(i,j,k,5)          
          call hll_flux_hd(v_l1,v_r1,flx_tmp,nvar_h,gm)
          F_h(i,j,k,1,2)=flx_tmp(1)
          F_h(i,j,k,2,2)=flx_tmp(3)
          F_h(i,j,k,3,2)=flx_tmp(2)
          F_h(i,j,k,4,2)=flx_tmp(4)
          F_h(i,j,k,5,2)=flx_tmp(5)
       enddo;enddo;enddo
    else 
       do k=1,kx-1; do j=1,jx ;do i=1,ix
          v_l1(1)=V_L(i,j,k,1);v_r1(1)=V_R(i,j,k,1)
          v_l1(2)=V_L(i,j,k,4);v_r1(2)=V_R(i,j,k,4)
          v_l1(3)=V_L(i,j,k,3);v_r1(3)=V_R(i,j,k,3)
          v_l1(4)=V_L(i,j,k,2);v_r1(4)=V_R(i,j,k,2)
          v_l1(5)=V_L(i,j,k,5);v_r1(5)=V_R(i,j,k,5)          
          call hll_flux_hd(v_l1,v_r1,flx_tmp,nvar_h,gm)
          F_h(i,j,k,1,3)=flx_tmp(1)
          F_h(i,j,k,2,3)=flx_tmp(4)
          F_h(i,j,k,3,3)=flx_tmp(3)
          F_h(i,j,k,4,3)=flx_tmp(2)
          F_h(i,j,k,5,3)=flx_tmp(5)
       enddo;enddo;enddo       
    endif

  end subroutine set_HLL_fluxes_hd


  subroutine set_HLL_fluxes_mhd(F,V_L,V_R,direction)
    integer,intent(in)::direction
    !    double precision,intent(out)::F(ix,jx,kx,nvar_m,3)
    double precision,intent(inout)::F(ix,jx,kx,nvar_m,3)    
    double precision,intent(in)::V_L(ix,jx,kx,nvar_m)
    double precision,intent(in)::V_R(ix,jx,kx,nvar_m)
    double precision flx_tmp(nvar_m)
    double precision v_l1(nvar_m)
    double precision v_r1(nvar_m)
    integer i,j,k

    select case(flag_divb)
    case(0) !! flag_divb=0
    if (direction.eq.1) then 
       do k=1,kx; do j=1,jx ;do i=1,ix-1
          v_l1=V_L(i,j,k,:)
          v_r1=V_R(i,j,k,:)
          call hll_flux_mhd(v_l1,v_r1,flx_tmp,nvar_m,gm)
          F(i,j,k,:,1)=flx_tmp
       enddo;enddo;enddo
    else if(direction.eq.2)then
       do k=1,kx; do j=1,jx-1 ;do i=1,ix
          v_l1(1)=V_L(i,j,k,1);v_r1(1)=V_R(i,j,k,1)
          v_l1(2)=V_L(i,j,k,3);v_r1(2)=V_R(i,j,k,3)
          v_l1(3)=V_L(i,j,k,2);v_r1(3)=V_R(i,j,k,2)
          v_l1(4)=V_L(i,j,k,4);v_r1(4)=V_R(i,j,k,4)
          v_l1(5)=V_L(i,j,k,5);v_r1(5)=V_R(i,j,k,5)          
          v_l1(6)=V_L(i,j,k,7);v_r1(6)=V_R(i,j,k,7)
          v_l1(7)=V_L(i,j,k,6);v_r1(7)=V_R(i,j,k,6)
          v_l1(8)=V_L(i,j,k,8);v_r1(8)=V_R(i,j,k,8)
          call hll_flux_mhd(v_l1,v_r1,flx_tmp,nvar_m,gm)
          F(i,j,k,1,2)=flx_tmp(1)
          F(i,j,k,2,2)=flx_tmp(3)
          F(i,j,k,3,2)=flx_tmp(2)
          F(i,j,k,4,2)=flx_tmp(4)
          F(i,j,k,5,2)=flx_tmp(5)
          F(i,j,k,6,2)=flx_tmp(7)
          F(i,j,k,7,2)=flx_tmp(6)
          F(i,j,k,8,2)=flx_tmp(8)
       enddo;enddo;enddo
    else 
       do k=1,kx-1; do j=1,jx ;do i=1,ix
          v_l1(1)=V_L(i,j,k,1);v_r1(1)=V_R(i,j,k,1)
          v_l1(2)=V_L(i,j,k,4);v_r1(2)=V_R(i,j,k,4)
          v_l1(3)=V_L(i,j,k,3);v_r1(3)=V_R(i,j,k,3)
          v_l1(4)=V_L(i,j,k,2);v_r1(4)=V_R(i,j,k,2)
          v_l1(5)=V_L(i,j,k,5);v_r1(5)=V_R(i,j,k,5)          
          v_l1(6)=V_L(i,j,k,8);v_r1(6)=V_R(i,j,k,8)
          v_l1(7)=V_L(i,j,k,7);v_r1(7)=V_R(i,j,k,7)
          v_l1(8)=V_L(i,j,k,6);v_r1(8)=V_R(i,j,k,6)          
          call hll_flux_mhd(v_l1,v_r1,flx_tmp,nvar_m,gm)
          F(i,j,k,1,3)=flx_tmp(1)
          F(i,j,k,2,3)=flx_tmp(4)
          F(i,j,k,3,3)=flx_tmp(3)
          F(i,j,k,4,3)=flx_tmp(2)
          F(i,j,k,5,3)=flx_tmp(5)
          F(i,j,k,6,3)=flx_tmp(8)
          F(i,j,k,7,3)=flx_tmp(7)
          F(i,j,k,8,3)=flx_tmp(6)
       enddo;enddo;enddo       
    endif
    case(1) !! flag_divb=1
    if (direction.eq.1) then 
       do k=1,kx; do j=1,jx ;do i=1,ix-1
          v_l1=V_L(i,j,k,:)
          v_r1=V_R(i,j,k,:)
          call hll_flux_mhd(v_l1,v_r1,flx_tmp,nvar_m,gm)
          F(i,j,k,:,1)=flx_tmp
       enddo;enddo;enddo
    else if(direction.eq.2)then
       do k=1,kx; do j=1,jx-1 ;do i=1,ix
          v_l1(1)=V_L(i,j,k,1);v_r1(1)=V_R(i,j,k,1)
          v_l1(2)=V_L(i,j,k,3);v_r1(2)=V_R(i,j,k,3)
          v_l1(3)=V_L(i,j,k,2);v_r1(3)=V_R(i,j,k,2)
          v_l1(4)=V_L(i,j,k,4);v_r1(4)=V_R(i,j,k,4)
          v_l1(5)=V_L(i,j,k,5);v_r1(5)=V_R(i,j,k,5)          
          v_l1(6)=V_L(i,j,k,7);v_r1(6)=V_R(i,j,k,7)
          v_l1(7)=V_L(i,j,k,6);v_r1(7)=V_R(i,j,k,6)
          v_l1(8)=V_L(i,j,k,8);v_r1(8)=V_R(i,j,k,8)
          v_l1(9)=V_L(i,j,k,9);v_r1(9)=V_R(i,j,k,9)          
          call hll_flux_mhd(v_l1,v_r1,flx_tmp,nvar_m,gm)
          F(i,j,k,1,2)=flx_tmp(1)
          F(i,j,k,2,2)=flx_tmp(3)
          F(i,j,k,3,2)=flx_tmp(2)
          F(i,j,k,4,2)=flx_tmp(4)
          F(i,j,k,5,2)=flx_tmp(5)
          F(i,j,k,6,2)=flx_tmp(7)
          F(i,j,k,7,2)=flx_tmp(6)
          F(i,j,k,8,2)=flx_tmp(8)
          F(i,j,k,9,2)=flx_tmp(9)          
       enddo;enddo;enddo
    else 
       do k=1,kx-1; do j=1,jx ;do i=1,ix
          v_l1(1)=V_L(i,j,k,1);v_r1(1)=V_R(i,j,k,1)
          v_l1(2)=V_L(i,j,k,4);v_r1(2)=V_R(i,j,k,4)
          v_l1(3)=V_L(i,j,k,3);v_r1(3)=V_R(i,j,k,3)
          v_l1(4)=V_L(i,j,k,2);v_r1(4)=V_R(i,j,k,2)
          v_l1(5)=V_L(i,j,k,5);v_r1(5)=V_R(i,j,k,5)          
          v_l1(6)=V_L(i,j,k,8);v_r1(6)=V_R(i,j,k,8)
          v_l1(7)=V_L(i,j,k,7);v_r1(7)=V_R(i,j,k,7)
          v_l1(8)=V_L(i,j,k,6);v_r1(8)=V_R(i,j,k,6)
          v_l1(9)=V_L(i,j,k,9);v_r1(9)=V_R(i,j,k,9)          
          call hll_flux_mhd(v_l1,v_r1,flx_tmp,nvar_m,gm)
          F(i,j,k,1,3)=flx_tmp(1)
          F(i,j,k,2,3)=flx_tmp(4)
          F(i,j,k,3,3)=flx_tmp(3)
          F(i,j,k,4,3)=flx_tmp(2)
          F(i,j,k,5,3)=flx_tmp(5)
          F(i,j,k,6,3)=flx_tmp(8)
          F(i,j,k,7,3)=flx_tmp(7)
          F(i,j,k,8,3)=flx_tmp(6)
          F(i,j,k,9,3)=flx_tmp(9)          
       enddo;enddo;enddo       
    endif
    end select

  end subroutine set_HLL_fluxes_mhd
  
  
  subroutine set_HLLC_fluxes(F_H,V_L,V_R,direction1)    
    !    double precision,intent(out)::F_H(ix,jx,kx,nvar_h,3)
    double precision,intent(inout)::F_H(ix,jx,kx,nvar_h,3)    
    double precision,intent(in)::V_L(ix,jx,kx,nvar_h)
    double precision,intent(in)::V_R(ix,jx,kx,nvar_h)
    integer,intent(in) :: direction1
    double precision flx_tmp(nvar_h)
    double precision v_l1(nvar_h)
    double precision v_r1(nvar_h)
    integer i,j,k

    if (direction1.eq.1) then 
       do k=1,kx; do j=1,jx ;do i=1,ix-1
          v_l1=V_L(i,j,k,:)
          v_r1=V_R(i,j,k,:)
          call hllc_flux(v_l1,v_r1,flx_tmp,nvar_h,gm)
          F_h(i,j,k,:,1)=flx_tmp
       enddo;enddo;enddo
    else if(direction1.eq.2)then
       do k=1,kx; do j=1,jx-1 ;do i=1,ix
          v_l1(1)=V_L(i,j,k,1);v_r1(1)=V_R(i,j,k,1)
          v_l1(2)=V_L(i,j,k,3);v_r1(2)=V_R(i,j,k,3)
          v_l1(3)=V_L(i,j,k,2);v_r1(3)=V_R(i,j,k,2)
          v_l1(4)=V_L(i,j,k,4);v_r1(4)=V_R(i,j,k,4)
          v_l1(5)=V_L(i,j,k,5);v_r1(5)=V_R(i,j,k,5)          
          call hllc_flux(v_l1,v_r1,flx_tmp,nvar_h,gm)
          F_h(i,j,k,1,2)=flx_tmp(1)
          F_h(i,j,k,2,2)=flx_tmp(3)
          F_h(i,j,k,3,2)=flx_tmp(2)
          F_h(i,j,k,4,2)=flx_tmp(4)
          F_h(i,j,k,5,2)=flx_tmp(5)
       enddo;enddo;enddo
    else 
       do k=1,kx-1; do j=1,jx ;do i=1,ix
          v_l1(1)=V_L(i,j,k,1);v_r1(1)=V_R(i,j,k,1)
          v_l1(2)=V_L(i,j,k,4);v_r1(2)=V_R(i,j,k,4)
          v_l1(3)=V_L(i,j,k,3);v_r1(3)=V_R(i,j,k,3)
          v_l1(4)=V_L(i,j,k,2);v_r1(4)=V_R(i,j,k,2)
          v_l1(5)=V_L(i,j,k,5);v_r1(5)=V_R(i,j,k,5)          
          call hllc_flux(v_l1,v_r1,flx_tmp,nvar_h,gm)
          F_h(i,j,k,1,3)=flx_tmp(1)
          F_h(i,j,k,2,3)=flx_tmp(4)
          F_h(i,j,k,3,3)=flx_tmp(3)
          F_h(i,j,k,4,3)=flx_tmp(2)
          F_h(i,j,k,5,3)=flx_tmp(5)
       enddo;enddo;enddo       
    endif

  end subroutine set_HLLC_fluxes

  subroutine set_HLLD_fluxes(F,V_L,V_R,direction)
    integer,intent(in)::direction
    !    double precision,intent(out)::F(ix,jx,kx,nvar_m,3)
    double precision,intent(inout)::F(ix,jx,kx,nvar_m,3)    
    double precision,intent(in)::V_L(ix,jx,kx,nvar_m)
    double precision,intent(in)::V_R(ix,jx,kx,nvar_m)
    double precision flx_tmp(nvar_m)
    double precision v_l1(nvar_m)
    double precision v_r1(nvar_m)
    integer i,j,k

    select case(flag_divb)
    case(0) !! flag_divb
    if (direction.eq.1) then 
       do k=1,kx; do j=1,jx ;do i=1,ix-1
          v_l1=V_L(i,j,k,:)
          v_r1=V_R(i,j,k,:)
          call hlld_flux(v_l1,v_r1,flx_tmp,nvar_m,gm)
          F(i,j,k,:,1)=flx_tmp
       enddo;enddo;enddo
    else if(direction.eq.2)then
       do k=1,kx; do j=1,jx-1 ;do i=1,ix
          v_l1(1)=V_L(i,j,k,1);v_r1(1)=V_R(i,j,k,1)
          v_l1(2)=V_L(i,j,k,3);v_r1(2)=V_R(i,j,k,3)
          v_l1(3)=V_L(i,j,k,2);v_r1(3)=V_R(i,j,k,2)
          v_l1(4)=V_L(i,j,k,4);v_r1(4)=V_R(i,j,k,4)
          v_l1(5)=V_L(i,j,k,5);v_r1(5)=V_R(i,j,k,5)          
          v_l1(6)=V_L(i,j,k,7);v_r1(6)=V_R(i,j,k,7)
          v_l1(7)=V_L(i,j,k,6);v_r1(7)=V_R(i,j,k,6)
          v_l1(8)=V_L(i,j,k,8);v_r1(8)=V_R(i,j,k,8)
          call hlld_flux(v_l1,v_r1,flx_tmp,nvar_m,gm)
          F(i,j,k,1,2)=flx_tmp(1)
          F(i,j,k,2,2)=flx_tmp(3)
          F(i,j,k,3,2)=flx_tmp(2)
          F(i,j,k,4,2)=flx_tmp(4)
          F(i,j,k,5,2)=flx_tmp(5)
          F(i,j,k,6,2)=flx_tmp(7)
          F(i,j,k,7,2)=flx_tmp(6)
          F(i,j,k,8,2)=flx_tmp(8)
       enddo;enddo;enddo
    else 
       do k=1,kx-1; do j=1,jx ;do i=1,ix
          v_l1(1)=V_L(i,j,k,1);v_r1(1)=V_R(i,j,k,1)
          v_l1(2)=V_L(i,j,k,4);v_r1(2)=V_R(i,j,k,4)
          v_l1(3)=V_L(i,j,k,3);v_r1(3)=V_R(i,j,k,3)
          v_l1(4)=V_L(i,j,k,2);v_r1(4)=V_R(i,j,k,2)
          v_l1(5)=V_L(i,j,k,5);v_r1(5)=V_R(i,j,k,5)          
          v_l1(6)=V_L(i,j,k,8);v_r1(6)=V_R(i,j,k,8)
          v_l1(7)=V_L(i,j,k,7);v_r1(7)=V_R(i,j,k,7)
          v_l1(8)=V_L(i,j,k,6);v_r1(8)=V_R(i,j,k,6)          
          call hlld_flux(v_l1,v_r1,flx_tmp,nvar_m,gm)
          F(i,j,k,1,3)=flx_tmp(1)
          F(i,j,k,2,3)=flx_tmp(4)
          F(i,j,k,3,3)=flx_tmp(3)
          F(i,j,k,4,3)=flx_tmp(2)
          F(i,j,k,5,3)=flx_tmp(5)
          F(i,j,k,6,3)=flx_tmp(8)
          F(i,j,k,7,3)=flx_tmp(7)
          F(i,j,k,8,3)=flx_tmp(6)
       enddo;enddo;enddo       
    endif
    case(1) !! flag_divb=1
    if (direction.eq.1) then 
       do k=1,kx; do j=1,jx ;do i=1,ix-1
          v_l1=V_L(i,j,k,:)
          v_r1=V_R(i,j,k,:)
          call hlld_flux(v_l1,v_r1,flx_tmp,nvar_m,gm)
          F(i,j,k,:,1)=flx_tmp
       enddo;enddo;enddo
    else if(direction.eq.2)then
       do k=1,kx; do j=1,jx-1 ;do i=1,ix
          v_l1(1)=V_L(i,j,k,1);v_r1(1)=V_R(i,j,k,1)
          v_l1(2)=V_L(i,j,k,3);v_r1(2)=V_R(i,j,k,3)
          v_l1(3)=V_L(i,j,k,2);v_r1(3)=V_R(i,j,k,2)
          v_l1(4)=V_L(i,j,k,4);v_r1(4)=V_R(i,j,k,4)
          v_l1(5)=V_L(i,j,k,5);v_r1(5)=V_R(i,j,k,5)          
          v_l1(6)=V_L(i,j,k,7);v_r1(6)=V_R(i,j,k,7)
          v_l1(7)=V_L(i,j,k,6);v_r1(7)=V_R(i,j,k,6)
          v_l1(8)=V_L(i,j,k,8);v_r1(8)=V_R(i,j,k,8)
          v_l1(9)=V_L(i,j,k,9);v_r1(9)=V_R(i,j,k,9)          
          call hlld_flux(v_l1,v_r1,flx_tmp,nvar_m,gm)
          F(i,j,k,1,2)=flx_tmp(1)
          F(i,j,k,2,2)=flx_tmp(3)
          F(i,j,k,3,2)=flx_tmp(2)
          F(i,j,k,4,2)=flx_tmp(4)
          F(i,j,k,5,2)=flx_tmp(5)
          F(i,j,k,6,2)=flx_tmp(7)
          F(i,j,k,7,2)=flx_tmp(6)
          F(i,j,k,8,2)=flx_tmp(8)
          F(i,j,k,9,2)=flx_tmp(9)          
       enddo;enddo;enddo
    else 
       do k=1,kx-1; do j=1,jx ;do i=1,ix
          v_l1(1)=V_L(i,j,k,1);v_r1(1)=V_R(i,j,k,1)
          v_l1(2)=V_L(i,j,k,4);v_r1(2)=V_R(i,j,k,4)
          v_l1(3)=V_L(i,j,k,3);v_r1(3)=V_R(i,j,k,3)
          v_l1(4)=V_L(i,j,k,2);v_r1(4)=V_R(i,j,k,2)
          v_l1(5)=V_L(i,j,k,5);v_r1(5)=V_R(i,j,k,5)          
          v_l1(6)=V_L(i,j,k,8);v_r1(6)=V_R(i,j,k,8)
          v_l1(7)=V_L(i,j,k,7);v_r1(7)=V_R(i,j,k,7)
          v_l1(8)=V_L(i,j,k,6);v_r1(8)=V_R(i,j,k,6)
          v_l1(9)=V_L(i,j,k,9);v_r1(9)=V_R(i,j,k,9)          
          call hlld_flux(v_l1,v_r1,flx_tmp,nvar_m,gm)
          F(i,j,k,1,3)=flx_tmp(1)
          F(i,j,k,2,3)=flx_tmp(4)
          F(i,j,k,3,3)=flx_tmp(3)
          F(i,j,k,4,3)=flx_tmp(2)
          F(i,j,k,5,3)=flx_tmp(5)
          F(i,j,k,6,3)=flx_tmp(8)
          F(i,j,k,7,3)=flx_tmp(7)
          F(i,j,k,8,3)=flx_tmp(6)
          F(i,j,k,9,3)=flx_tmp(9)          
       enddo;enddo;enddo       
    endif
    end select
       
  end subroutine set_HLLD_fluxes




  subroutine hll_flux_hd(v_l,v_r,flx_tmp,nvar,gm)
    implicit none
    integer,intent(in)::nvar    
    double precision,intent(inout)::v_l(nvar),v_r(nvar),gm      
    double precision,intent(out):: flx_tmp(nvar)
    double precision :: u_l(nvar),u_r(nvar),f_l(nvar),f_r(nvar),f_hll(nvar)
    double precision :: roiL,cs2L,csL
    double precision :: roiR,cs2R,csR    
    double precision :: SL,SR
    double precision :: s1,s2,s3
    integer :: n
    
    call v2u_1(v_l,u_l,nvar,0)
    call v2u_1(v_r,u_r,nvar,0)

    f_l(1)=u_l(2) ; f_r(1)=u_r(2)
    f_l(2)=v_l(5) +u_l(2)**2/u_l(1); f_r(2)=v_r(5) +u_r(2)**2/u_r(1)
    f_l(3)=u_l(2)*u_l(3)/u_l(1) ; f_r(3)=u_r(2)*u_r(3)/u_r(1)
    f_l(4)=u_l(2)*u_l(4)/u_l(1) ; f_r(4)=u_r(2)*u_r(4)/u_r(1)
    f_l(5)=u_l(2)*(u_l(5)+v_l(5))/u_l(1) ; f_r(5)=u_r(2)*(u_r(5)+v_r(5))/u_r(1)

    ! !! Using v_l,v_r
    ! f_l(1)=v_l(2) ; f_r(1)=v_r(2)
    ! f_l(2)=v_l(5) +v_l(1)*v_l(2)**2; f_r(2)=v_r(5) +v_r(1)*v_r(2)**2
    ! f_l(3)=v_l(1)*v_l(2)*v_l(3) ; f_r(3)=v_r(1)*v_r(2)*v_r(3)
    ! f_l(4)=v_l(1)*v_l(2)*v_l(4) ; f_r(4)=v_r(1)*v_r(2)*v_r(4)
    ! f_l(5)=v_l(1)*v_l(2)*(u_l(5)+v_l(5)) ; f_r(5)=v_r(1)*v_r(2)*(u_r(5)+v_r(5))
    
    
    !! Left
    roiL = 1.d0/max(v_l(1),tiny)

    cs2L = gm*v_l(5)*roiL
    csL = sqrt(max(cs2L,tiny))

    !! Right
    roiR = 1.d0/max(v_r(1),tiny)

    cs2R = gm*v_r(5)*roiR
    csR = sqrt(max(cs2R,tiny))

    !! Wave speed estimation
    SL = min(v_l(2),v_r(2)) - max(csL,csR)
    SR = max(v_l(2),v_r(2)) + max(csL,csR)

    !! HLL flux
    !    f_hll = ( SR*f_l - SL*f_r + SL*SR*(u_r-u_l))/(SR-SL)
    f_hll = ( SR*f_l - SL*f_r + SL*SR*(u_r-u_l))/max((SR-SL),tiny)

    !! switching function
    s1 = 0.5d0 + sign(0.5d0,SL)
    s3 = 0.5d0 - sign(0.5d0,SR)
    s2 = -sign(0.5d0,SL) + sign(0.5d0,SR)    

    !! Get flux
    do n=1,nvar
       flx_tmp(n) = s1*f_l(n) + s2*f_hll(n) + s3*f_r(n)
    enddo
    
  end subroutine hll_flux_hd



  subroutine hll_flux_mhd(v_l,v_r,flx_tmp,nvar,gm)
    implicit none
    integer,intent(in)::nvar
    double precision,intent(inout)::v_l(nvar),v_r(nvar),gm      
    double precision,intent(out):: flx_tmp(nvar)
    double precision :: u_l(nvar),u_r(nvar),f_l(nvar),f_r(nvar),f_hll(nvar)
    double precision :: roiL,cs2L,csL,bbL,ca2L,cb2L,cfast2L,cfastL
    double precision :: roiR,cs2R,csR,bbR,ca2R,cb2R,cfast2R,cfastR    
    double precision :: SL,SR
    double precision :: s1,s2,s3
    double precision :: ch
    integer :: n

    call v2u_1(v_l,u_l,nvar,1)
    call v2u_1(v_r,u_r,nvar,1)
    
    !! Left
    roiL = 1.d0/max(v_l(1),tiny)
    bbL  = v_l(6)**2 + v_l(7)**2 + v_l(8)**2

    cs2L = gm*v_l(5)*roiL
    ca2L = v_l(6)**2*roiL
    cb2L = cs2L + bbL*roiL
    cfast2L = 0.5d0*( cb2L + sqrt( abs(cb2L**2 - 4.d0*cs2L*ca2L) ) )
    cfastL = sqrt(max(cfast2L,tiny))

    !! Right
    roiR = 1.d0/max(v_r(1),tiny)
    bbR  = v_r(6)**2 + v_r(7)**2 + v_r(8)**2

    cs2R = gm*v_r(5)*roiR
    ca2R = v_r(6)**2*roiR
    cb2R = cs2R + bbR*roiR
    cfast2R = 0.5d0*( cb2R + sqrt( abs(cb2R**2 - 4.d0*cs2R*ca2R) ) )
    cfastR = sqrt(max(cfast2R,tiny))

    !! Wave speed estimation
    SL = min(v_l(2),v_r(2)) - max(cfastL,cfastR)
    SR = max(v_l(2),v_r(2)) + max(cfastL,cfastR)


    !! flux at L/R
    f_l(1)=u_l(2) 
    f_l(2)=v_l(5) +u_l(2)**2/u_l(1) + 0.5d0*bbL - u_l(6)**2
    f_l(3)=u_l(2)*u_l(3)/u_l(1) - u_l(7)*u_l(6)
    f_l(4)=u_l(2)*u_l(4)/u_l(1) - u_l(8)*u_l(6)
    f_l(5)=u_l(2)*(u_l(5)+v_l(5)+0.5d0*bbL)/u_l(1) &
         - u_l(6)*( u_l(6)*v_l(2)+u_l(7)*v_l(3)+u_l(8)*v_l(4) )
    f_l(6)=0.d0
    f_l(7) = v_l(2)*u_l(7) - v_l(3)*u_l(6)
    f_l(8) = v_l(2)*u_l(8) - v_l(4)*u_l(6)    
    
    f_r(1)=u_r(2)    
    f_r(2)=v_r(5) +u_r(2)**2/u_r(1) + 0.5d0*bbR - u_r(6)**2        
    f_r(3)=u_r(2)*u_r(3)/u_r(1) - u_r(7)*u_r(6)
    f_r(4)=u_r(2)*u_r(4)/u_r(1) - u_r(8)*u_r(6)
    f_r(5)=u_r(2)*(u_r(5)+v_l(5)+0.5d0*bbR)/u_r(1) &
         - u_r(6)*( u_r(6)*v_r(2)+u_r(7)*v_r(3)+u_r(8)*v_r(4) )
    f_r(6)=0.d0
    f_r(7) = v_r(2)*u_r(7) - v_r(3)*u_r(6)
    f_r(8) = v_r(2)*u_r(8) - v_r(4)*u_r(6)    
    
    !! HLL flux
    f_hll = ( SR*f_l - SL*f_r + SL*SR*(u_r-u_l))/(SR-SL)

    !! switching function
    s1 = 0.5d0 + sign(0.5d0,SL)
    s3 = 0.5d0 - sign(0.5d0,SR)
    s2 = -sign(0.5d0,SL) + sign(0.5d0,SR)    

    !! Get flux 
    do n=1,8
       flx_tmp(n) = s1*f_l(n) + s2*f_hll(n) + s3*f_r(n)
    enddo

    !! divB cleaning
    if(flag_divb.eq.1) then
       ch = cmax
       flx_tmp(6) = 0.5d0 * ( u_l(9)+u_r(9) - ch * (v_r(6)-v_l(6)) )
       flx_tmp(9) = 0.5d0 * ( ch*ch*(v_l(6)+v_r(6)) - ch * (v_r(9)-v_l(9)) )       
    endif
    
  end subroutine hll_flux_mhd

  
  !  subroutine hllc_flux(u_l,u_r,flx_tmp,nvar,gm)
  subroutine hllc_flux(v_l,v_r,flx_tmp,nvar,gm)  
    implicit none
    integer,intent(in)::nvar
    double precision,intent(in):: gm
    double precision,intent(inout)::v_l(nvar),v_r(nvar)
    double precision,intent(out):: flx_tmp(nvar)
    double precision u_l(nvar),u_r(nvar)
    double precision flx_l(nvar),flx_r(nvar)
    double precision pr_l,pr_r,vn_l,vn_r
    double precision cs_l,cs_r,sl,sm,sr,pr_m,st_l(5),st_r(5),thl,thm,thr

    call v2u_1(v_l,u_l,nvar,0)
    call v2u_1(v_r,u_r,nvar,0)

    ! call v2f_1(v_l,flx_l,nvar,0)
    ! call v2f_1(v_r,flx_r,nvar,0)    
    
    pr_l = v_l(5)
    pr_r = v_r(5)

    flx_l(1)=u_l(2) ; flx_r(1)=u_r(2)
    flx_l(2)=pr_l +u_l(2)**2/u_l(1); flx_r(2)=pr_r +u_r(2)**2/u_r(1)
    flx_l(3)=u_l(2)*u_l(3)/u_l(1) ; flx_r(3)=u_r(2)*u_r(3)/u_r(1)
    flx_l(4)=u_l(2)*u_l(4)/u_l(1) ; flx_r(4)=u_r(2)*u_r(4)/u_r(1)
    flx_l(5)=u_l(2)*(u_l(5)+pr_l)/u_l(1) ; flx_r(5)=u_r(2)*(u_r(5)+pr_r)/u_r(1)

    cs_l=sqrt(gm*pr_l/u_l(1))
    cs_r=sqrt(gm*pr_r/u_r(1))
    vn_l=u_l(2)/u_l(1)
    vn_r=u_r(2)/u_r(1)
    sl=acc*min(vn_l-cs_l,vn_r-cs_r)
    sr=acc*max(vn_l+cs_l,vn_r+cs_r)

    !    q(1)=u_l(1)*(sl-vn_l)
    !    r(1)=u_r(1)*(sr-vn_r)
    !    q(2)=u_l(1)*vn_l*(sl-vn_l)-pr_l
    !    r(2)=u_r(1)*vn_r*(sr-vn_r)-pr_r
    !    q(3)=u_l(3)*(sl-vn_l)-pr_l
    !    r(3)=u_r(3)*(sr-vn_r)-pr_r
    !    sm=(r(2)-q(2))/(r(1)-q(1))

    sm=(u_r(2)*(sr-vn_r)-u_l(2)*(sl-vn_l)  &
         +pr_l-pr_r)/(u_r(1)*(sr-vn_r)-u_l(1)*(sl-vn_l))

    !    pr_m1=sm*u_l(1)*(sl-vn_l)-u_l(2)*(sl-vn_l)+pr_l
    !    pr_m2=sm*u_r(1)*(sr-vn_r)-u_r(2)*(sr-vn_r)+pr_r
    !    pr_m=0.5d0*((pr_l+pr_r)+4.0d0*(vn_l-vn_r)/((u_l(1)+u_r(1))*(cs_l+cs_r)))
    pr_m=sm*u_l(1)*(sl-vn_l)-u_l(2)*(sl-vn_l)+pr_l
    !    print *,pr_l,pr_m1,pr_m2,pr_m,pr_r

    st_l(1)=u_l(1)*(sl-vn_l)/(sl-sm)
    st_r(1)=u_r(1)*(sr-vn_r)/(sr-sm)
    st_l(2)=st_l(1)*sm
    st_r(2)=st_r(1)*sm
    st_l(3)=st_l(1)*u_l(3)/u_l(1)
    st_r(3)=st_r(1)*u_r(3)/u_r(1)
    st_l(4)=st_l(1)*u_l(4)/u_l(1)
    st_r(4)=st_r(1)*u_r(4)/u_r(1)    
    !    st_l(5)=st_l(1)*(u_l(5)/u_l(1)+(sm-vn_l)*(sm+pr_l/u_l(1)*(sl-vn_l)))
    !    st_r(5)=st_r(1)*(u_r(5)/u_r(1)+(sm-vn_r)*(sm+pr_r/u_r(1)*(sr-vn_r)))
    st_l(5)=((sl-vn_l)*u_l(5)-pr_l*vn_l+pr_m*sm)/(sl-sm)
    st_r(5)=((sr-vn_r)*u_r(5)-pr_r*vn_r+pr_m*sm)/(sr-sm)


    thl=(abs(sl)+sl)/max(2.0d0*sl,epsi)
    thm=(abs(sm)+sm)/max(2.0d0*sm,epsi)
    thr=1.0d0-(abs(sr)+sr)/max(2.0d0*sr,epsi)


    flx_tmp=thr*flx_r+thl*flx_l+(1.0d0-thl-thr)* &
         (thm*(flx_l+sl*(st_l-u_l)) + (1.0d0-thm)*(flx_r+sr*(st_r-u_r)))

  end subroutine hllc_flux


  !! Ref: OpenMHD code developped by S. Zenitani (ver 20150325)
  subroutine hlld_flux(v_l,v_r,flx_tmp,nvar,gm)  
    implicit none
    integer,intent(in)::nvar
    double precision,intent(in)::v_l(nvar),v_r(nvar),gm      
    double precision,intent(out):: flx_tmp(nvar)
    double precision :: u_l(nvar),b2_l,pt_l,vbL,flx_l(nvar),vf_l,sl,ro_l1,sal &
         ,u_l1(nvar),vy_l1,vz_l1,sl1,ro_ls
    double precision :: u_r(nvar),b2_r,pt_r,vbR,flx_r(nvar),vf_r,sr,ro_r1,sar &
         ,u_r1(nvar),vy_r1,vz_r1,sr1,ro_rs
    double precision :: f1,f2,u_hll(nvar),sm,pt,vy_2,vz_2,u2(nvar)
    double precision, parameter :: hllg_factor = 1.001d0
    double precision :: ch
    
   flx_tmp = 0.d0 ; flx_l = 0.d0 ; flx_r = 0.d0

    call v2u_1(v_l,u_l,nvar_m,1)
    call v2u_1(v_r,u_r,nvar_m,1)    
    
    b2_l=(v_l(6)**2+v_l(7)**2+v_l(8)**2)
    pt_l=v_l(5)+0.5d0*b2_l
    vbL = v_l(2)*v_l(6) + v_l(3)*v_l(7) + v_l(4)*v_l(8)

    flx_l(1) = v_l(1)*v_l(2)
    flx_l(2) = v_l(1)*v_l(2)**2 + v_l(5) +0.5d0*b2_l - v_l(6)**2
    flx_l(3) = v_l(1)*v_l(2)*v_l(3) - v_l(6)*v_l(7)
    flx_l(4) = v_l(1)*v_l(2)*v_l(4) - v_l(6)*v_l(8)
    flx_l(5) = (u_l(5) + pt_l)*v_l(2) - vbL*v_l(6)
    flx_l(6) = 0.d0
    flx_l(7) = v_l(2)*v_l(7) - v_l(6)*v_l(3)
    flx_l(8) = v_l(2)*v_l(8) - v_l(6)*v_l(4)    

    b2_r=(v_r(6)**2+v_r(7)**2+v_r(8)**2)
    pt_r=v_r(5)+0.5d0*b2_r
    vbR = v_r(2)*v_r(6) + v_r(3)*v_r(7) + v_r(4)*v_r(8)

    flx_r(1) = v_r(1)*v_r(2)
    flx_r(2) = v_r(1)*v_r(2)**2 + v_r(5) +0.5d0*b2_r - v_r(6)**2
    flx_r(3) = v_r(1)*v_r(2)*v_r(3) - v_r(6)*v_r(7)
    flx_r(4) = v_r(1)*v_r(2)*v_r(4) - v_r(6)*v_r(8)
    flx_r(5) = (u_r(5) + pt_r)*v_r(2) - vbR*v_r(6)
    flx_r(6) = 0.d0
    flx_r(7) = v_r(2)*v_r(7) - v_r(6)*v_r(3)
    flx_r(8) = v_r(2)*v_r(8) - v_r(6)*v_r(4)    

    !! fast mode wave speed
    f1 = gm * v_l(5)
    f2 = 4.d0*f1*v_l(6)**2
    vf_l = sqrt( ( (f1+b2_l) &
         + sqrt( max((f1+b2_l)**2-f2,0.d0))) / (2.d0*v_l(1)) )

    f1 = gm * v_r(5)
    f2 = 4.d0*f1*v_r(6)**2
    vf_r = sqrt( ( (f1+b2_r) &
         + sqrt( max((f1+b2_r)**2-f2,0.d0))) / (2.d0*v_r(1)) )
    

    !! Wave speed estimation
    sl = min(v_l(2),v_r(2)) - max(vf_l,vf_r)
    sr = max(v_l(2),v_r(2)) + max(vf_l,vf_r)
    
    
    !If super-sonic flow ----------------------------------------------------
    if (sl.ge.0.0d0) then
       flx_tmp = flx_l
    else if (sr .le. 0.0d0) then
       flx_tmp = flx_r
       !-----------------------------------------------------------------------
       !If sub-sonic flow ---------------------------------------------------- 
    else
       f1 = 1.d0/(sr-sl)
       u_hll(1:8) = f1*( sr*u_r(1:8)-sl*u_l(1:8)-flx_r(1:8)+flx_l(1:8) ) 

       !! entropy wave speed
       sm = u_hll(2)/u_hll(1)

       pt = pt_l + v_l(1)*(sl-v_l(2))*(sm-v_l(2))
       ro_l1 = v_l(1) * (sl-v_l(2)) / (sl-sm) ! ro(L*)
       ro_r1 = v_r(1) * (sr-v_r(2)) / (sr-sm) ! ro(R*)       

       ! For logical consistency, 
       ! we employ bx_hll as b_x in the intermediate states
       sal = abs(u_hll(6))/sqrt(ro_l1) !! Alfven wave (L)
       sar = abs(u_hll(6))/sqrt(ro_r1) !! Alfven wave (R)

       !=== revert to HLLC-G ===
       if( (sl.ge.(sm-hllg_factor*sal)) .or. &
            (( sm+hllg_factor*sar) .ge. sr )) then
          vy_2 = u_hll(3) / u_hll(1) ! vy_HLL
          vz_2 = u_hll(4) / u_hll(1) ! vz_HLL          

          ! F = F(L*)
          if (sm .ge. 0.d0) then
             u2(5) = ( (sl-v_l(2))*u_l(5) - pt_l*v_l(2) + pt*sm &
                  + u_hll(6)*( vbL &
                  - sm*u_hll(6) - vy_2*u_hll(7) - vz_2*u_hll(8) )) &
                  / (sl-sm)
             u2(2) = ro_l1 * sm
             u2(3) = ro_l1 * vy_2
             u2(4) = ro_l1 * vz_2
             
             flx_tmp(1)    = sl * ( ro_l1      - v_l(1)   ) + flx_l(1)
             flx_tmp(2:5)  = sl * ( u2(2:5)    - u_l(2:5) ) + flx_l(2:5)
             flx_tmp(6:8)  = sl * ( u_hll(6:8) - v_l(6:8) ) + flx_l(6:8)
          ! F = F(R*)
          else
             u2(5) = ( (sr-v_r(2))*u_r(5) - pt_r*v_r(2) + pt*sm &
                  + u_hll(6)*( vbR &
                  - sm*u_hll(6) - vy_2*u_hll(7) - vz_2*u_hll(8) )) &
                  / (sr-sm)
             u2(2) = ro_r1 * sm
             u2(3) = ro_r1 * vy_2
             u2(4) = ro_r1 * vz_2
             
             flx_tmp(1)    = sr * ( ro_r1      - v_r(1) )   + flx_r(1)
             flx_tmp(2:5)  = sr * ( u2(2:5)    - u_r(2:5) ) + flx_r(2:5)
             flx_tmp(6:8)  = sr * ( u_hll(6:8) - v_r(6:8) ) + flx_r(6:8)
          endif

          !! HLLD flux
       else
          !! intermediate state (L*)
          ! HLLD denominator
          f1 = 1.d0 / ( v_l(1)*(sl-v_l(2))*(sl-sm) - u_hll(6)**2 )
          u_l1(6) = u_hll(6)
          u_l1(7) = v_l(7) * f1 * ( v_l(1)*(sl-v_l(2))**2 - u_hll(6)**2 )
          u_l1(8) = v_l(8) * f1 * ( v_l(1)*(sl-v_l(2))**2 - u_hll(6)**2 )
          vy_l1   = v_l(3) - f1 * u_hll(6)*v_l(7)*(sm-v_l(2))
          vz_l1   = v_l(4) - f1 * u_hll(6)*v_l(8)*(sm-v_l(2))          

          u_l1(1) = ro_l1
          u_l1(5) = ( ( sl-v_l(2) )*u_l(5) - pt_l*v_l(2) + pt*sm &
               + u_hll(6)*( vbL - sm*u_hll(6)-vy_l1*u_l1(7)-vz_l1*u_l1(8) ) ) &
               / (sl-sm)
          u_l1(2) = ro_l1 * sm
          u_l1(3) = ro_l1 * vy_l1
          u_l1(4) = ro_l1 * vz_l1

          !! intermediate state (R*)
          !! HLLD denominator
          f1 = 1.d0 / ( v_r(1)*(sr-v_r(2))*(sr-sm) - u_hll(6)**2 )
          u_r1(6) = u_hll(6)
          u_r1(7) = v_r(7) * f1 * ( v_r(1)*(sr-v_r(2))**2 - u_hll(6)**2 )
          u_r1(8) = v_r(8) * f1 * ( v_r(1)*(sr-v_r(2))**2 - u_hll(6)**2 )
          vy_r1   = v_r(3) - f1 * u_hll(6)*v_r(7)*(sm-v_r(2))
          vz_r1   = v_r(4) - f1 * u_hll(6)*v_r(8)*(sm-v_r(2))          

          u_r1(1) = ro_r1
          u_r1(5) = ( ( sr-v_r(2) )*u_r(5) - pt_r*v_r(2) + pt*sm &
               + u_hll(6)*( vbR - sm*u_hll(6)-vy_r1*u_r1(7)-vz_r1*u_r1(8) ) ) &
               / (sr-sm)
          u_r1(2) = ro_r1 * sm
          u_r1(3) = ro_r1 * vy_r1
          u_r1(4) = ro_r1 * vz_r1

          !! rotational discontinuity
          sl1 = sm - sal
          sr1 = sm + sar
          !! F = F(L*)
          if(sl1 .ge. 0.d0) then
             flx_tmp(1)   = sl * (u_l1(1)   - v_l(1)  ) + flx_l(1)
             flx_tmp(2:5) = sl * (u_l1(2:5) - u_l(2:5)) + flx_l(2:5)
             flx_tmp(6:8) = sl * (u_l1(6:8) - v_l(6:8)) + flx_l(6:8) 
          !! F = F(R*)
          elseif (sr1 .le. 0.d0) then
             flx_tmp(1)   = sr * (u_r1(1)   - v_r(1)  ) + flx_r(1)
             flx_tmp(2:5) = sr * (u_r1(2:5) - u_r(2:5)) + flx_r(2:5)
             flx_tmp(6:8) = sr * (u_r1(6:8) - v_r(6:8)) + flx_r(6:8) 

          !! Central states
          !! Question: can we really rule out U_hll(bx) == 0 ?
          else
             ro_ls = sqrt(ro_l1)
             ro_rs = sqrt(ro_r1)
             f1 = 1.d0 / (ro_ls + ro_rs)
             f2 = dsign(1.d0, u_hll(6))
             vy_2  = f1 * ( ro_ls*vy_l1 + ro_rs*vy_r1 + (u_r1(7)-u_l1(7))*f2 )
             vz_2  = f1 * ( ro_ls*vz_l1 + ro_rs*vz_r1 + (u_r1(8)-u_l1(8))*f2 )
             u2(6) = u_hll(6)
             u2(7) = f1 * ( ro_ls*u_r1(7) + ro_rs*u_l1(7) &
                  + ro_ls*ro_rs*( vy_r1-vy_l1 )*f2 )
             u2(8) = f1 * ( ro_ls*u_r1(8) + ro_rs*u_l1(8) &
                  + ro_ls*ro_rs*( vz_r1-vz_l1 )*f2 )

             !! F = F(L**)
             if(sm .ge. 0.d0) then
                u2(1) = ro_l1
                u2(5) = u_l1(5) - ro_ls * (vy_l1*u_l1(7) + vz_l1*u_l1(8) &
                     - vy_2*u2(7) - vz_2*u2(8) ) * f2
                u2(2) = ro_l1 * sm
                u2(3) = ro_l1 * vy_2
                u2(4) = ro_l1 * vz_2 
                
                flx_tmp(1)   = sl1*u2(1)   - (sl1-sl)*u_l1(1)   - sl*v_l(1)   + flx_l(1)
                flx_tmp(2:5) = sl1*u2(2:5) - (sl1-sl)*u_l1(2:5) - sl*u_l(2:5) + flx_l(2:5)
                flx_tmp(6:8) = sl1*u2(6:8) - (sl1-sl)*u_l1(6:8) - sl*v_l(6:8) + flx_l(6:8)
             !! F = F(R**)
             else
                u2(1) = ro_r1
                u2(5) = u_r1(5) + ro_rs * (vy_r1*u_r1(7) + vz_r1*u_r1(8) &
                     - vy_2*u2(7) - vz_2*u2(8) ) * f2
                u2(2) = ro_r1 * sm
                u2(3) = ro_r1 * vy_2
                u2(4) = ro_r1 * vz_2 
                
                flx_tmp(1)   = sr1*u2(1)   - (sr1-sr)*u_r1(1)   - sr*v_r(1)   + flx_r(1)
                flx_tmp(2:5) = sr1*u2(2:5) - (sr1-sr)*u_r1(2:5) - sr*u_r(2:5) + flx_r(2:5)
                flx_tmp(6:8) = sr1*u2(6:8) - (sr1-sr)*u_r1(6:8) - sr*v_r(6:8) + flx_r(6:8)
             endif
          endif
       endif
    endif
    
    !! divB cleaning
    if(flag_divb.eq.1) then
       ch = cmax
       flx_tmp(6) = 0.5d0 * ( u_l(9)+u_r(9) - ch * (v_r(6)-v_l(6)) )
       flx_tmp(9) = 0.5d0 * ( ch*ch*(v_l(6)+v_r(6)) - ch * (v_r(9)-v_l(9)) )       
    endif
                
  end subroutine hlld_flux




  ! subroutine hlld_flux(v_l,v_r,flx_tmp,nvar,gm)  
  !   implicit none
  !   integer,intent(in)::nvar
  !   double precision,intent(in)::v_l(nvar),v_r(nvar),gm      
  !   double precision,intent(out):: flx_tmp(nvar)
  !   double precision pr_l,pr_r!,flx_l(nvar),flx_r(nvar)
  !   double precision vn_l,vn_r,sl,sls,sm,srs,sr,pt_m
  !   double precision vt1_l,vt2_l,bt1_l,bt2_l,vt1_s_l,vt2_s_l,bt1_s_l,bt2_s_l
  !   double precision vt1_r,vt2_r,bt1_r,bt2_r,vt1_s_r,vt2_s_r,bt1_s_r,bt2_s_r
  !   double precision de_l,de_s_l,pt_l,cs_l,va_l,vf_l,e_l,e_s_l,e_ss_l
  !   double precision de_r,de_s_r,pt_r,cs_r,va_r,vf_r,e_r,e_s_r,e_ss_r
  !   double precision tmp,tmp_l,tmp_r,fact,al1,ar1,al2,ar2,al3,ar3,bn_pml,num3
  !   double precision bn_l,vt1_ss_l,bt1_ss_l,vt2_ss_l,bt2_ss_l,b2_l
  !   double precision bn_r,vt1_ss_r,bt1_ss_r,vt2_ss_r,bt2_ss_r,b2_r,bn
    

  !   ! pr_l=(gm-1.0d0)*(u_l(5) &
  !   !      -0.5d0*(u_l(2)**2+u_l(3)**2+u_l(4)**2)/u_l(1)  &
  !   !      -0.5d0*(u_l(6)**2+u_l(7)**2+u_l(8)**2))
  !   ! pr_r=(gm-1.0d0)*(u_r(5) &
  !   !      -0.5d0*(u_r(2)**2+u_r(3)**2+u_r(4)**2)/u_r(1)  &
  !   !      -0.5d0*(u_r(6)**2+u_r(7)**2+u_r(8)**2))    

  !   ! de_l=u_l(1)
  !   ! de_r=u_r(1)

  !   ! vn_l=u_l(2)/de_l
  !   ! vn_r=u_r(2)/de_r

  !   ! vt1_l=u_l(3)/de_l ; vt1_r=u_r(3)/de_r
  !   ! vt2_l=u_l(4)/de_l ; vt2_r=u_r(4)/de_r
  !   ! !    bn=0.5d0*(u_l(6)+u_r(6))
  !   ! bn_l=u_l(6) ;bn_r=u_r(6) ; bn=0.5d0*(bn_l+bn_r)
  !   ! bt1_l=u_l(7) ;bt1_r=u_r(7)
  !   ! bt2_l=u_l(8) ;bt2_r=u_r(8)
  !   ! e_l=u_l(5)
  !   ! e_r=u_r(5)

  !   de_l = v_l(1)
  !   vn_l = v_l(2)
  !   vt1_l= v_l(3)
  !   vt2_l= v_l(4)
  !   pr_l = v_l(5)
  !   bn_l = v_l(6)
  !   bt1_l= v_l(7)
  !   bt2_l= v_l(8)
  !   e_l  = pr_l/(gm-1.d0) &
  !        + 0.5d0*de_l*(vn_l**2+vt1_l**2+vt2_l**2) &
  !        + 0.5d0*(bn_l**2+bt1_l**2+bt2_l**2)

  !   de_r = v_r(1)
  !   vn_r = v_r(2)
  !   vt1_r= v_r(3)
  !   vt2_r= v_r(4)
  !   pr_r = v_r(5)
  !   bn_r = v_r(6)
  !   bt1_r= v_r(7)
  !   bt2_r= v_r(8)
  !   e_r  = pr_r/(gm-1.d0) &
  !        + 0.5d0*de_r*(vn_r**2+vt1_r**2+vt2_r**2) &
  !        + 0.5d0*(bn_r**2+bt1_r**2+bt2_r**2)
    
  !   bn = 0.5d0*(bn_l+bn_r)
    
  !   cs_l=gm*pr_l/de_l
  !   cs_r=gm*pr_r/de_r

  !   b2_l=(bn_l**2+bt1_l**2+bt2_l**2)
  !   b2_r=(bn_r**2+bt1_r**2+bt2_r**2)

  !   va_l=b2_l/de_l
  !   va_r=b2_r/de_r

  !   pt_l=pr_l+0.5d0*b2_l
  !   pt_r=pr_r+0.5d0*b2_r

  !   vf_l=sqrt((cs_l+va_l + &
  !        sqrt(abs((cs_l+va_l)**2-4.0d0*cs_l*bn_l**2/de_l)))/2.0d0)
  !   vf_r=sqrt((cs_r+va_r + &
  !        sqrt(abs((cs_r+va_r)**2-4.0d0*cs_r*bn_r**2/de_r)))/2.0d0)


  !   sl=acc*min(vn_l-vf_l,vn_r-vf_r)
  !   sr=acc*max(vn_l+vf_l,vn_r+vf_r)

  !   !If super-sonic flow ----------------------------------------------------
  !   if (sl.gt.0.0d0) then
  !      flx_tmp(1)=de_l*vn_l
  !      flx_tmp(2)=de_l*vn_l**2+pt_l-bn**2
  !      flx_tmp(3)=de_l*vn_l*vt1_l-bn*bt1_l
  !      flx_tmp(4)=de_l*vn_l*vt2_l-bn*bt2_l
  !      flx_tmp(5)=vn_l*(e_l+pt_l)-bn*(vn_l*bn_l+vt1_l*bt1_l+vt2_l*bt2_l)
  !      flx_tmp(6)=0.0d0
  !      flx_tmp(7)=(vn_l*bt1_l-vt1_l*bn)
  !      flx_tmp(8)=(vn_l*bt2_l-vt2_l*bn)
  !   else if (sr .lt. 0.0d0) then
  !      flx_tmp(1)=de_r*vn_r
  !      flx_tmp(2)=de_r*vn_r**2+pt_r-bn**2
  !      flx_tmp(3)=de_r*vn_r*vt1_r-bn*bt1_r
  !      flx_tmp(4)=de_r*vn_r*vt2_r-bn*bt2_r
  !      flx_tmp(5)=vn_r*(e_r+pt_r)-bn*(vn_r*bn_r+vt1_r*bt1_r+vt2_r*bt2_r)
  !      flx_tmp(6)=0.0d0
  !      flx_tmp(7)=(vn_r*bt1_r-vt1_r*bn)
  !      flx_tmp(8)=(vn_r*bt2_r-vt2_r*bn)       
  !      !-----------------------------------------------------------------------
  !      !If sub-sonic flow ---------------------------------------------------- 
  !   else
  !      tmp=((sr-vn_r)*de_r-(sl-vn_l)*de_l)
  !      ar1=de_r*(sr-vn_r)
  !      al1=de_l*(sl-vn_l)       
  !      sm=(ar1*vn_r-al1*vn_l+(pt_l-pt_r))/tmp       
  !      pt_m=((ar1*pt_l-al1*pt_r)+de_l*de_r*(sr-vn_r)*(sl-vn_l)*(vn_r-vn_l))/tmp

  !      de_s_l=de_l*(sl-vn_l)/(sl-sm)
  !      sls=sm-abs(bn_l)/de_s_l

  !      de_s_r=de_r*(sr-vn_r)/(sr-sm)
  !      srs=sm+abs(bn_r)/de_s_r

  !      tmp_l=(de_l*(sl-vn_l)*(sl-sm)-bn_l**2)       
  !      if(tmp_l .ne. 0.0d0) then
  !         fact=abs(tanh(tmp_l/epsi))**2
  !         !       fact=1.0d0
  !         al2=(sm-vn_l)*bn_l/tmp_l
  !         vt1_s_l=vt1_l-al2*bt1_l*fact
  !         bt1_s_l=bt1_l*(de_l*(sl-vn_l)**2-bn_l**2)/tmp_l*fact
  !         vt2_s_l=vt2_l-al2*bt2_l*fact
  !         bt2_s_l=bt2_l*(de_l*(sl-vn_l)**2-bn_l**2)/tmp_l*fact
  !      else
  !         vt1_s_l=vt1_l;bt1_s_l=0.0d0
  !         vt2_s_l=vt2_l;bt2_s_l=0.0d0
  !      endif
  !      tmp_r=(de_r*(sr-vn_r)*(sr-sm)-bn_r**2)
  !      if(tmp_r .ne. 0.0d0) then
  !         fact=abs(tanh(tmp_r/epsi))**2
  !         ar2=(sm-vn_r)*bn_r/tmp_r
  !         vt1_s_r=vt1_r-ar2*bt1_r*fact
  !         bt1_s_r=bt1_r*(de_r*(sr-vn_r)**2-bn_r**2)/tmp_r*fact
  !         vt2_s_r=vt2_r-ar2*bt2_r*fact
  !         bt2_s_r=bt2_r*(de_r*(sr-vn_r)**2-bn_r**2)/tmp_r*fact
  !      else
  !         vt1_s_r=vt1_l;bt1_s_r=0.0d0
  !         vt2_s_r=vt2_l;bt2_s_l=0.0d0
  !      endif
  !      e_s_l=((sl-vn_l)*e_l-pt_l*vn_l+pt_m*sm)/(sl-sm) &
  !           +bn_l*(vn_l*bn_l+vt1_l*bt1_l+vt2_l*bt2_l &
  !           -(sm*bn_l+vt1_s_l*bt1_s_l+vt2_s_l*bt2_s_l))/(sl-sm)
  !      e_s_r=((sr-vn_r)*e_r-pt_r*vn_r+pt_m*sm)/(sr-sm) &
  !           +bn_r*(vn_r*bn_r+vt1_r*bt1_r+vt2_r*bt2_r &
  !           -(sm*bn_r+vt1_s_r*bt1_s_r+vt2_s_r*bt2_s_r))/(sr-sm)


  !      !If rotational discon ----------------------------------------------
  !      if (sls .gt. 0.0d0) then
  !         flx_tmp(1)=de_s_l*sm
  !         flx_tmp(2)=de_s_l*sm**2+pt_m-bn**2
  !         flx_tmp(3)=de_s_l*sm*vt1_s_l-bn*bt1_s_l
  !         flx_tmp(4)=de_s_l*sm*vt2_s_l-bn*bt2_s_l
  !         flx_tmp(5)=sm*(e_s_l+pt_m)-bn_l*(sm*bn+vt1_s_l*bt1_s_l+vt2_s_l*bt2_s_l)
  !         flx_tmp(6)=0.0d0
  !         flx_tmp(7)=(sm*bt1_s_l-vt1_s_l*bn)
  !         flx_tmp(8)=(sm*bt2_s_l-vt2_s_l*bn)
  !      else if (srs .lt. 0.0d0) then
  !         flx_tmp(1)=de_s_r*sm
  !         flx_tmp(2)=de_s_r*sm**2+pt_m-bn**2
  !         flx_tmp(3)=de_s_r*sm*vt1_s_r-bn*bt1_s_r
  !         flx_tmp(4)=de_s_r*sm*vt2_s_r-bn*bt2_s_r
  !         flx_tmp(5)=sm*(e_s_r+pt_m)-bn_r*(sm*bn+vt1_s_r*bt1_s_r+vt2_s_r*bt2_s_r)
  !         flx_tmp(6)=0.0d0
  !         flx_tmp(7)=(sm*bt1_s_r-vt1_s_r*bn)
  !         flx_tmp(8)=(sm*bt2_s_r-vt2_s_r*bn)
  !         !----------------------------------------------------------
  !      else          
  !         num3=1.0d0/(sqrt(de_s_l)+sqrt(de_s_r))
  !         al3=sqrt(de_s_l)
  !         ar3=sqrt(de_s_r)

  !         bn_pml=sign(1.0d0,bn)

  !         vt1_ss_l=(al3*vt1_s_l+ar3*vt1_s_r+(bt1_s_r-bt1_s_l)*bn_pml)*num3
  !         bt1_ss_l=(al3*bt1_s_l+ar3*bt1_s_r+ &
  !              sqrt(de_s_l*de_s_r)*(vt1_s_r-vt1_s_l)*bn_pml)*num3

  !         vt2_ss_l=(al3*vt2_s_l+ar3*vt2_s_r+(bt2_s_r-bt2_s_l)*bn_pml)*num3
  !         bt2_ss_l=(al3*bt2_s_l+ar3*bt2_s_r+ &
  !              sqrt(de_s_l*de_s_r)*(vt2_s_r-vt2_s_l)*bn_pml)*num3

  !         vt1_ss_r=vt1_ss_l
  !         vt2_ss_r=vt2_ss_l
  !         bt1_ss_r=bt1_ss_l
  !         bt2_ss_r=bt2_ss_l

  !         e_ss_l=e_s_l-sqrt(de_s_l)*(sm*bn_l+vt1_s_l*bt1_s_l+vt2_s_l*bt2_s_l &
  !              -(sm*bn_l+vt1_ss_l*bt1_ss_l+vt2_ss_l*bt2_ss_l))*bn_pml
  !         e_ss_r=e_s_r+sqrt(de_s_r)*(sm*bn_r+vt1_s_r*bt1_s_r+vt2_s_r*bt2_s_r &
  !              -(sm*bn_r+vt1_ss_r*bt1_ss_r+vt2_ss_r*bt2_ss_r))*bn_pml

  !         if (sm .gt. 0.0d0) then
  !            flx_tmp(1)=de_s_l*sm
  !            flx_tmp(2)=de_s_l*sm**2+pt_m-bn**2
  !            flx_tmp(3)=de_s_l*sm*vt1_ss_l-bn*bt1_ss_l
  !            flx_tmp(4)=de_s_l*sm*vt2_ss_l-bn*bt2_ss_l
  !            flx_tmp(5)=sm*(e_ss_l+pt_m)-bn*(sm*bn+vt1_ss_l*bt1_ss_l+vt2_ss_l*bt2_ss_l)
  !            flx_tmp(6)=0.0d0
  !            flx_tmp(7)=(sm*bt1_ss_l-vt1_ss_l*bn)
  !            flx_tmp(8)=(sm*bt2_ss_l-vt2_ss_l*bn)
  !         else
  !            flx_tmp(1)=de_s_r*sm
  !            flx_tmp(2)=de_s_r*sm**2+pt_m-bn**2
  !            flx_tmp(3)=de_s_r*sm*vt1_ss_r-bn*bt1_ss_r
  !            flx_tmp(4)=de_s_r*sm*vt2_ss_r-bn*bt2_ss_r
  !            flx_tmp(5)=sm*(e_ss_r+pt_m)-bn*(sm*bn+vt1_ss_r*bt1_ss_r+vt2_ss_r*bt2_ss_r)
  !            flx_tmp(6)=0.0d0
  !            flx_tmp(7)=(sm*bt1_ss_r-vt1_ss_r*bn)
  !            flx_tmp(8)=(sm*bt2_ss_r-vt2_ss_r*bn)
  !         endif
  !      endif
  !   endif
  ! end subroutine hlld_flux
  
  
  subroutine u2v_hd(V_h,U_h)
    double precision, intent(in) :: U_h(ix,jx,kx,nvar_h)
    double precision, intent(out) :: V_h(ix,jx,kx,nvar_h)
    double precision :: roi(ix,jx,kx)
    
    !! ro
    V_h(:,:,:,1) = max(U_h(:,:,:,1),ro_lim)
    roi = 1.d0 / V_h(:,:,:,1)
    !! vx
    V_h(:,:,:,2) = U_h(:,:,:,2)*roi
    !! vy
    V_h(:,:,:,3) = U_h(:,:,:,3)*roi
    !! vz
    V_h(:,:,:,4) = U_h(:,:,:,4)*roi
    !! pr
    V_h(:,:,:,5) = (gm-1.d0)*( U_h(:,:,:,5) &
         -0.5d0*( U_h(:,:,:,2)**2+U_h(:,:,:,3)**2+U_h(:,:,:,4)**2 )*roi )
    
    where(V_h(:,:,:,5).lt.pr_lim) V_h(:,:,:,5) = pr_lim

  end subroutine u2v_hd


  subroutine u2v_mhd(V_m,U_m)
    double precision, intent(in) :: U_m(ix,jx,kx,nvar_m)
    double precision, intent(out) :: V_m(ix,jx,kx,nvar_m)
    double precision :: roi(ix,jx,kx)
    
    !! ro
    V_m(:,:,:,1) = max(U_m(:,:,:,1),ro_lim)
    roi = 1.d0 / V_m(:,:,:,1)
    !! vx
    V_m(:,:,:,2) = U_m(:,:,:,2)*roi
    !! vy
    V_m(:,:,:,3) = U_m(:,:,:,3)*roi
    !! vz
    V_m(:,:,:,4) = U_m(:,:,:,4)*roi
    !! pr
    V_m(:,:,:,5) = (gm-1.d0)*( U_m(:,:,:,5) &
         -0.5d0*( U_m(:,:,:,2)**2+U_m(:,:,:,3)**2+U_m(:,:,:,4)**2 )*roi &
         -0.5d0*( U_m(:,:,:,6)**2+U_m(:,:,:,7)**2+U_m(:,:,:,8)**2 ) )
    
    where(V_m(:,:,:,5).lt.pr_lim) V_m(:,:,:,5) = pr_lim

    !! bx,by,bz
    V_m(:,:,:,6:8) = U_m(:,:,:,6:8)

    !! psi
    if(flag_divb.eq.1) V_m(:,:,:,9) = U_m(:,:,:,9)
    
  end subroutine u2v_mhd


  subroutine v2u_1(v1,u1,nvar,flag0)
    integer, intent(in) :: nvar,flag0
    double precision, intent(in) :: v1(nvar)
    double precision, intent(out) :: u1(nvar)
    
    select case(flag0)
    case(0) !! HD
       u1(1) = v1(1)
       u1(2) = v1(1)*v1(2)
       u1(3) = v1(1)*v1(3)
       u1(4) = v1(1)*v1(4)
       u1(5) = v1(5)/(gm-1.d0) &
            + 0.5d0*v1(1)*(v1(2)**2 + v1(3)**2 + v1(4)**2)
    case(1) !! MHD
       u1(1) = v1(1)
       u1(2) = v1(1)*v1(2)
       u1(3) = v1(1)*v1(3)
       u1(4) = v1(1)*v1(4)
       u1(6:8) = v1(6:8)
       u1(5) = v1(5)/(gm-1.d0) &
            + 0.5d0*v1(1)*(v1(2)**2 + v1(3)**2 + v1(4)**2) &
            + 0.5d0*(v1(6)**2 + v1(7)**2 + v1(8)**2)
       !! psi
       if(flag_divb.eq.1) u1(9) = v1(9)

    end select
  end subroutine v2u_1

  subroutine v2f_1(v1,f1,nvar,flag0)
    integer, intent(in) :: nvar,flag0
    double precision, intent(in) :: v1(nvar)
    double precision, intent(out) :: f1(nvar)
    double precision :: en1,bb1

    select case(flag0)
    case(0) !! HD
       en1 = v1(5)/(gm-1.d0) + 0.5d0*v1(1)*( v1(2)**2 + v1(3)**2 + v1(4)**2 )
       f1(1) = v1(1)*v1(2)
       f1(2) = v1(1)*v1(2)**2 
       f1(3) = v1(1)*v1(3)*v1(2)
       f1(4) = v1(1)*v1(4)*v1(2)
       f1(5) = ( en1 + v1(5) ) * v1(2) 
    case(1) !! MHD
       bb1 = v1(6)**2 + v1(7)**2 + v1(8)**2
       en1 = v1(5)/(gm-1.d0) + 0.5d0*v1(1)*( v1(2)**2 + v1(3)**2 + v1(4)**2 ) &
            + 0.5d0*bb1
       f1(1) = v1(1)*v1(2)
       f1(2) = v1(1)*v1(2)**2 + v1(5)+0.5d0*bb1 - v1(6)**2
       f1(3) = v1(1)*v1(3)*v1(2) - v1(7)*v1(6)
       f1(4) = v1(1)*v1(4)*v1(2) - v1(8)*v1(6)
       f1(6) = 0.d0
       f1(7) = v1(2)*v1(7) - v1(3)*v1(6)
       f1(8) = v1(2)*v1(8) - v1(4)*v1(6)
       !       f1(5) = ( en1 + v1(5) + 0.5d0*bb1 ) * v1(2) &
       f1(5) = ( en1 + v1(5) + bb1 ) * v1(2) &       
            - v1(6) * ( v1(6)*v1(2)+v1(7)*v1(3)+v1(8)*v1(4) )
    end select
  end subroutine v2f_1

  
end module HLL_rot

