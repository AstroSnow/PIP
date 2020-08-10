module Visc_rot
use globalvar,only:ix,jx,kx,nvar_h,nvar_m,x,y,z,debug_parameter,&
ndim,s_order,dxc,dyc,dzc,my_rank,j_cri,vd_cri,nu_0,flag_visc,flag_mhd,&
flag_pip,flag_bnd,neighbor,visc
use Util_rot,only:get_laplace,get_grad,get_divergence
implicit none  
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize_visc(flag_visc)
	integer,intent(inout)::flag_visc
	if (flag_visc.eq.0) return
	       allocate(visc(ix,jx,kx,3))
!	if(my_rank.eq.0) print*,'visc test '
end subroutine initialize_visc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_visc(U,nu)
	double precision,intent(in)::U(ix,jx,kx,nvar_m) 
	double precision,intent(out)::nu(ix,jx,kx,3)
	double precision L(ix,jx,kx,3),L2(ix,jx,kx,3)
	double precision V(ix,jx,kx,3),L3(ix,jx,kx)
	double precision tmp
	integer dir,upper_lower
	integer i,j,k    

	V(:,:,:,1)=U(:,:,:,2)/U(:,:,:,1)
	V(:,:,:,2)=U(:,:,:,3)/U(:,:,:,1)
	V(:,:,:,3)=U(:,:,:,4)/U(:,:,:,1)

!	print*,V(54,1:5,1,1)

	!get laplacian of velocity
	call get_laplace(ix,jx,kx,ndim,s_order,dxc,dyc,dzc,V,L)
!if (my_rank.eq.0) print*, dyc(1:5)
!if (my_rank.eq.0) print*, L(55,1:5,:,2)

	call get_divergence(ix,jx,kx,ndim,s_order,dxc,dyc,dzc,V,L3)
	call get_grad(ix,jx,kx,ndim,s_order,dxc,dyc,dzc,L3,L2)

	L=L+1.0d0/3.0d0*L2
!if (my_rank.eq.0) print*, L(55,1:5,:,2)
!stop

	do i=1,2*ndim
		dir=(i+1)/2
		upper_lower=mod(i+1,2) 
		if(neighbor(i).eq.-1) then
			call visc_bc(L,V,dir,upper_lower)
		endif
		if(neighbor(i).ne.-1) then
			call visc_ibc(L,V,dir,upper_lower)
		endif
	enddo

!!!Set the third viscosity to zero
if (ndim .eq. 1) then
L(:,:,:,2)=0.0d0
L(:,:,:,3)=0.0d0
endif
if (ndim .eq. 2) then
L(:,:,:,3)=0.0d0
endif


	select case(flag_visc)
	case(1)!uniform viscosity 
		nu(:,:,:,:)=nu_0*L
	case(2)!artificially localized viscosity
		print*,'NOT PROGRAMMED YET'    
		stop                            
	end select
end subroutine set_visc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine source_visc(S_h,S_m,U_h,U_m)
	double precision,intent(inout)::S_h(ix,jx,kx,nvar_h),S_m(ix,jx,kx,nvar_m)
	double precision,intent(in)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)
!	double precision visc(ix,jx,kx,3)
	integer n
!NEED FRACTIONAL VISCOSITY!!!
	double precision::xif(ix,jx,kx)
	xif=U_m(:,:,:,1)/(U_m(:,:,:,1)+U_h(:,:,:,1))

	if(flag_mhd.eq.0.or.flag_pip.eq.1) then
		call set_visc(U_h,visc)

!if (my_rank.eq.0) print*, visc(55,1:5,:,2)
!stop
		S_h(:,:,:,2) = S_h(:,:,:,2) &
		    + (1.0d0-xif)*visc(:,:,:,1)
		S_h(:,:,:,3) = S_h(:,:,:,3) &
		    + (1.0d0-xif)*visc(:,:,:,2)
		S_h(:,:,:,4) = S_h(:,:,:,4) &
		    + (1.0d0-xif)*visc(:,:,:,3)
	endif
	if(flag_mhd.eq.1) then
		call set_visc(U_m,visc)

!if (my_rank.eq.0) print*, visc(55,1:5,:,2)
!stop
		S_m(:,:,:,2) = S_m(:,:,:,2) &
		    + xif*visc(:,:,:,1)
		S_m(:,:,:,3) = S_m(:,:,:,3) &
		    + xif*visc(:,:,:,2)
		S_m(:,:,:,4) = S_m(:,:,:,4) &
		    + xif*visc(:,:,:,3)
	endif
end subroutine source_visc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine visc_bc(L,V,dir,upper_lower)
	double precision,intent(inout)::L(ix,jx,kx,3)
	double precision,intent(in)::V(ix,jx,kx,3)
	integer, intent(in)::dir,upper_lower
	if(s_order.eq.4) then
		print*,'NOT PROGRAMMED'
		stop
	else
		if (dir .eq.1) then
			if (upper_lower .eq. 0) then
			L(1,:,:,1)=L(3,:,:,1);L(2,:,:,1)=L(3,:,:,1) 
			L(1,:,:,2)=L(3,:,:,2);L(2,:,:,2)=L(3,:,:,2)
			L(1,:,:,3)=L(3,:,:,3);L(2,:,:,3)=L(3,:,:,3)
			elseif (upper_lower .eq. 1) then
			L(ix,:,:,1)=L(ix-2,:,:,1);L(ix-1,:,:,1)=L(ix-2,:,:,1) 
			L(ix,:,:,2)=L(ix-2,:,:,2);L(ix-1,:,:,2)=L(ix-2,:,:,2)
			L(ix,:,:,3)=L(ix-2,:,:,3);L(ix-1,:,:,3)=L(ix-2,:,:,3)
			endif
		endif
		if (dir .eq. 2) then
			if (upper_lower .eq. 0) then
			L(:,1,:,1)=L(:,3,:,1);L(:,2,:,1)=L(:,3,:,1) 
			L(:,1,:,2)=L(:,3,:,2);L(:,2,:,2)=L(:,3,:,2)
			L(:,1,:,3)=L(:,3,:,3);L(:,2,:,3)=L(:,3,:,3)
			elseif (upper_lower .eq. 1) then
			L(:,jx,:,1)=L(:,jx-2,:,1);L(:,jx-1,:,1)=L(:,jx-2,:,1)
			L(:,jx,:,2)=L(:,jx-2,:,2);L(:,jx-1,:,2)=L(:,jx-2,:,2)
			L(:,jx,:,3)=L(:,jx-2,:,3);L(:,jx-1,:,3)=L(:,jx-2,:,3)
			endif
		endif
	endif
end subroutine visc_bc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine visc_ibc(L,V,dir,upper_lower)
	double precision,intent(inout)::L(ix,jx,kx,3)
	double precision,intent(in)::V(ix,jx,kx,3)
	integer, intent(in)::dir,upper_lower
	if(s_order.eq.4) then
		print*,'NOT PROGRAMMED'
		stop
	else
		if (dir .eq.1) then
			if (upper_lower .eq. 0) then
			L(1,:,:,1)=0.0d0; L(1,:,:,2)=L(2,:,:,2); L(1,:,:,3)=L(2,:,:,3)
			elseif (upper_lower .eq. 1) then
			L(ix,:,:,1)=0.0d0; L(ix,:,:,2)=L(ix-1,:,:,2); L(ix,:,:,3)=L(ix-1,:,:,3)
			endif
		endif
		if (dir .eq. 2) then
			if (upper_lower .eq. 0) then
			L(:,1,:,1)=L(:,2,:,1); L(:,1,:,2)=L(:,1,:,2); L(:,1,:,3)=L(:,1,:,3)
			elseif (upper_lower .eq. 1) then
			L(:,jx,:,1)=L(:,jx-1,:,1); L(:,jx,:,2)=L(:,jx-1,:,2); L(:,jx,:,3)=L(:,jx-1,:,3)
			endif
		endif
	endif
end subroutine visc_ibc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module Visc_rot
