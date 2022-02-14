module rad_cooling

use globalvar,only:ix,jx,kx,nvar_h,nvar_m,ndim,flag_pip,flag_rad,edref,radrhoref,rad_ts


contains

subroutine source_rad_cooling(S_h,S_m,U_h,U_m)
	double precision,intent(inout)::S_h(ix,jx,kx,nvar_h),S_m(ix,jx,kx,nvar_m)
	double precision,intent(inout)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)
	double precision::rad_time(ix,jx,kx)

	!Newton Cooling
	!https://warwick.ac.uk/fac/sci/physics/research/cfsa/people/pastmembers/leakej/publications/4099.pdf
	if (flag_rad .eq. 1) then

        rad_time=(U_m(:,:,:,1)/radrhoref)**(-1.7)

		if (flag_pip .eq. 1) then
            rad_time=((U_m(:,:,:,1)+U_h(:,:,:,1))/radrhoref)**(-1.7)
            S_h(:,:,:,5)=S_h(:,:,:,5)-(U_h(:,:,:,5)-edref(:,:,:,2))/rad_time/rad_ts
		endif

		S_m(:,:,:,5)=S_m(:,:,:,5)-(U_m(:,:,:,5)-edref(:,:,:,1))/rad_time/rad_ts
		
	endif

end subroutine source_rad_cooling

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initialize_radloss(flag_rad)
	integer,intent(in)::flag_rad

	if (flag_rad .ge. 1) then
		allocate(edref(ix,jx,kx,2))
	endif

end subroutine initialize_radloss

subroutine set_ref_rad_ed(U_h,U_m)
	double precision,intent(in)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)	
		
	!Set the reference internal energy
	edref(:,:,:,1)=U_m(:,:,:,5)
    if (flag_pip .eq. 1) edref(:,:,:,2)=U_h(:,:,:,5)

    print*,'setting edref'

endsubroutine set_ref_rad_ed


end module rad_cooling
