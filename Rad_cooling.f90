module rad_cooling
!Radiative cooling module
!Parameters set in setting.txt:
!           flag_rad    - 0 (no cooling), 1 (energy based cooling), 2 (sech^2 type cooling) 
!           radrhoref   - reference density for the radiative time
!           rad_ts      - reference time scale for cooling

use parameters,only:pi
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

	!Sech-type profile for cooling of intermediate densities
	if (flag_rad .eq. 2) then

        edref(:,:,:,1)=1.d0-dtanh((U_m(:,:,:,1)-1.25d0)*pi/0.1)**2  

        rad_time=1.d0!(U_m(:,:,:,1)/radrhoref)**(-1.7)

!		S_m(:,:,:,5)=S_m(:,:,:,5)-(U_m(:,:,:,5)-edref(:,:,:,1))/rad_time/rad_ts
		S_m(:,:,:,5)=S_m(:,:,:,5)-(edref(:,:,:,1))/rad_time/rad_ts

		if (flag_pip .eq. 1) then
            edref(:,:,:,2)=1.d0-dtanh((U_h(:,:,:,1)-1.25d0)*pi/0.1)**2
            rad_time=1.d0!((U_h(:,:,:,1))/radrhoref)**(-1.7)
!            S_h(:,:,:,5)=S_h(:,:,:,5)-(U_h(:,:,:,5)-edref(:,:,:,2))/rad_time/rad_ts
            S_h(:,:,:,5)=S_h(:,:,:,5)-(edref(:,:,:,2))/rad_time/rad_ts
		endif
		
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
		

    if (flag_rad .eq. 2) then
	    !Set the initial array to be zero
    	edref(:,:,:,:)=0.0d0
    else
    	!Set the reference internal energy
    	edref(:,:,:,1)=U_m(:,:,:,5)
        if (flag_pip .eq. 1) edref(:,:,:,2)=U_h(:,:,:,5)
    endif

    print*,'setting edref'

endsubroutine set_ref_rad_ed


end module rad_cooling
