module rad_cooling
!Radiative cooling module
!Parameters set in setting.txt:
!           flag_rad    - 0 (no cooling), 1 (energy based cooling), 2 (sech^2 type cooling) 
!           radrhoref   - reference density for the radiative time
!           rad_ts      - reference time scale for cooling

use parameters,only:pi
use globalvar,only:ix,jx,kx,nvar_h,nvar_m,ndim,flag_pip,flag_rad,edref,radrhoref,rad_ts,radlossfun,T0
use scheme_rot,only:get_Te_HD,get_Te_MHD

contains

subroutine source_rad_cooling(S_h,S_m,U_h,U_m)
	double precision,intent(inout)::S_h(ix,jx,kx,nvar_h),S_m(ix,jx,kx,nvar_m)
	double precision,intent(inout)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)
	double precision::rad_time(ix,jx,kx),Te_m(ix,jx,kx)
	double precision::deltatrad,Tdimindc,weight1,weight2,Tdim
	integer::Tdimindl,Tdimindu

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

!Cooling based of density
!        edref(:,:,:,1)=1.d0-dtanh((U_m(:,:,:,1)-1.25d0)*pi/0.1)**2
!        edref(:,:,:,1)=1.d0-dtanh(dlog10(U_m(:,:,:,1)/100.0)/2.d0*pi/0.5)**2  

!Set rad_time to 1. May want to change
        rad_time=1.0!(U_m(:,:,:,1)/radrhoref)**(-1.7)

!Cool based of temperature. Centered on a temperature of 0.1
        call get_Te_MHD(U_m,Te_m) !This assumes a factor of 1/2 in temperature so put 2*Te in next line
        edref(:,:,:,1)=(1.d0-dtanh(dlog10(2.0*Te_m(:,:,:)/0.1)/2.d0*pi/0.2)**2 )/rad_time/rad_ts
!        edref=max(edref-1.0e-5,0.d0)+1.0e-5 !This line doesn't work for some reason
!print*,maxval(edref(:,:,:,1)),minval(Te_m)
!Apply the energy sink of the form ro^2*Lambda
!		S_m(:,:,:,5)=S_m(:,:,:,5)-(U_m(:,:,:,5)-edref(:,:,:,1))/rad_time/rad_ts
		S_m(:,:,:,5)=S_m(:,:,:,5)-U_m(:,:,:,1)**2*(edref(:,:,:,1))

		if (flag_pip .eq. 1) then
            edref(:,:,:,2)=1.d0-dtanh((U_h(:,:,:,1)-1.25d0)*pi/0.1)**2
            rad_time=1.d0!((U_h(:,:,:,1))/radrhoref)**(-1.7)
!            S_h(:,:,:,5)=S_h(:,:,:,5)-(U_h(:,:,:,5)-edref(:,:,:,2))/rad_time/rad_ts
            S_h(:,:,:,5)=S_h(:,:,:,5)-(edref(:,:,:,2))/rad_time/rad_ts
		endif
		
	endif

	!Cooling using the radiative loss functions from Chianti
	if (flag_rad .eq. 3) then

        rad_time=1.0
		deltatrad=radlossfun(2,1)-radlossfun(1,1)

		!Get temperature
        call get_Te_MHD(U_m,Te_m) !This assumes a factor of 1/2 in temperature so put 2*Te in next line
		Te_m=2.d0*Te_m
!print*,radlossfun(2,1),radlossfun(1,1),radlossfun(2,1)-radlossfun(1,1),deltatrad
!stop
		!Interpolate the losses
		do i=1,ix
		do j=1,jx
		do k=1,kx
			Tdim=dlog10(Te_m(i,j,k)*T0)
			
			Tdimindl=floor((Tdim-3.d0)/deltatrad+1) !NEED TO INCLUDE THE DELTA
			Tdimindc=(Tdim-3.d0)
			Tdimindu=ceiling((Tdim-3.d0)/deltatrad+1)

!			print*,deltatrad,(Tdim-3.d0),Te_m(i,j,k)*T0,(Tdim-3.d0)/deltatrad
!			print*,dlog10(Te_m(i,j,k)*T0),radlossfun(floor((Tdim-3.d0)/deltatrad+1.0),1),radlossfun(ceiling((Tdim-3.d0)/deltatrad+1.0),1)

!			weight1=(radlossfun(Tdimindu,1)-Tdimindc)/(radlossfun(Tdimindu,1)-radlossfun(Tdimindl,1))
!			weight2=(Tdimindc-radlossfun(Tdimindl,1))/(radlossfun(Tdimindu,1)-radlossfun(Tdimindl,1))

			if (Tdim .LE. 3.5d0) then
				edref(i,j,k,1)=0.d0
			elseif (Tdimindl .eq. Tdimindu) then
				edref(i,j,k,1)=radlossfun(Tdimindl,2)
			else 
!				edref(i,j,k,1)=radlossfun(Tdimindl,2)*weight1+radlossfun(Tdimindu,2)*weight2
				edref(i,j,k,1)=radlossfun(Tdimindl,2)+(radlossfun(Tdimindu,2)-radlossfun(Tdimindl,2))*&
							(Tdim-radlossfun(Tdimindl,1))/(radlossfun(Tdimindu,1)-radlossfun(Tdimindl,1))
			endif
!			print*,Tdim,edref(i,j,k,1),rad_ts
		enddo
		enddo
		enddo
!stop
        edref(:,:,:,1)=edref(:,:,:,1)/rad_ts
        
		S_m(:,:,:,5)=S_m(:,:,:,5)-U_m(:,:,:,1)**2*(edref(:,:,:,1))

		
	endif

end subroutine source_rad_cooling

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initialize_radloss(flag_rad)
USE HDF5

	integer,intent(in)::flag_rad
	INTEGER :: ErrorFlag	
	INTEGER(HID_T) :: file_id
    INTEGER(HID_T) :: dset_id       ! Dataset identifier
    INTEGER(HID_T) :: space_id       ! Dataspace identifier
    INTEGER(HID_T) :: dtype_id       ! Dataspace identifier
    INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
    INTEGER(HSIZE_T), DIMENSION(2) :: max_dims
	Character(len=65),parameter::filename='lossfunc.h5'
	CHARACTER(LEN=65), PARAMETER :: dset1name = "temperature"  ! Dataset name
	CHARACTER(LEN=65), PARAMETER :: dset2name = "rad_loss"     ! Dataset name
	INTEGER::nelements,i

	if (flag_rad .ge. 1) then
		allocate(edref(ix,jx,kx,2))
	endif
	
	if (flag_rad .eq. 3) then
	print*,'Reading file'
		CALL h5open_f(ErrorFlag)
		CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, ErrorFlag)
		CALL h5dopen_f(file_id, dset1name, dset_id, ErrorFlag)
		CALL h5dget_space_f(dset_id, space_id,ErrorFlag)
		
		!Get dataspace dims
		CALL h5sget_simple_extent_dims_f(space_id,data_dims, max_dims, ErrorFlag)

		nelements = data_dims(1)

		!Allocate dimensions to dset_data for reading
		ALLOCATE(radlossfun(nelements,2))


		!Get data
		CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, radlossfun(:,1), data_dims, ErrorFlag)
		
		CALL h5dopen_f(file_id, dset2name, dset_id, ErrorFlag)
		CALL h5dget_space_f(dset_id, space_id,ErrorFlag)
		CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, radlossfun(:,2), data_dims, ErrorFlag)
		
		CALL h5close_f(ErrorFlag)

		radlossfun(:,2)=radlossfun(:,2)/maxval(radlossfun(:,2))
		print*,'file read'
!		do i=1,nelements
!			print*,radlossfun(i,1),radlossfun(i,2)!/maxval(radlossfun(:,2))
!		enddo

	endif

end subroutine initialize_radloss

subroutine set_ref_rad_ed(U_h,U_m)
	double precision,intent(in)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)	
	integer:: nrlf
	character*200::tmp

	edref(:,:,:,:)=0.0d0

!    if (flag_rad .eq. 2) then
!	    !Set the initial array to be zero
!    	edref(:,:,:,:)=0.0d0
!    else
!    	!Set the reference internal energy
!    	edref(:,:,:,1)=U_m(:,:,:,5)
!        if (flag_pip .eq. 1) edref(:,:,:,2)=U_h(:,:,:,5)
!    endif

    print*,'setting edref'

endsubroutine set_ref_rad_ed


end module rad_cooling
