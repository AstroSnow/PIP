module WENO_rot

  !====================================================================
  ! WENO Solver
  ! ONLY 1D AT THE MOMENT!!!!
  !====================================================================
  use globalvar,only:ndim,ix,jx,kx,nvar_h,nvar_m,s_order,gm &
       ,ro_lim,pr_lim,tiny,flag_divb,cmax,dxc
  use scheme_rot,only:hd_fluxes,mhd_fluxes
  use HLL_rot,only:v2u_1
  implicit none
  double precision,parameter::epsi=1.0d-6
  !  use Util_rot,only:hllc_flux,hlld_flux,MinMod_func
contains 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WENO_fluxes(F,U)
!Main controller for the WENO solver
    double precision,intent(inout)::U(ix,jx,kx,nvar_m)
    double precision,intent(out)::F(ix,jx,kx,nvar_m,3)
    double precision::F_P(ix,jx,kx,nvar_m,ndim),F_M(ix,jx,kx,nvar_m,ndim),F_temp(ix,jx,kx,nvar_m,ndim)
    integer::i
!print*,'WENO'
!    call WENO_flux_split(F_P,F_M,U)
    call WENO_flux_rusanov_MHD(F_P,F_M,U)

    F(:,:,:,:,:)=0.0d0
    F_temp(:,:,:,:,:)=0.0d0
    call WENO_flux_recon_5(F_temp,F_P,F_M,1,5,0,0) !x flux
    F(:,:,:,:,1)=F(:,:,:,:,1)+F_temp(:,:,:,:,1)
    if (ndim .ge. 2) then
    call WENO_flux_recon_5(F_temp,F_P,F_M,2,0,5,0) !y flux
    F(:,:,:,:,2)=F(:,:,:,:,2)+F_temp(:,:,:,:,2)
    else if (ndim .ge. 3) then
    call WENO_flux_recon_5(F_temp,F_P,F_M,3,0,0,5) !z flux
    F(:,:,:,:,3)=F(:,:,:,:,3)+F_temp(:,:,:,:,3)
    end if

end subroutine WENO_fluxes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine WENO_flux_split()
!Split the initial fluxes

!if (flag_weno_fs .eq. 1) then
!Local Lax-Friedrix (Rusanov) flux splitting
!    call WENO_flux_rusanov_HD(F_P,F_m,U)
!endif

!end WENO_flux_split
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WENO_flux_recon_5(flux,F_P,F_M,recdir,ixr,jxr,kxr)
!5th order WENO solver
! ONLY 1D AT THE MOMENT!
    double precision,intent(in)::F_P(ix,jx,kx,nvar_m),F_M(ix,jx,kx,nvar_m)
    integer,intent(in)::recdir,ixr,jxr,kxr
    double precision,intent(out)::flux(ix,jx,kx,nvar_m)
double precision::flux_temp(ix-ixr,jx-jxr,kx-kxr,nvar_m)
    double precision::vmm(ix-ixr,jx-jxr,kx-kxr,nvar_m),vm(ix-ixr,jx-jxr,kx-kxr,nvar_m),vo(ix-ixr,jx-jxr,kx-kxr,nvar_m)
    double precision::vpp(ix-ixr,jx-jxr,kx-kxr,nvar_m),vp(ix-ixr,jx-jxr,kx-kxr,nvar_m),B0n(ix-ixr,jx-jxr,kx-kxr,nvar_m)
    double precision::B1n(ix-ixr,jx-jxr,kx-kxr,nvar_m),B2n(ix-ixr,jx-jxr,kx-kxr,nvar_m)
    double precision::alpha0n(ix-ixr,jx-jxr,kx-kxr,nvar_m),alpha1n(ix-ixr,jx-jxr,kx-kxr,nvar_m),alpha2n(ix-ixr,jx-jxr,kx-kxr,nvar_m)
    double precision::alphasumn(ix-ixr,jx-jxr,kx-kxr,nvar_m)
    double precision::w0n(ix-ixr,jx-jxr,kx-kxr,nvar_m),w1n(ix-ixr,jx-jxr,kx-kxr,nvar_m),w2n(ix-ixr,jx-jxr,kx-kxr,nvar_m)
    double precision::umm(ix-ixr,jx-jxr,kx-kxr,nvar_m),um(ix-ixr,jx-jxr,kx-kxr,nvar_m),uo(ix-ixr,jx-jxr,kx-kxr,nvar_m)
    double precision::upp(ix-ixr,jx-jxr,kx-kxr,nvar_m),up(ix-ixr,jx-jxr,kx-kxr,nvar_m),B0p(ix-ixr,jx-jxr,kx-kxr,nvar_m)
    double precision::B1p(ix-ixr,jx-jxr,kx-kxr,nvar_m),B2p(ix-ixr,jx-jxr,kx-kxr,nvar_m)
    double precision::alpha0p(ix-ixr,jx-jxr,kx-kxr,nvar_m),alpha1p(ix-ixr,jx-jxr,kx-kxr,nvar_m),alpha2p(ix-ixr,jx-jxr,kx-kxr,nvar_m)
    double precision::alphasump(ix-ixr,jx-jxr,kx-kxr,nvar_m),phin(ix-ixr,jx-jxr,kx-kxr,nvar_m)
    double precision::w0p(ix-ixr,jx-jxr,kx-kxr,nvar_m),w1p(ix-ixr,jx-jxr,kx-kxr,nvar_m),w2p(ix-ixr,jx-jxr,kx-kxr,nvar_m)
integer::i

! Left Extrapolation
if (recdir .eq. 1) then
vmm = F_P(1:ix-5,:,:,:)
vm  = F_P(2:ix-4,:,:,:)
vo  = F_P(3:ix-3,:,:,:)
vp  = F_P(4:ix-2,:,:,:)
vpp = F_P(5:ix-1,:,:,:)
else if (recdir .eq. 2) then
vmm = F_P(:,1:jx-5,:,:)
vm  = F_P(:,2:jx-4,:,:)
vo  = F_P(:,3:jx-3,:,:)
vp  = F_P(:,4:jx-2,:,:)
vpp = F_P(:,5:jx-1,:,:)
else if (recdir .eq. 3) then
vmm = F_P(:,:,1:kx-5,:)
vm  = F_P(:,:,2:kx-4,:)
vo  = F_P(:,:,3:kx-3,:)
vp  = F_P(:,:,4:kx-2,:)
vpp = F_P(:,:,5:kx-1,:)
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Old version
! Smooth Indicators
!B0n = 13.d0/12.d0*(vmm-2.d0*vm+vo )**2 + 0.25d0*(vmm-4.d0*vm+3.d0*vo)**2
!B1n = 13.d0/12.d0*(vm -2.d0*vo +vp)**2 + 0.25d0*(vm-vp)**2
!B2n = 13.d0/12.d0*(vo -2.d0*vp+vpp)**2 + 0.25d0*(3.d0*vo-4.d0*vp+vpp)**2

!Alpha weights 
!alpha0n = 0.1d0/(max(epsi,B0n))**2
!alpha1n = 0.6d0/(max(epsi,B1n))**2
!alpha2n = 0.3d0/(max(epsi,B2n))**2
!alphasumn = alpha0n + alpha1n + alpha2n

!ENO stencils weigths
!w0n = alpha0n/alphasumn
!w1n = alpha1n/alphasumn
!w2n = alpha2n/alphasumn

!Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
!flux_temp(:,:,:,:) = w0n*(2.d0*vmm - 7.d0*vm + 11.d0*vo)/6.d0 &
!     + w1n*( -vm  + 5.d0*vo  + 2.d0*vp)/6.d0 &
!     + w2n*(2.d0*vo   + 5.d0*vp - vpp )/6.d0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Jiang & Wu flux reconstruction

!Moved to subroutine
!IS0=13.d0*(vmm-vm)**2+3.d0*(vmm-3.d0*vm)**2
!IS1=13.d0*(vm-vo)**2+3.d0*(vm+vo)**2
!IS2=13.d0*(vo-vp)**2+3.d0*(3.d0*vo-vp)**2

!alpha0n=1.d0/(tiny+IS0)**2
!alpha1n=6.d0/(tiny+IS1)**2
!alpha2n=3.d0/(tiny+IS2)**2
!alphasumn = alpha0n + alpha1n + alpha2n

!ENO stencils weigths
!w0n = alpha0n/alphasumn
!w2n = alpha2n/alphasumn

!Positive flux f_i+1/2
phin(:,:,:,:)=0.0d0
call WENO_weights(phin,vm-vmm,vo-vm,vp-vo,vpp-vp,ixr,jxr,kxr)
flux_temp(:,:,:,:)=(-vm+7.d0*vo+7.d0*vp-vpp)/12.d0-phin


! Right Extrapolation
if (recdir .eq. 1) then
umm = F_M(2:ix-4,:,:,:)
um  = F_M(3:ix-3,:,:,:)
uo  = F_M(4:ix-2,:,:,:)
up  = F_M(5:ix-1,:,:,:)
upp = F_M(6:ix  ,:,:,:)
elseif (recdir .eq. 2) then
umm = F_M(:,2:jx-4,:,:)
um  = F_M(:,3:jx-3,:,:)
uo  = F_M(:,4:jx-2,:,:)
up  = F_M(:,5:jx-1,:,:)
upp = F_M(:,6:jx  ,:,:)
elseif (recdir .eq. 3) then
umm = F_M(:,:,2:kx-4,:)
um  = F_M(:,:,3:kx-3,:)
uo  = F_M(:,:,4:kx-2,:)
up  = F_M(:,:,5:kx-1,:)
upp = F_M(:,:,6:kx  ,:)
endif

!print*,F_M(ix-10,1,1,:)
!Smooth Indicators
!B0p = 13.d0/12.d0*(umm-2.d0*um+uo )**2 + 0.25d0*(umm-4.d0*um+3.d0*uo)**2
!B1p = 13.d0/12.d0*(um -2.d0*uo +up)**2 + 0.25d0*(um-up)**2
!B2p = 13.d0/12.d0*(uo -2.d0*up+upp)**2 + 0.25d0*(3.d0*uo -4.d0*up+upp)**2

!Alpha weights 
!alpha0p = 0.1d0/(epsi + B0p)**2
!alpha1p = 0.6d0/(epsi + B1p)**2
!alpha2p = 0.3d0/(epsi + B2p)**2
!alphasump = alpha0p + alpha1p + alpha2p

!ENO stencils weigths
!w0p = alpha0p/alphasump;
!w1p = alpha1p/alphasump;
!w2p = alpha2p/alphasump;

!Numerical Flux at cell boundary, $u_{i+1/2}^{+}$;
!flux_temp(:,:,:,:) = flux_temp(:,:,:,:) + w0p*( -umm + 5.d0*um + 2.d0*uo  )/6.d0 &
!            + w1p*( 2.d0*um + 5.d0*uo  - up   )/6.d0 &
!            + w2p*(11.d0*uo  - 7.d0*up + 2.d0*upp)/6.d0

!Negative flux f_i-1/2
phin(:,:,:,:)=0.0d0
!call WENO_weights(phin,um-umm,uo-um,up-uo,upp-up,ixr,jxr,kxr)
!flux_temp(:,:,:,:)=flux_temp(:,:,:,:)+(-um+7.d0*uo+7.d0*up-upp)/12.d0-phin
call WENO_weights(phin,upp-up,up-uo,uo-um,um-umm,ixr,jxr,kxr)
flux_temp(:,:,:,:)=flux_temp(:,:,:,:)+(-umm+7.d0*um+7.d0*uo-up)/12.d0+phin



!Calculate residual
flux(:,:,:,:)=0.0d0
if (recdir .eq. 1) then
!    flux(3,:,:,:)=flux(3,:,:,:)-flux_temp(1,:,:,:) !Left cell
!    do i=1,ix-6
!        flux(i+2,:,:,:)=flux(i+2,:,:,:)+flux_temp(i,:,:,:)
!        flux(i+3,:,:,:)=flux(i+3,:,:,:)-flux_temp(i,:,:,:)
!    enddo
!    flux(ix-3,:,:,:)=flux(ix-3,:,:,:)+flux_temp(ix-6,:,:,:) !Right cell
!print*,size(flux(4:ix-3,1,1,1)),size(flux_temp(1:ix-6,1,1,1)),size(flux_temp(2:ix-5,1,1,1))
    flux(4:ix-3,:,:,:)=-flux_temp(1:ix-6,:,:,:)+flux_temp(2:ix-5,:,:,:)
else if (recdir .eq. 2) then
print*,'WRITE THIS IS VECTOR FORM'
stop
    flux(:,4,:,:)=flux(:,4,:,:)-flux_temp(:,1,:,:) !Left cell
    do i=2,jx-6
        flux(:,i+2,:,:)=flux(:,i+2,:,:)+flux_temp(:,i,:,:)
        flux(:,i+3,:,:)=flux(:,i+3,:,:)-flux_temp(:,i,:,:)
    enddo
    flux(:,jx-3,:,:)=flux(:,jx-3,:,:)+flux_temp(:,jx-6,:,:) !Right cell
else if (recdir .eq. 3) then

end if

!Temporary flux fix
!flux(ix-3:ix,:,:,:)=flux(ix-6:ix-3,:,:,:)
!flux(1:3,:,:,:)=flux(3:5,:,:,:)
!print*,size(flux(3:ix-3,1,1,1)),size(vmm(:,1,1,1))
!print*,'Break'
!do i=1,ix-5
!print*,i,flux(i,1,1,1,1)
!enddo
!print*,alpha2n(:,1,1,2)
!print*,alpha2n(:,1,1,5)
!stop
end subroutine WENO_flux_recon_5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WENO_flux_rusanov_HD(F_P,F_M,u)
!Split the initial fluxes
    double precision,intent(inout)::u(ix,jx,kx,nvar_m)
    double precision,intent(out)::F_P(ix,jx,kx,nvar_m,ndim),F_M(ix,jx,kx,nvar_m,ndim)
    double precision::flx(ix,jx,kx,nvar_m,ndim),pr(ix,jx,kx),a(ix,jx,kx)
!    double precision::atemp(5)
    integer::i
print*,'NEED TO DO THE HD FLUX PROPERLY'
stop
	! call v2f_1(v_l,flx_l,nvar,0)
	! call v2f_1(v_r,flx_r,nvar,0)    

!	pr = (gm-1.d0)*(u(:,:,:,5) &
!            - 0.5d0*u(:,:,:,1)*(u(:,:,:,2)**2 + u(:,:,:,3)**2 + u(:,:,:,4)**2))

!	flx(:,:,:,1)=u(:,:,:,2) 
!	flx(:,:,:,2)=pr +u(:,:,:,2)**2/u(:,:,:,1)
!	flx(:,:,:,3)=u(:,:,:,2)*u(:,:,:,3)/u(:,:,:,1) 
!	flx(:,:,:,4)=u(:,:,:,2)*u(:,:,:,4)/u(:,:,:,1)
!	flx(:,:,:,5)=u(:,:,:,2)*(u(:,:,:,5)+pr)/u(:,:,:,1) 

    call hd_fluxes(flx,U)

	a=dsqrt(gm*pr/u(:,:,:,1)) !Sound speed (might need an adjustment)
!    a=dsqrt(u(:,:,:,6)**2/u(:,:,:,1))
!    do i=3,ix-3
!        atemp(1)=abs(flx(i+1,1,1,1)-flx(i-1,1,1,1))/(dxc(i+1)+dxc(i-1))
!        atemp(2)=abs(flx(i+1,1,1,2)-flx(i-1,1,1,2))/(dxc(i+1)+dxc(i-1))
!        atemp(3)=abs(flx(i+1,1,1,3)-flx(i-1,1,1,3))/(dxc(i+1)+dxc(i-1))
!        atemp(4)=abs(flx(i+1,1,1,4)-flx(i-1,1,1,4))/(dxc(i+1)+dxc(i-1))
!        atemp(5)=abs(flx(i+1,1,1,5)-flx(i-1,1,1,5))/(dxc(i+1)+dxc(i-1))
!        a(i,:,:)=maxval(atemp)
!    enddo


    do i=1,5
	F_P(:,:,:,i,1)=0.5d0*(flx(:,:,:,i,1)+a(:,:,:)*u(:,:,:,i)) !can include a limiter on a
	F_M(:,:,:,i,1)=0.5d0*(flx(:,:,:,i,1)-a(:,:,:)*u(:,:,:,i))
    enddo
!print*,f_p(26:31,1,1,1),f_m(26:31,1,1,1),flx(26:31,1,1,1)
!stop
!print*,f_m(:,1,1,2)
!stop
end subroutine WENO_flux_rusanov_HD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WENO_flux_rusanov_MHD(F_P,F_M,u)
!Split the initial fluxes
    double precision,intent(inout)::u(ix,jx,kx,nvar_m)
    double precision,intent(out)::F_P(ix,jx,kx,nvar_m,ndim),F_M(ix,jx,kx,nvar_m,ndim)
    double precision::flx(ix,jx,kx,nvar_m,ndim),pr(ix,jx,kx),a(ix,jx,kx)
    double precision::b2(ix,jx,kx),pt(ix,jx,kx),vb(ix,jx,kx)
    double precision::roi(ix,jx,kx),cs2(ix,jx,kx),ca2(ix,jx,kx)
    double precision::cb2(ix,jx,kx),cfast2(ix,jx,kx)
!    double precision::atemp(5)
    integer::i 

!Gas Pressure
	pr = (gm-1.d0)*(u(:,:,:,5) &
            - 0.5d0*u(:,:,:,1)*(u(:,:,:,2)**2 + u(:,:,:,3)**2 + u(:,:,:,4)**2) &
            -0.5d0*( u(:,:,:,6)**2+u(:,:,:,7)**2+u(:,:,:,8)**2 ))


!B^2
!    b2 = u(:,:,:,6)**2+u(:,:,:,7)**2+u(:,:,:,8)**2
!Total pressure
!    pt = pr + 0.5d0*b2
!v dot B
!    vb = (u(:,:,:,2)*u(:,:,:,6)+u(:,:,:,3)*u(:,:,:,7)+u(:,:,:,4)*u(:,:,:,8))/u(:,:,:,1)

!THIS IS FOR x ONLY! NEED AN IF LOOP
!	flx(:,:,:,1)=u(:,:,:,2)  !rho*v_x
!	flx(:,:,:,2)=pr +u(:,:,:,2)**2/u(:,:,:,1) &
!        +0.5d0*b2 -u(:,:,:,6)**2 !pr + rho v_x^2 + 0.5*b^2 -B_x^2
!	flx(:,:,:,3)=u(:,:,:,2)*u(:,:,:,3)/u(:,:,:,1) -u(:,:,:,6)*u(:,:,:,7) !rho vx vy - Bx By
!	flx(:,:,:,4)=u(:,:,:,2)*u(:,:,:,4)/u(:,:,:,1) -u(:,:,:,6)*u(:,:,:,8) !rho vx vz - Bx Bz
!	flx(:,:,:,5)=u(:,:,:,2)*(u(:,:,:,5)+pt)/u(:,:,:,1) &
!        -vb*u(:,:,:,6) !vx(E+pt) - vb*B_x
!    flx(:,:,:,6)=0.0d0
!    flx(:,:,:,7)=u(:,:,:,2)*u(:,:,:,7)/u(:,:,:,1)-u(:,:,:,3)*u(:,:,:,6)/u(:,:,:,1)
!    flx(:,:,:,8)=u(:,:,:,2)*u(:,:,:,8)/u(:,:,:,1)-u(:,:,:,4)*u(:,:,:,6)/u(:,:,:,1)

    call mhd_fluxes(flx,U)
!print*,minval(flx)

!	a=dsqrt(gm*pr/u(:,:,:,1)) !Sound speed (SHOULD BE MAXIMUM OF WAVE SPEEDS?)
!a=dsqrt(u(:,:,:,6)**2/u(:,:,:,1))+dsqrt(gm*pr/u(:,:,:,1)) 
    !! fast mode wave speed
    roi = 1.d0/max(u(:,:,:,1),tiny)
    b2  = u(:,:,:,6)**2 + u(:,:,:,7)**2 + u(:,:,:,8)**2

    cs2 = gm*pr*roi !sound speed squared
    ca2 = u(:,:,:,6)**2*roi !Alfven speed squared
    cb2 = cs2 + b2*roi
    cfast2 = 0.5d0*( cb2 + sqrt( abs(cb2**2 - 4.d0*cs2*ca2) ) )
    a = sqrt(max(cfast2,tiny)) !fast wave speed

!print*,sqrt(cfast2(4:9,1,1))
!print*,sqrt(cs2(4:9,1,1))
!print*,sqrt(ca2(4:9,1,1))
!    do i=3,ix-3
!        atemp(1)=abs(flx(i+1,1,1,1)-flx(i-1,1,1,1))/(dxc(i+1)+dxc(i-1))
!        atemp(2)=abs(flx(i+1,1,1,2)-flx(i-1,1,1,2))/(dxc(i+1)+dxc(i-1))
!        atemp(3)=abs(flx(i+1,1,1,3)-flx(i-1,1,1,3))/(dxc(i+1)+dxc(i-1))
!        atemp(4)=abs(flx(i+1,1,1,4)-flx(i-1,1,1,4))/(dxc(i+1)+dxc(i-1))
!        atemp(5)=abs(flx(i+1,1,1,5)-flx(i-1,1,1,5))/(dxc(i+1)+dxc(i-1))
!        a(i,:,:)=maxval(atemp)
!    enddo


    do i=1,nvar_m
	F_P(:,:,:,i,1)=0.5d0*(flx(:,:,:,i,1)+a(:,:,:)*u(:,:,:,i)) !can include a limiter on a
	F_M(:,:,:,i,1)=0.5d0*(flx(:,:,:,i,1)-a(:,:,:)*u(:,:,:,i))
    enddo
!print*,f_p(26:31,1,1,1),f_m(26:31,1,1,1),flx(26:31,1,1,1)
!stop
!print*,f_m(:,1,1,2)
!stop
end subroutine WENO_flux_rusanov_MHD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WENO_weights(phin,a,b,c,d,ixr,jxr,kxr)
    double precision,intent(out)::phin(ix-ixr,jx-jxr,kx-kxr,nvar_m)
    integer,intent(in)::ixr,jxr,kxr
    double precision,intent(in)::a(ix-ixr,jx-jxr,kx-kxr,nvar_m),b(ix-ixr,jx-jxr,kx-kxr,nvar_m)
    double precision,intent(in)::c(ix-ixr,jx-jxr,kx-kxr,nvar_m),d(ix-ixr,jx-jxr,kx-kxr,nvar_m)
    double precision::IS0(ix-ixr,jx-jxr,kx-kxr,nvar_m),IS1(ix-ixr,jx-jxr,kx-kxr,nvar_m),IS2(ix-ixr,jx-jxr,kx-kxr,nvar_m)
    double precision::alpha0n(ix-ixr,jx-jxr,kx-kxr,nvar_m),alpha1n(ix-ixr,jx-jxr,kx-kxr,nvar_m),alpha2n(ix-ixr,jx-jxr,kx-kxr,nvar_m)
    double precision::alphasumn(ix-ixr,jx-jxr,kx-kxr,nvar_m), w0n(ix-ixr,jx-jxr,kx-kxr,nvar_m), w2n(ix-ixr,jx-jxr,kx-kxr,nvar_m)
    IS0=13.d0*(a-b)**2+3.d0*(a-3.d0*b)**2
    IS1=13.d0*(b-c)**2+3.d0*(b+c)**2
    IS2=13.d0*(c-d)**2+3.d0*(3.d0*c-d)**2

    alpha0n=1.d0/(tiny+IS0)**2
    alpha1n=6.d0/(tiny+IS1)**2
    alpha2n=3.d0/(tiny+IS2)**2
    alphasumn = alpha0n + alpha1n + alpha2n

    !ENO stencils weigths
    w0n = alpha0n/alphasumn
    w2n = alpha2n/alphasumn

    phin=w0n*(a-2.d0*b+c)/3.d0+(w2n-0.5d0)*(b-2.d0*c+d)/6.d0


end subroutine WENO_weights
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module WENO_rot
