module solver_rot
!====================================================================
! This is the solver module for the PIP code.
! Module for calling the different routines
! Author A.Hillier
! First version - 2013/04/23
! Second version - 2013/06/02
! Third version - 2013/06/04 NN
!====================================================================
  use globalvar,only:ix,jx,kx,nvar_h,nvar_m,ndim,ac,xi_n,dtout,&
       flag_mpi,my_rank,nout,time,flag_sch,flag_mhd,flag_pip,eta,&
       flag_divb,flag_grav,flag_amb,flag_artvis,flag_resi,flag_heat,&
       flag_visc,dt,dx,dy,dz,t_order,margin,dsc,dxc,dyc,dzc,ig,nt,x,flag_col, &
       flag_ir,no_advection,x,y,z,s_order,flag_hll,flag_time,hc_integ, &
       flag_hc_test,flag_cyl,gra,flag_damp
  ! This is the HLLD solver routines
  use HLL_rot,only:hll_fluxes_ideal_hd,hll_fluxes_ideal_mhd, &
       hllc_fluxes_ideal,hlld_fluxes_ideal
!  use MLW_rot        ! This is the MLW solver routines
!  use SLW_rot        ! This is the SLW solver routines
  use boundary_rot,only:PIPbnd    ! here all the boundary routines are called
!  use in_out      ! this is the library for the routines to save the output and also read files to give initial conditions for the simulations
 ! this contains all the routines need for each scheme to calculate cfl condition, converting between physcial and conserved quantities, calculate 
  use scheme_rot,only:set_heat_conductivity,&
       set_viscosity,&
       artvis,mhd_fluxes,hd_fluxes, &
       divb_cleaning_Dedner_source,divb_cleaning_Dedner_iter,derivative,vel_damp!,divb_flux_tmp
  use Res_rot,only:set_resistivity
  use PIP_rot,only:set_xin,set_IR,set_collisional,source_PIP
  use Gra_rot,only:set_gravity,source_gravity
  use HC_rot,only:HC,HC_STS,HC_SUBCYCLE
  use Util_rot,only:get_rotation

  implicit none
  integer n,order,ii

contains 
!------------------------------------------------------------------------
! add solver select here

  subroutine run_solver(U_h,U_m)
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)
    double precision ::U_h0(ix,jx,kx,nvar_h),U_m0(ix,jx,kx,nvar_m)
    double precision :: F_h(ix,jx,kx,nvar_h,3)
    double precision :: F_m(ix,jx,kx,nvar_m,3)
    double precision :: S_h(ix,jx,kx,nvar_h)
    double precision :: S_m(ix,jx,kx,nvar_m)
    double precision dt_sub,dt_coll_i
    double precision :: dt_coll(t_order)
    integer :: istep

    
    U_m0=U_m
    U_h0=U_h

    call set_dt_coll(dt_coll,t_order)
    
    do istep=1,t_order
       dt_coll_i=dt_coll(t_order-istep+1)
       ! set coefficients
       call set_coefficients(U_h,U_m,1)
       ! set numerical flux
       call set_flux(F_h,F_m,U_h,U_m)    
       ! set source term       
       call set_source(S_h,S_m,U_h0,U_m0,U_h,U_m,F_h,F_m,dt_coll_i)       
       ! time integral  
       call time_integral(F_h,F_m,S_h,S_m,U_h0,U_m0,U_h,U_m,dt_coll_i,istep)
       ! fix boundaries so that the substeps are updated, not the main variable
       call PIPbnd(U_h,U_m)
    enddo
    if (flag_damp.eq.1)  call vel_damp(U_h,U_m)

    
    call post_step(U_h,U_m,dt)

    call PIPbnd(U_h,U_m)


          
  end subroutine run_solver

  subroutine set_dt_coll(dt_coll,t_order)
    integer,intent(in)::t_order
    double precision,intent(out)::dt_coll(t_order)
    if (t_order.eq.4) then
       dt_coll(:)=dt*(/0.5d0,1.0d0/6.0d0,1.d0/12.d0,0.25d0/)   
    else if(t_order.eq.2) then
       dt_coll(:)=dt*(/0.5d0,0.5d0/)
    else
       dt_coll(:)=dt
    end if
  end subroutine set_dt_coll

  subroutine set_coefficients(U_h,U_m,flag)
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)
    integer,intent(in)::flag
    if(flag_resi.ge.2.or.(flag.eq.0.and.flag_resi.eq.1)) then
!       print*,'set_resistivity'
       call set_resistivity(U_h,U_m)       
    endif
    if(flag_grav.ge.2.or.(flag.eq.0.and.flag_grav.eq.1)) then
       call set_gravity(U_h,U_m)       
    endif
    if(flag_col.ge.2.or.(flag.eq.0.and.flag_col.eq.1)) then
       call set_collisional(U_h,U_m)       
    endif
    if(flag_amb.ge.2.or.(flag.eq.0.and.flag_amb.eq.1)) then
       call set_xin(U_h,U_m)       
    endif
    if(flag_IR.ge.2.or.(flag.eq.0.and.flag_IR.eq.1)) then
       call set_IR(U_h,U_m)       
    endif
  end subroutine set_coefficients

  !  subroutine set_flux(F_h,F_m,U_h,U_m,dt_sub,dt_coll_i)
  subroutine set_flux(F_h,F_m,U_h,U_m)  
    double precision,intent(inout)::F_h(ix,jx,kx,nvar_h,3),F_m(ix,jx,kx,nvar_m,3)
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)
    !    double precision,intent(in)::dt_sub,dt_coll_i
    integer dir        
    !Set flux for HD and MHD part----------------------------
    select case(flag_sch)
    case(0)!HLL-branch
       select case(flag_hll)
       case(0) !! HLL
          if(flag_mhd.eq.1) call hll_fluxes_ideal_mhd(F_m,U_m)
          if(flag_mhd.eq.0.or.flag_pip.eq.1) call hll_fluxes_ideal_hd(F_h,U_h)       
       case(1) !! Rusanov
          print*,'Rusanov flux has not been implemented yet...'
          stop
       case(2) !! HLLD
          if(flag_mhd.eq.1) call hlld_fluxes_ideal(F_m,U_m)
          if(flag_mhd.eq.0.or.flag_pip.eq.1) call hllc_fluxes_ideal(F_h,U_h)       
       end select
    case(1)!SLW 
       if(flag_mhd.eq.1) call mhd_fluxes(F_m,U_m)
       if(flag_mhd.eq.0.or.flag_pip.eq.1) call hd_fluxes(F_h,U_h)       
    end select
    !-------------------------------------------------------
    if(flag_mhd.eq.1) then
       if (no_advection.eq.1) then
!          F_m(:,:,:,6:8,1:3)=0.0d0
          F_m(:,:,:,:,1:3)=0.0d0
       endif
       if (flag_resi.ge.1.or.flag_artvis.eq.2) &
            call resistive_flux(F_m,U_m)

       if (flag_amb.ge.1)  &
            call ambipolar_flux(F_m,U_m)
    endif
  end subroutine set_flux


  subroutine set_source(S_h,S_m,U_h0,U_m0,U_h,U_m,F_h,F_m,dt_coll_i)  
    double precision,intent(inout):: S_h(ix,jx,kx,nvar_h),S_m(ix,jx,kx,nvar_m)
    double precision,intent(inout):: U_h0(ix,jx,kx,nvar_h),U_m0(ix,jx,kx,nvar_m)
    double precision,intent(inout):: U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)
    double precision,intent(inout):: dt_coll_i
    double precision,intent(inout) :: F_h(ix,jx,kx,nvar_h,3),F_m(ix,jx,kx,nvar_m,3)
    
    S_h=0.0d0
    S_m=0.0d0
    if(flag_grav.eq.1) call source_gravity(S_h,S_m,U_h,U_m)
    if(flag_pip.eq.1)  then
       call source_PIP(U_h0,U_m0,U_h,U_m,dt_coll_i,S_h,S_m)              
    endif
    if(flag_cyl.eq.1) call source_cyl(S_h,S_m,U_h,U_m,F_h,F_m) 
  end subroutine set_source
    
  subroutine time_integral(F_h,F_m,S_h,S_m,U_h0,U_m0, U_h,U_m,dt_coll_i,istep)
    double precision,intent(inout)::F_h(ix,jx,kx,nvar_h,3),F_m(ix,jx,kx,nvar_m,3)
    double precision,intent(inout)::S_h(ix,jx,kx,nvar_h),S_m(ix,jx,kx,nvar_m)
    double precision,intent(inout)::U_h0(ix,jx,kx,nvar_h),U_m0(ix,jx,kx,nvar_m)
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)
    double precision,intent(in) :: dt_coll_i
    integer,intent(in) :: istep
    double precision dt_sub,c0,c1
    integer i,j

    !! Heat conduction using super timestepping: first half step
    if(istep.eq.1 .and. flag_heat.eq.1 .and. hc_integ.eq.2) then
!       print*,'first step'
       select case(flag_mhd)
       case(0)
          call HC_STS(U_h,nvar_h,dt,0,0)
       case(1)
          call HC_STS(U_m,nvar_m,dt,1,0)
       end select
    endif

    !! MHD/HD step
    if (flag_hc_test.ne.1) then !! For heat conduction test
    select case(flag_time)
    case(0) !! standard RK
       dt_sub=dt/float(t_order-istep+1)       
       
       if(flag_mhd.eq.0 .or. flag_pip.eq.1) then
          call add_flux_rk(F_h,U_h0,U_h,nvar_h,dt_sub)
          call add_source(S_h,U_h,nvar_h,dt_sub)
       endif
       if(flag_mhd.eq.1) then 
          call add_flux_rk(F_m,U_m0,U_m,nvar_m,dt_sub)
          call add_source(S_m,U_m,nvar_m,dt_sub)
          if(flag_mhd.eq.1.and.flag_divb.eq.1) then
             call divb_cleaning_Dedner_source(dt_sub,U_m) !! original Dedner's 9-wave method
          endif                    
       endif
    case(1) !! optimal Strong Stability Preserving RK (SSPRK) Ref:Gottlieb+2009
       if(istep.eq.1) then
          c1 = 1.d0
          c0 = 1.d0
          dt_sub = dt
       else
          c1 = float(istep-1)/float(2*t_order-istep)          
          c0 = 1.d0-c1
          dt_sub = c1*dt
       endif
       
       if(flag_mhd.eq.0 .or. flag_pip.eq.1) then
          call add_flux_ssprk(F_h,U_h0,U_h,nvar_h,dt_sub,c0,c1,istep)
          call add_source(S_h,U_h,nvar_h,dt_sub)
       endif

       if(flag_mhd.eq.1) then 
          call add_flux_ssprk(F_m,U_m0,U_m,nvar_m,dt_sub,c0,c1,istep)
          call add_source(S_m,U_m,nvar_m,dt_sub)
          if(flag_mhd.eq.1.and.flag_divb.eq.1) then
             call divb_cleaning_Dedner_source(dt_sub,U_m) !! original Dedner's 9-wave method
          endif
       endif
    end select
    endif

    !! Heat conduction using super timestepping: second half step
    !! Heat conduction using subcycling
    if(istep.eq.t_order .and. flag_heat.eq.1 .and. hc_integ.ge.2) then
!       print*,'second step'
       select case(hc_integ)
       ! case(1)
       !    select case(flag_mhd)
       !    case(0)
       !       call HC(U_m,nvar_m,1)
       !    case(1)
       !       call HC(U_h,nvar_h,0)
       !    end select
       case(2)
          select case(flag_mhd)
          case(0)
             call HC_STS(U_h,nvar_h,dt,0,1)
          case(1)
             call HC_STS(U_m,nvar_m,dt,1,1)
          end select
       case(3)
          select case(flag_mhd)
          case(0)
             call HC_SUBCYCLE(U_h,nvar_h,dt,0)
          case(1)
             call HC_SUBCYCLE(U_m,nvar_m,dt,1)
          end select
       end select
    endif
    
  end subroutine time_integral
  
  subroutine post_step(U_h,U_m,dt_sub)  
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)
    double precision,intent(in)::dt_sub    
    integer dir,i,j,k

    select case(flag_sch)
    case(0) ! HLL
       !! nothing done       
    case(1) ! SLW
       call set_artvis(U_h,U_m,dt_sub)       

       !! Iterative divB cleaning should be applied AFTER artificial diffusion is applied
       if(flag_mhd.eq.1.and.flag_divb.eq.2) then
          call divb_cleaning_Dedner_iter(dt_sub,U_m)
       endif
    end select

    if(flag_heat.eq.1 .and. hc_integ.eq.1) then
       select case(flag_mhd)
       case(0)
          call HC(U_h,nvar_h,0)
       case(1)
          call HC(U_m,nvar_m,1)
       end select
    endif
  end subroutine post_step

  
  subroutine set_artvis(U_h,U_m,dt)
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)
    double precision, intent(in) :: dt
    integer dir

    if(flag_artvis.eq.1) then
       do dir=1,ndim
          if(flag_mhd.eq.1) then
             call artvis(dir,U_m,dt,nvar_m,1,nvar_m)
             call PIPbnd(U_h,U_m)
          endif
          if(flag_mhd.eq.0.or.flag_pip.eq.1) then
             call artvis(dir,U_h,dt,nvar_h,0,nvar_h)
             call PIPbnd(U_h,U_m)
          endif
       enddo
    else if(flag_artvis.eq.2) then
       do dir=1,ndim
          if(flag_mhd.eq.1) then
             call artvis(dir,U_m,dt,nvar_m,1,5)
             call PIPbnd(U_h,U_m)
          endif
          if(flag_mhd.eq.0.or.flag_pip.eq.1) then               
             call artvis(dir,U_h,dt,nvar_h,0,nvar_h)
             call PIPbnd(U_h,U_m)
          endif
       enddo
    end if

  end subroutine set_artvis
  
  subroutine add_flux_rk(F,U0,U,nvar,dt)  
    double precision,intent(inout)::U0(ix,jx,kx,nvar),U(ix,jx,kx,nvar)
    double precision,intent(inout)::F(ix,jx,kx,nvar,3)
    double precision,intent(inout)::dt
    integer,intent(in)::nvar
    integer i,j,dir!k,l 
    !for SLW centrer difference
    if(flag_sch.eq.1) then       
       dir=1
!       do i=margin(dir)+1,ig(dir)-margin(dir)
       do i=3,ig(dir)-2
          U(i,:,:,:)=U0(i,:,:,:) &
               - (dt/(dsc(i,dir)*12.0d0))* &
               (      -F(i+2,:,:,:,dir) & 
               + 8.0d0*F(i+1,:,:,:,dir) &
               - 8.0d0*F(i-1,:,:,:,dir) &
               +       F(i-2,:,:,:,dir))          
       enddo         
 
       if (ndim .ge. 2) then
          dir=2
!          do i=margin(dir)+1,ig(dir)-margin(dir)
          do i=3,ig(dir)-2
             U(:,i,:,:)=U(:,i,:,:) &
                  - (dt/(dsc(i,dir)*12.0d0))* &
                  (      -F(:,i+2,:,:,dir) & 
                  + 8.0d0*F(:,i+1,:,:,dir) &
                  - 8.0d0*F(:,i-1,:,:,dir) &
                  +       F(:,i-2,:,:,dir))          
          enddo

          if (ndim .ge. 3) then
             dir=3
!             do i=margin(dir)+1,ig(dir)-margin(dir)
             do i=3,ig(dir)-2
                U(:,:,i,:)=U(:,:,i,:) &
                     - (dt/(dsc(i,dir)*12.0d0))* &
                  (      -F(:,:,i+2,:,dir) & 
                  + 8.0d0*F(:,:,i+1,:,dir) &
                  - 8.0d0*F(:,:,i-1,:,dir) &
                  +       F(:,:,i-2,:,dir))          
             enddo
          endif
       endif
    else
       dir=1
!       do i=margin(dir)+1,ig(dir)-margin(dir)
       do i=2,ig(dir)
          U(i,:,:,:)=U0(i,:,:,:) &
               - (dt/dsc(i,dir))*(F(i,:,:,:,dir) -F(i-1,:,:,:,dir))
       enddo
       if (ndim .ge. 2) then
          dir=2
!          do i=margin(dir)+1,ig(dir)-margin(dir)
          do i=2,ig(dir)
             U(:,i,:,:)=U(:,i,:,:) &
                  - (dt/dsc(i,dir))*(F(:,i,:,:,dir) -F(:,i-1,:,:,dir))       
          enddo
          if (ndim .ge. 3) then
             dir=3
!             do i=margin(dir)+1,ig(dir)-margin(dir)
             do i=2,ig(dir)
                U(:,:,i,:)=U(:,:,i,:) &
                     - (dt/dsc(i,dir))*(F(:,:,i,:,dir) -F(:,:,i-1,:,dir))   
             enddo
          endif
       endif
    endif
  end subroutine add_flux_rk


  subroutine add_flux_ssprk(F,U0,U,nvar,dt,c0,c1,istep)  
    double precision,intent(inout)::U0(ix,jx,kx,nvar),U(ix,jx,kx,nvar)
    double precision,intent(inout)::F(ix,jx,kx,nvar,3)
    double precision,intent(in)::dt,c0,c1
    integer,intent(in)::nvar,istep
    integer i,j,dir!k,l 
    !for SLW centrer difference
    if(flag_sch.eq.1) then   !! SLW
       select case(istep)
       case(1)
          dir=1
          do i=3,ig(dir)-2
             U(i,:,:,:)=U0(i,:,:,:) &
                  - (dt/(dsc(i,dir)*12.0d0))* &
                  (      -F(i+2,:,:,:,dir) & 
                  + 8.0d0*F(i+1,:,:,:,dir) &
                  - 8.0d0*F(i-1,:,:,:,dir) &
                  +       F(i-2,:,:,:,dir))          
          enddo
       case(2:)
          dir=1
          do i=3,ig(dir)-2
             U(i,:,:,:)=c0*U0(i,:,:,:) + c1*U(i,:,:,:) &
                  - (dt/(dsc(i,dir)*12.0d0))* &
                  (      -F(i+2,:,:,:,dir) & 
                  + 8.0d0*F(i+1,:,:,:,dir) &
                  - 8.0d0*F(i-1,:,:,:,dir) &
                  +       F(i-2,:,:,:,dir))          
          enddo
       end select

       if (ndim .ge. 2) then
          dir=2
          do i=3,ig(dir)-2
             U(:,i,:,:)=U(:,i,:,:) &
                  - (dt/(dsc(i,dir)*12.0d0))* &
                  (      -F(:,i+2,:,:,dir) & 
                  + 8.0d0*F(:,i+1,:,:,dir) &
                  - 8.0d0*F(:,i-1,:,:,dir) &
                  +       F(:,i-2,:,:,dir))          
          enddo
          
          if (ndim .ge. 3) then
             dir=3
             do i=3,ig(dir)-2
                U(:,:,i,:)=U(:,:,i,:) &
                     - (dt/(dsc(i,dir)*12.0d0))* &
                     (      -F(:,:,i+2,:,dir) & 
                     + 8.0d0*F(:,:,i+1,:,dir) &
                     - 8.0d0*F(:,:,i-1,:,dir) &
                     +       F(:,:,i-2,:,dir))          
             enddo
          endif
       endif
    else !! HLL-branch
       select case(istep)
       case(1)
          dir=1
          do i=2,ig(dir)
             U(i,:,:,:)=U0(i,:,:,:) &
                  - (dt/dsc(i,dir))*(F(i,:,:,:,dir) -F(i-1,:,:,:,dir))
          enddo
       case(2:)
          dir=1
          do i=2,ig(dir)
             U(i,:,:,:)=c0*U0(i,:,:,:) + c1*U(i,:,:,:) &
                  - (dt/dsc(i,dir))*(F(i,:,:,:,dir) -F(i-1,:,:,:,dir))
          enddo          
       end select
       
       if (ndim .ge. 2) then
          dir=2
          do i=2,ig(dir)
             U(:,i,:,:)=U(:,i,:,:) &
                  - (dt/dsc(i,dir))*(F(:,i,:,:,dir) -F(:,i-1,:,:,dir))       
          enddo
          if (ndim .ge. 3) then
             dir=3
             do i=2,ig(dir)
                U(:,:,i,:)=U(:,:,i,:) &
                     - (dt/dsc(i,dir))*(F(:,:,i,:,dir) -F(:,:,i-1,:,dir))   
             enddo
          endif
       endif
    endif
  end subroutine add_flux_ssprk

  
  subroutine add_source(S,U,nvar,dt)
    double precision,intent(inout)::S(ix,jx,kx,nvar),U(ix,jx,kx,nvar)
    double precision,intent(inout)::dt    
    integer,intent(in)::nvar
    U=U+S*dt
  end subroutine add_source

subroutine resistive_flux(F_m,U_m)
  double precision,intent(inout):: F_m(ix,jx,kx,nvar_m,3),U_m(ix,jx,kx,nvar_m)
  double precision bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
!  double precision cx(ix,jx,kx),cy(ix,jx,kx),cz(ix,jx,kx),cc(ix,jx,kx)
  double precision cur(ix,jx,kx,3),cc(ix,jx,kx)
  double precision ex(ix,jx,kx),ey(ix,jx,kx),ez(ix,jx,kx),et_art(ix,jx,kx)
  double precision art_res,max_cc,deriv(ix,jx,kx),tmp,tmp2,bb(ix,jx,kx)
  double precision dgrid,limit,bmax
  integer i,j,k,n_art
  call get_rotation(ix,jx,kx,ndim,s_order,dxc,dyc,dzc,U_m(:,:,:,6:8),cur)
  if(flag_sch.eq.1) then
!     call get_current4th(cx,cy,cz,U_m)
     bx=U_m(:,:,:,6)
     by=U_m(:,:,:,7)
     bz=U_m(:,:,:,8)     
     if(flag_artvis.eq.2) then
!        art_res=0.00002d0
!        call set_resistivity(U_h,U_m)       
!        art_res=1.00d-3
        art_res=1.00d0
        et_art=art_res
        limit=0.05d0
!        et_art=0.0d0
!        cc=cx*cx+cy*cy+cz*cz
        cc=cur(:,:,:,1)*cur(:,:,:,1) &
             +cur(:,:,:,2)*cur(:,:,:,2) &
             +cur(:,:,:,3)*cur(:,:,:,3)
!        bb=max(bx*bx+by*by+bz*bz,1.0d-5)
        bmax=max(maxval(bx*bx+by*by+bz*bz),1.0d-5)
        max_cc=maxval(abs(cc))
        if(ndim.eq.2) then 
           dgrid=min(minval(dxc),minval(dyc))           
           do k=1,kx;do j=2,jx-1;do i=2,ix-1
              deriv(i,j,k)=abs(&
                   ((cc(i+1,j,k)-2.0d0*cc(i,j,k)+cc(i-1,j,k))/(2.0d0*dxc(i)**2) &
                   + (cc(i,j+1,k)-2.0d0*cc(i,j,k)+cc(i,j-1,k)) &
                   /(2.0d0*dyc(j)**2)) &
!                   /bb(i,j,k)*dgrid**4)              
                   /bmax*dgrid**5) 
!              if(deriv(i,j,k).gt.1) then
!              et_art(i,j,k)=art_res* &
!                   exp(-(abs(min(deriv(i,j,k),limit)-limit))/1.0d0)
!
           enddo;enddo;enddo
           et_art(3:ix-2,3:jx-2,1)=art_res* &
                min((0.125d0*(4.0d0*deriv(3:ix-2,3:jx-2,1)+&
                deriv(2:ix-3,3:jx-2,1)+deriv(4:ix-1,3:jx-2,1)+&
                deriv(3:ix-2,2:jx-3,1)+deriv(3:ix-2,4:jx-1,1)) &
                )**2,limit)/limit
        else if(ndim.eq.1) then
           dgrid=minval(dxc)
           do k=1,kx;do j=1,jx;do i=2,ix-1
              deriv(i,j,k)=abs(&
                   ((cc(i+1,j,k)-2.0d0*cc(i,j,k)+cc(i-1,j,k))/&
                   (2.0d0*dxc(i)**2)) &
!                   /bb(i,j,k)*dgrid**4)              
                   /bmax*dgrid**5) 
           enddo;enddo;enddo           
           et_art(3:ix-2,1,1)=art_res*min(0.25d0*( &
                deriv(2:ix-3,1,1)+deriv(3:ix-2,1,1)+deriv(3:ix-1,1,1)),&
                limit)/limit
        endif        
        eta=eta+et_art
     endif
!     ex=cx*eta
!     ey=cy*eta
!     ez=cz*eta
     ex=cur(:,:,:,1)*eta
     ey=cur(:,:,:,2)*eta
     ez=cur(:,:,:,3)*eta

     
     F_m(:,:,:,5,1)=F_m(:,:,:,5,1)+ey*bz-ez*by
     F_m(:,:,:,7,1)=F_m(:,:,:,7,1)-ez
     F_m(:,:,:,8,1)=F_m(:,:,:,8,1)+ey
     if (ndim.ge.2) then
        F_m(:,:,:,5,2)=F_m(:,:,:,5,2)+ez*bx-ex*bz
        F_m(:,:,:,6,2)=F_m(:,:,:,6,2)+ez
        F_m(:,:,:,8,2)=F_m(:,:,:,8,2)-ex
        if (ndim.ge.3) then
           F_m(:,:,:,5,3)=F_m(:,:,:,5,3)+ex*by-ey*bx
           F_m(:,:,:,6,3)=F_m(:,:,:,6,3)-ey
           F_m(:,:,:,7,3)=F_m(:,:,:,7,3)+ex
        endif
     endif     
  else     
!     call get_current(cx,cy,cz,U_m)
     bx=U_m(:,:,:,6)
     by=U_m(:,:,:,7)
     bz=U_m(:,:,:,8)
!     ex=cx*eta
!     ey=cy*eta
!     ez=cz*eta
     ex=cur(:,:,:,1)*eta
     ey=cur(:,:,:,2)*eta
     ez=cur(:,:,:,3)*eta
     
     F_m(1:ix-1,:,:,5,1)=F_m(1:ix-1,:,:,5,1)+ &
          0.25d0*((ey(2:ix,:,:)+ey(1:ix-1,:,:))*(bz(1:ix-1,:,:)+bz(2:ix,:,:)) &
          -(ez(1:ix-1,:,:)+ez(2:ix,:,:))*(by(1:ix-1,:,:)+by(2:ix,:,:)))
     F_m(1:ix-1,:,:,7,1)=F_m(1:ix-1,:,:,7,1)-&
          0.5d0*(ez(1:ix-1,:,:)+ez(2:ix,:,:))
     F_m(1:ix-1,:,:,8,1)=F_m(1:ix-1,:,:,8,1)+& 
          0.5d0*(ey(1:ix-1,:,:)+ey(2:ix,:,:))
     
     if (ndim.ge.2) then
        !        F_m(:,:,:,5,2)=F_m(:,:,:,5,2)+&!ez*bx-ex*bz
        F_m(:,1:jx-1,:,5,2)=F_m(:,1:jx-1,:,5,2)+&!ez*bx-ex*bz        
             0.25*((ez(:,2:jx,:)+ez(:,1:jx-1,:))*(bx(:,1:jx-1,:)+bx(:,2:jx,:)) &
             -(ex(:,1:jx-1,:)+ex(:,2:jx,:))*(bz(:,1:jx-1,:)+bz(:,2:jx,:)))   
        !        F_m(:,:,:,6,2)=F_m(:,:,:,6,2)+&!ez
        F_m(:,1:jx-1,:,6,2)=F_m(:,1:jx-1,:,6,2)+&!ez        
             0.5*(ez(:,1:jx-1,:)+ez(:,2:jx,:))
        !        F_m(:,:,:,8,2)=F_m(:,:,:,8,2)-&!ex
        F_m(:,1:jx-1,:,8,2)=F_m(:,1:jx-1,:,8,2)-&!ex        
             0.5*(ex(:,1:jx-1,:)+ex(:,2:jx,:))
        if (ndim.ge.3) then
           !           F_m(:,:,:,5,3)=F_m(:,:,:,5,3)+&!ex*by-ey*bx
           F_m(:,:,1:kx-1,5,3)=F_m(:,:,1:kx-1,5,3)+&!ex*by-ey*bx           
             0.25*((ex(:,:,2:kx)+ex(:,:,1:kx-1))*(by(:,:,1:kx-1)+by(:,:,2:kx)) &
             -(ey(:,:,1:kx-1)+ey(:,:,2:kx))*(bx(:,:,1:kx-1)+bx(:,:,2:kx)))   
           !           F_m(:,:,:,6,3)=F_m(:,:,:,6,3)-&!
           F_m(:,:,1:kx-1,6,3)=F_m(:,:,1:kx-1,6,3)-&!           
                0.5*(ey(:,:,2:kx)+ey(:,:,1:kx-1))
           !           F_m(:,:,:,7,3)=F_m(:,:,:,7,3)+&!
           F_m(:,:,1:kx-1,7,3)=F_m(:,:,1:kx-1,7,3)+&!           
                0.5*(ex(:,:,1:kx-1)+ex(:,:,2:kx))
        endif
     endif     


  end if

end subroutine resistive_flux


subroutine ambipolar_flux(F_m,U_m)
  double precision,intent(inout):: F_m(ix,jx,kx,nvar_m,3),U_m(ix,jx,kx,nvar_m)
  double precision bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
!  double precision cx(ix,jx,kx),cy(ix,jx,kx),cz(ix,jx,kx)
  double precision cur(ix,jx,kx,3),cc(ix,jx,kx)
  double precision ex(ix,jx,kx),ey(ix,jx,kx),ez(ix,jx,kx)
  double precision b2(ix,jx,kx),bc(ix,jx,kx),et_amb(ix,jx,kx)
  integer i  
  call get_rotation(ix,jx,kx,ndim,s_order,dxc,dyc,dzc,U_m(:,:,:,6:8),cur)
  bx=U_m(:,:,:,6)
  by=U_m(:,:,:,7)
  bz=U_m(:,:,:,8)
  b2=(bx**2+by**2+bz**2)
  et_amb=xi_n/(ac*(1.0d0-xi_n)*U_m(:,:,:,1)**2)
  bc=bx*cur(:,:,:,1)+by*cur(:,:,:,2)+bz*cur(:,:,:,3)
  ex=(b2*cur(:,:,:,1)-bc*bx)*et_amb
  ey=(b2*cur(:,:,:,2)-bc*by)*et_amb
  ez=(b2*cur(:,:,:,3)-bc*bz)*et_amb

  if(flag_sch.eq.1) then
     F_m(:,:,:,5,1)=F_m(:,:,:,5,1)+ey*bz-ez*by
     F_m(:,:,:,7,1)=F_m(:,:,:,7,1)-ez
     F_m(:,:,:,8,1)=F_m(:,:,:,8,1)+ey
     if (ndim.ge.2) then
        F_m(:,:,:,5,2)=F_m(:,:,:,5,2)+ez*bx-ex*bz
        F_m(:,:,:,6,2)=F_m(:,:,:,6,2)+ez
        F_m(:,:,:,8,2)=F_m(:,:,:,8,2)-ex
        if (ndim.ge.3) then
           F_m(:,:,:,5,3)=F_m(:,:,:,5,3)+ex*by-ey*bx
           F_m(:,:,:,6,3)=F_m(:,:,:,6,3)-ey
           F_m(:,:,:,7,3)=F_m(:,:,:,7,3)+ex
        endif
     endif
  else
     F_m(1:ix-1,:,:,5,1)=F_m(1:ix-1,:,:,5,1)+ &
          0.25d0*((ey(2:ix,:,:)+ey(1:ix-1,:,:))*(bz(1:ix-1,:,:)+bz(2:ix,:,:)) &
          -(ez(1:ix-1,:,:)+ez(2:ix,:,:))*(by(1:ix-1,:,:)+by(2:ix,:,:)))
     F_m(1:ix-1,:,:,7,1)=F_m(1:ix-1,:,:,7,1)-&
          0.5d0*(ez(1:ix-1,:,:)+ez(2:ix,:,:))
     F_m(1:ix-1,:,:,8,1)=F_m(1:ix-1,:,:,8,1)+& 
          0.5d0*(ey(1:ix-1,:,:)+ey(2:ix,:,:))     
     if (ndim.ge.2) then
        F_m(:,:,:,5,2)=F_m(:,:,:,5,2)+&!ez*bx-ex*bz
             0.25*((ez(:,2:jx,:)+ez(:,1:jx-1,:))*(bx(:,1:jx-1,:)+bx(:,2:jx,:)) &
             -(ex(:,1:jx-1,:)+ex(:,2:jx,:))*(bz(:,1:jx-1,:)+bz(:,2:jx,:)))   
        F_m(:,:,:,6,2)=F_m(:,:,:,6,2)+&!ez
             0.5*(ez(:,1:jx-1,:)+ez(:,2:jx,:))
        F_m(:,:,:,8,2)=F_m(:,:,:,8,2)-&!ex
             0.5*(ex(:,1:jx-1,:)+ex(:,2:jx,:))
        if (ndim.ge.3) then
           F_m(:,:,:,5,3)=F_m(:,:,:,5,3)+&!ex*by-ey*bx
             0.25*((ex(:,:,2:kx)+ex(:,:,1:kx-1))*(by(:,:,1:kx-1)+by(:,:,2:kx)) &
             -(ey(:,:,1:kx-1)+ey(:,:,2:kx))*(bx(:,:,1:kx-1)+bx(:,:,2:kx)))   
           F_m(:,:,:,6,3)=F_m(:,:,:,6,3)-&!
                0.5*(ey(:,:,2:kx)+ey(:,:,1:kx-1))
           F_m(:,:,:,7,3)=F_m(:,:,:,7,3)+&!
                0.5*(ex(:,:,1:kx-1)+ex(:,:,2:kx))
        endif
     endif          
  endif      
     !Old version of ambipolar diffusion
     ! This formula is stable, but numerical diffusion at the separatrix 
     ! is large 
!     bx=U_m(:,:,:,6)
!     by=U_m(:,:,:,7)
!     bz=U_m(:,:,:,8)
!     b2=(bx**2+by**2+bz**2)
!     bc=bx*cx+by*cy+bz*cz
!     ex=((b2*cx-bc*bx))/(ac*(1-xi_n)/xi_n*U_m(:,:,:,1)**2)
!     ey=((b2*cy-bc*by))/(ac*(1-xi_n)/xi_n*U_m(:,:,:,1)**2)
!     ez=((b2*cz-bc*bz))/(ac*(1-xi_n)/xi_n*U_m(:,:,:,1)**2)     
!     F_m(1:ix-1,:,:,5,1)=F_m(1:ix-1,:,:,5,1)+ &
!          0.25d0*((ey(2:ix,:,:)+ey(1:ix-1,:,:))*(bz(1:ix-1,:,:)+bz(2:ix,:,:)) &
!          -(ez(1:ix-1,:,:)+ez(2:ix,:,:))*(by(1:ix-1,:,:)+by(2:ix,:,:)))
!     F_m(1:ix-1,:,:,7,1)=F_m(1:ix-1,:,:,7,1)-&
!          0.5d0*(ez(1:ix-1,:,:)+ez(2:ix,:,:))
!     F_m(1:ix-1,:,:,8,1)=F_m(1:ix-1,:,:,8,1)+& 
!          0.5d0*(ey(1:ix-1,:,:)+ey(2:ix,:,:))     
!     if (ndim.ge.2) then
!        F_m(:,:,:,5,2)=F_m(:,:,:,5,2)+&!ez*bx-ex*bz
!             0.25*((ez(:,2:jx,:)+ez(:,1:jx-1,:))*(bx(:,1:jx-1,:)+bx(:,2:jx,:)) !&
!             -(ex(:,1:jx-1,:)+ex(:,2:jx,:))*(bz(:,1:jx-1,:)+bz(:,2:jx,:)))   
!        F_m(:,:,:,6,2)=F_m(:,:,:,6,2)+&!ez
!             0.5*(ez(:,1:jx-1,:)+ez(:,2:jx,:))
!        F_m(:,:,:,8,2)=F_m(:,:,:,8,2)-&!ex
!             0.5*(ex(:,1:jx-1,:)+ex(:,2:jx,:))
!        if (ndim.ge.3) then
!           F_m(:,:,:,5,3)=F_m(:,:,:,5,3)+&!ex*by-ey*bx
!             0.25*((ex(:,:,2:kx)+ex(:,:,1:kx-1))*(by(:,:,1:kx-1)+by(:,:,2:kx)) &
!             -(ey(:,:,1:kx-1)+ey(:,:,2:kx))*(bx(:,:,1:kx-1)+bx(:,:,2:kx)))   
!           F_m(:,:,:,6,3)=F_m(:,:,:,6,3)-&!
!                0.5*(ey(:,:,2:kx)+ey(:,:,1:kx-1))
!           F_m(:,:,:,7,3)=F_m(:,:,:,7,3)+&!
!                0.5*(ex(:,:,1:kx-1)+ex(:,:,2:kx))
!        endif
!     endif          
!  endif
end subroutine ambipolar_flux

!! NOTE: x=>r, y=>z, z=>phi
subroutine source_cyl(S_h,S_m,U_h,U_m,F_h,F_m)
  double precision,intent(inout):: S_h(ix,jx,kx,nvar_h),S_m(ix,jx,kx,nvar_m)
  double precision,intent(in):: U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)
  double precision,intent(in):: F_h(ix,jx,kx,nvar_h,3),F_m(ix,jx,kx,nvar_m,3)

  if(flag_mhd.eq.0.or.flag_pip.eq.1) then
     call source_cyl_sub(S_h,U_h,F_h,nvar_h,0)
  endif
  if(flag_mhd.eq.1) then
     call source_cyl_sub(S_m,U_m,F_m,nvar_m,1)     
  endif
  
end subroutine source_cyl


!! NOTE: x=>r, y=>z, z=>phi
subroutine source_cyl_sub(S,U,F,nvar,mhd)
  integer, intent(in) :: nvar,mhd
  double precision, intent(inout) :: S(ix,jx,kx,nvar)
  double precision, intent(in) :: U(ix,jx,kx,nvar)
  double precision, intent(in) :: F(ix,jx,kx,nvar,3)
  integer :: i
  
  select case(flag_sch)
  case(0) ! HLL-branch  
     !! continuity eq
     do i=2,ix
        S(i,:,:,1) = S(i,:,:,1) - 0.5d0*(F(i-1,:,:,1,1)+F(i,:,:,1,1))/x(i)
     enddo
     
     !! r (x)-component of the momentum eq
     do i=2,ix
        S(i,:,:,2) = S(i,:,:,2) &
             -  ( (U(i,:,:,2)**2-U(i,:,:,4)**2)/U(i,:,:,1) &
             + (U(i,:,:,8)**2-U(i,:,:,6)**2) ) / x(i)
     enddo
     
     !! z (y)-component of the momentum eq
     do i=2,ix
        S(i,:,:,3) = S(i,:,:,3) - 0.5d0*(F(i-1,:,:,3,1)+F(i,:,:,3,1))/x(i)
     enddo
     
     !! phi (z)-component of the momentum eq
     do i=2,ix
        S(i,:,:,4) = S(i,:,:,4) - (F(i-1,:,:,4,1)+F(i,:,:,4,1))/x(i)
     enddo
     
     !! enegy eq
     do i=2,ix
        S(i,:,:,5) = S(i,:,:,5) - 0.5d0*(F(i-1,:,:,5,1)+F(i,:,:,5,1))/x(i)
     enddo

     if(mhd.eq.1) then
     !! r and phi (x and z) component of the induction eq = 0
     !! z (y)-component of the induction eq
     do i=2,ix
        S(i,:,:,7) = S(i,:,:,7) - 0.5d0*(F(i-1,:,:,7,1)+F(i,:,:,7,1))/x(i)
     enddo
     
     !! psi
     if(flag_divb.eq.1) then
        do i=2,ix
           S(i,:,:,9) = S(i,:,:,9) - 0.5d0*(F(i-1,:,:,9,1)+F(i,:,:,9,1))/x(i)
        enddo
     endif
     endif
  case(1) ! SLW

     !! continuity eq
     do i=1,ix
        S(i,:,:,1) = S(i,:,:,1) - F(i,:,:,1,1)/x(i)
     enddo

     !! r (x)-component of the momentum eq
     do i=1,ix
        S(i,:,:,2) = S(i,:,:,2) - (F(i,:,:,2,1)-F(i,:,:,4,3))/x(i)
     enddo

     !! z (y)-component of the momentum eq
     do i=1,ix
        S(i,:,:,3) = S(i,:,:,3) - F(i,:,:,3,1)/x(i)
     enddo

     !! phi (z)-component of the momentum eq
     do i=1,ix
        S(i,:,:,4) = S(i,:,:,4) - 2.d0*F(i,:,:,4,1)/x(i)
     enddo

     !! energy eq
     do i=1,ix
        S(i,:,:,5) = S(i,:,:,5) - F(i,:,:,5,1)/x(i)
     enddo

     if(mhd.eq.1) then
     !! r and phi (x and z) component of the induction eq = 0
     !! z (y)-component of the induction eq
     do i=1,ix
        S(i,:,:,7) = S(i,:,:,7) - F(i,:,:,7,1)/x(i)
     enddo

     if(flag_divb.eq.1) then
        do i=1,ix
           S(i,:,:,9) = S(i,:,:,9) - F(i,:,:,9,1)/x(i)
        enddo
     endif
     endif
  end select

  
end subroutine source_cyl_sub

end module solver_rot

