module scheme_rot
!====================================================================
! This is the scheme_routines module for the PIP code.
! Module for subroutines for schemes 
! contents
! pv2cq  ---- convert physical variables to conserved quantity
! cq2pv  ---- convert conserved quantity to physical variables 
! cfl    ---- set CFL condition
! First version - 2013/04/26 NN
! Second version - 2013/06/03 changed cfl to make it useable for 2d AH
! Third  version - 2013/06/27 - rountines added that are reuired to be used in the substeps of the SLW routine
!====================================================================

  use globalvar,only:ix,jx,kx,nvar_h,nvar_m,eta,ac,safety,&
       pip_imp_factor,dt,db_clean,flag_b_stg,dx,dy,dz,theta,ndim,flag_divb,&
       flag_amb,flag_mhd,flag_pip,flag_amb,flag_mpi,flag_resi,margin,gm,&
       flag_bnd,xi_n,flag_pip_imp,nt,eta_0,gra,scl_height,s_order,flag_sch,&
       x,y,z,col,n_fraction,flag_ir,gm_rec,gm_ion,gm_rec_rad,gm_ion_rad,t_ir,mpi_pos,my_rank,&
       debug_parameter,ro_lim,pr_lim,tiny,cmax,dsc,b_cr,damp_time,flag_damp,&
       oldke_damp,flag_rad,ion_pot,flag_IR_type,nexcite,colrat,radrat,gm_rec_ref
  use MPI_rot,only:mpi_double_interface
  use Boundary_rot,only:bnd_divb
  !use PIP_rot,only:hydrogen_excitation_timestep
  implicit none
  integer i,j,k,ierr

contains 
!------------------------------------------------------------------------
  subroutine pv2cq_hd(de,vx,vy,vz,pr,U_h)
    double precision,intent(inout):: de(ix,jx,kx),pr(ix,jx,kx)
    double precision,intent(inout)::vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h)
    u_h = 0.d0 ! initialization
    u_h(:,:,:,1)=de
    u_h(:,:,:,2)=de*vx
    u_h(:,:,:,3)=de*vy
    u_h(:,:,:,4)=de*vz
    u_h(:,:,:,5)=pr/(gm-1.0d0)+0.5d0*de*(vx**2+vy**2+vz**2)    
  end subroutine pv2cq_hd

  !Convert physical variables to conserved quantity for MHD
  subroutine pv2cq_mhd (de,vx,vy,vz,pr,bx,by,bz,U_m)
    double precision,intent(inout):: de(ix,jx,kx),vx(ix,jx,kx)
    double precision,intent(inout)::vy(ix,jx,kx),vz(ix,jx,kx)
    double precision,intent(inout):: pr(ix,jx,kx),bx(ix+flag_b_stg,jx,kx)
    double precision,intent(inout):: by(ix,jx+flag_b_stg,kx),bz(ix,jx,kx+flag_b_stg)
    double precision,intent(inout)::U_m(ix,jx,kx,nvar_m)
!    u_m = 0.d0 ! initialization
    if(flag_b_stg.eq.1) then
       u_m(:,:,:,1)=de
       u_m(:,:,:,2)=de*vx
       u_m(:,:,:,3)=de*vy
       u_m(:,:,:,4)=de*vz
       u_m(:,:,:,6)=0.5d0*(bx(1:ix,:,:)+bx(2:ix+1,:,:))
       if (ndim.ge.2) then 
          u_m(:,:,:,7)=0.5d0*(by(:,1:jx,:)+by(:,2:jx+1,:))       
       else
          u_m(:,:,:,7)=by(:,1:jx,:)
       endif
       if (ndim.ge.3) then 
          u_m(:,:,:,8)=0.5d0*(bz(:,:,1:kx)+bz(:,:,2:kx+1))       
       else
          u_m(:,:,:,8)=bz(:,:,1:kx)
       endif
       u_m(:,:,:,5)=pr/(gm-1.0d0)+0.5d0*de*(vx**2+vy**2+vz**2) &
            +0.5d0*(u_m(:,:,:,6)**2+u_m(:,:,:,7)**2+u_m(:,:,:,8)**2)
    else
       u_m(:,:,:,1)=de
       u_m(:,:,:,2)=de*vx
       u_m(:,:,:,3)=de*vy
       u_m(:,:,:,4)=de*vz
       u_m(:,:,:,5)=pr/(gm-1.0d0)+0.5d0*de*(vx**2+vy**2+vz**2) &
            +0.5d0*(bx**2+by**2+bz**2)
       u_m(:,:,:,6)=bx
       u_m(:,:,:,7)=by
       u_m(:,:,:,8)=bz       
    endif
  end subroutine pv2cq_mhd

!------------------------------------------------------------------------
!------------------------------------------------------------------------
  subroutine get_Te_HD(U,Te)
    double precision,intent(inout)::U(ix,jx,kx,nvar_h)
    double precision,intent(inout)::Te(ix,jx,kx)
    double precision pr(ix,jx,kx)
    call get_Pr_HD(U,pr)
    Te=gm*(pr/U(:,:,:,1))    
  end subroutine get_Te_HD
  subroutine get_Pr_HD(U,Pr)
    double precision,intent(inout)::U(ix,jx,kx,nvar_m)
    double precision,intent(inout)::Pr(ix,jx,kx)
    pr=(gm-1.0d0)*(U(:,:,:,5)&
         -0.5d0*(U(:,:,:,2)*U(:,:,:,2) &
         +U(:,:,:,3)*U(:,:,:,3) &
         +U(:,:,:,4)*U(:,:,:,4))/U(:,:,:,1))
!HD PRESSURE FIXES
   if(minval(pr).le.pr_lim) pr=max(pr,pr_lim)  
    U(:,:,:,5)=pr/(gm-1.0d0)+0.5d0*((U(:,:,:,2)*U(:,:,:,2)&
         +U(:,:,:,3)*U(:,:,:,3)&
         +U(:,:,:,4)*U(:,:,:,4))/U(:,:,:,1))
  end subroutine get_Pr_HD

  subroutine get_Te_MHD(U,Te)
    double precision,intent(inout)::U(ix,jx,kx,nvar_m)
    double precision,intent(inout)::Te(ix,jx,kx)
    double precision pr(ix,jx,kx)
    call get_Pr_MHD(U,pr)
    Te=0.5d0*gm*(pr/U(:,:,:,1))    
  end subroutine get_Te_MHD
  subroutine get_Pr_MHD(U,pr)
    double precision,intent(inout)::U(ix,jx,kx,nvar_m)
    double precision,intent(inout)::Pr(ix,jx,kx)
    pr=(gm-1.0d0)*(u(:,:,:,5)&
         -0.5d0*((U(:,:,:,2)*U(:,:,:,2)&
         +U(:,:,:,3)*U(:,:,:,3)&
         +U(:,:,:,4)*U(:,:,:,4))/U(:,:,:,1)&
         +U(:,:,:,6)*U(:,:,:,6)&
         +U(:,:,:,7)*U(:,:,:,7)&
         +U(:,:,:,8)*U(:,:,:,8)))
   if(minval(pr).le.pr_lim) pr=max(pr,pr_lim)  
    u(:,:,:,5)=pr/(gm-1.0d0)+0.5d0*((U(:,:,:,2)*U(:,:,:,2)&
         +U(:,:,:,3)*U(:,:,:,3)&
         +U(:,:,:,4)*U(:,:,:,4))/U(:,:,:,1)&
         +U(:,:,:,6)*U(:,:,:,6)&
         +U(:,:,:,7)*U(:,:,:,7)&
         +U(:,:,:,8)*U(:,:,:,8)) 
  end subroutine get_Pr_MHD

  
  !Convert physical variables to conserved quantity for HD
  subroutine cq2pv_hd(de,vx,vy,vz,pr,U_h)
    double precision,intent(inout):: de(ix,jx,kx),pr(ix,jx,kx)
    double precision,intent(inout)::vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h)
    de=u_h(:,:,:,1)
    vx=u_h(:,:,:,2)/de
    vy=u_h(:,:,:,3)/de
    vz=u_h(:,:,:,4)/de
    pr=(gm-1.0d0)*(u_h(:,:,:,5)-0.5d0*de*(vx**2+vy**2+vz**2))
!!!!!!FIXES FOR HD DENSITY, PRESSURE AND ENERGY !!!!!!!  
    if(minval(pr).le.pr_lim) pr=max(pr,pr_lim)
    if(minval(de).le.ro_lim) de=max(de,ro_lim)
    u_h(:,:,:,1)=de
    u_h(:,:,:,2)=vx*de
    u_h(:,:,:,3)=vy*de
    u_h(:,:,:,4)=vz*de
    u_h(:,:,:,5)=pr/(gm-1.0d0)+0.5d0*de*(vx**2+vy**2+vz**2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  end subroutine cq2pv_hd

  !Convert physical variables to conserved quantity for MHD
  subroutine cq2pv_mhd (de,vx,vy,vz,pr,bx,by,bz,U_m)
    double precision,intent(inout):: de(ix,jx,kx),vx(ix,jx,kx)
    double precision,intent(inout)::vy(ix,jx,kx),vz(ix,jx,kx)
    double precision,intent(inout):: pr(ix,jx,kx)
    double precision,intent(inout)::bx(ix,jx,kx)
    double precision,intent(inout):: by(ix,jx,kx),bz(ix,jx,kx)
    double precision,intent(inout)::U_m(ix,jx,kx,nvar_m)

   
    de=u_m(:,:,:,1)
    vx=u_m(:,:,:,2)/de
    vy=u_m(:,:,:,3)/de
    vz=u_m(:,:,:,4)/de
    bx=u_m(:,:,:,6)
    by=u_m(:,:,:,7)
    bz=u_m(:,:,:,8)
    pr=(gm-1.0d0)*(u_m(:,:,:,5)-0.5d0*de*(vx**2+vy**2+vz**2)  &
         -0.5d0*(bx**2+by**2+bz**2))
    if(minval(pr).le.pr_lim) pr=max(pr,pr_lim)
    if(minval(de).le.ro_lim) de=max(de,ro_lim)
    u_m(:,:,:,1)=de
    u_m(:,:,:,2)=vx*de
    u_m(:,:,:,3)=vy*de
    u_m(:,:,:,4)=vz*de
    u_m(:,:,:,5)=pr/(gm-1.0d0)+0.5d0*de*(vx**2+vy**2+vz**2)  &
         +0.5d0*(bx**2+by**2+bz**2)
  end subroutine cq2pv_mhd


  !set time step with CFL condition : return dt
  subroutine cfl(U_h,U_m,dttemp)
    double precision,intent(inout)::U_m(ix,jx,kx,nvar_m),U_h(ix,jx,kx,nvar_h)
    double precision,intent(out)::dttemp
    if(flag_pip.eq.1) then
       call cfl_mhd(U_m,dttemp)
       call cfl_hd(U_h,dttemp)
       call cfl_pn_col(U_m,U_h,dttemp)
       if(flag_ir.ge.1) call cfl_pip_ir(U_m,U_h,dttemp)
       if(flag_ir.eq.4) call cfl_pip_ir_nexcite(U_m,U_h,dttemp)
    else
       if(flag_mhd.eq.1) then
          call cfl_mhd(U_m,dttemp)
       else
          call cfl_hd(U_h,dttemp)
       endif
    endif
    if(flag_resi.ge.1) call cfl_resi
    if(flag_mpi.eq.1) then
       dttemp=mpi_double_interface(dttemp,1)
    endif

    !! for divB cleaning
    if(flag_divb.ge.1) then
!       cmax = safety*minval(dsc)/dt
       !       cmax = min(minval(dx),minval(dy),minval(dz))/dt
       cmax = safety*min(minval(dx),minval(dy),minval(dz))/dt       
       if(flag_mpi.eq.1) &
            cmax = mpi_double_interface(cmax,3)
    endif
    
  end subroutine cfl
  !set time step for HD with CFL condition 
  subroutine cfl_hd(U_h,dttemp)
  double precision,intent(inout)::dttemp
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h)
    double precision dt_min,cs2,v2,gmin,cs,vabs
    double precision de(ix,jx,kx),pr(ix,jx,kx)
    double precision vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
    call cq2pv_hd(de,vx,vy,vz,pr,U_h)

    if (flag_pip.eq.1) then 
       dt_min=dttemp
    else
       dt_min=10.0    
    endif

    do k=margin(3)+1,kx-margin(3)
       do j=margin(2)+1,jx-margin(2) 
          do i=margin(1)+1,ix-margin(1)
             gmin=min(dx(i),dy(j),dz(k))
!             cs2=gm*pr(i,j,k)/de(i,j,k) 
!             v2=vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2
!             dt_min=min(dt_min,safety*gmin/sqrt(v2+cs2))
             cs=sqrt(gm*pr(i,j,k)/de(i,j,k) )
             vabs=sqrt(max(vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2,tiny))
             dt_min=min(dt_min,safety*gmin/(vabs+cs))
          enddo
       enddo
    enddo
    dttemp=dt_min
  end subroutine cfl_hd
  
  subroutine cfl_mhd(U_m,dttemp)
    double precision,intent(inout)::dttemp
    double precision,intent(inout)::U_m(ix,jx,kx,nvar_m)
    double precision dt_min,gmin
    double precision de(ix,jx,kx),pr(ix,jx,kx)
    double precision vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
    double precision bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
    double precision etmin,maxet
    double precision cs2,valf2,vabs

    call cq2pv_mhd(de,vx,vy,vz,pr,bx,by,bz,U_m)
    dt_min=10.0

    ! !set dt using cfl condition
    ! !This is CFL for uniform grid 
    ! !If you use this non-uniform grid this dt may be under estimate of dt
    ! gmin=min(minval(dx),minval(dy),minval(dz))    
    ! dt_min=min(dt_min,safety*gmin/maxval( &
    !      sqrt(vx*vx+vy*vy+vz*vz)+ &
    !      sqrt(gm*pr/de) + &
    !      sqrt((bx*bx+by*by+bz*bz)/de)))

    do k=margin(3)+1,kx-margin(3)
       do j=margin(2)+1,jx-margin(2) 
          do i=margin(1)+1,ix-margin(1)
             gmin=min(dx(i),dy(j),dz(k))
             cs2=gm*pr(i,j,k)/max(de(i,j,k),tiny)
             valf2=(bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2)/max(de(i,j,k),tiny)
             vabs=sqrt(vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2)
             dt_min=min(dt_min,safety*gmin/(vabs+sqrt(cs2+valf2))) !! fast mode speed+fluid speed
	     !print*,'dt: ',dt_min
          enddo
       enddo
    enddo
    
    !! This part needs to be modified for non-uniform grids    
    dttemp=dt_min
    if(flag_amb.eq.1.and.flag_pip.ne.1) then 
       etmin=1.0d-8
       maxet=0.0d0
       dttemp=min(dttemp,safety*gmin**2/ &
            maxval(xi_n*(bx*bx+by*by+bz*bz)/&
            (ac*de*de*(1.0-xi_n))))
    endif
  end subroutine cfl_mhd
  
  subroutine cfl_resi
    double precision etmin,gmin
    etmin=1.0d-8
    do k=margin(3)+1,kx-margin(3) 
       do j=margin(2)+1,jx-margin(2) 
          do i=margin(1)+1,ix-margin(1)
             gmin=min(dx(i),dy(j),dz(k))
             dt=min(dt,safety*gmin**2/max(eta(i,j,k),etmin))
          enddo
       enddo
    enddo
  end subroutine cfl_resi


  subroutine cfl_pn_col(U_m,U_h,dttemp)
    double precision,intent(inout)::U_m(ix,jx,kx,nvar_m),U_h(ix,jx,kx,nvar_h)
      double precision,intent(inout)::dttemp
    if(flag_pip_imp.eq.1) then
       dttemp=min(dttemp,pip_imp_factor/max(maxval(u_m(:,:,:,1)*ac(:,:,:)),maxval(u_h(:,:,:,1)*ac(:,:,:))))
    else
       dttemp=min(dttemp,min(safety,0.3d0)/max(maxval(u_m(:,:,:,1)*ac(:,:,:)),maxval(u_h(:,:,:,1)*ac(:,:,:))))
    endif
   end subroutine cfl_pn_col

   subroutine cfl_pip_ir(U_m,U_h,dttemp)
  double precision,intent(inout)::dttemp
     double precision,intent(inout)::U_m(ix,jx,kx,nvar_m),U_h(ix,jx,kx,nvar_h)
     double precision::pr(ix,jx,kx)
     dttemp=min(dttemp,min(safety,0.1d0)/max(maxval(gm_rec),maxval(gm_ion),1.0d-5))   !,maxval(gm_rec/U_m(:,:,:,1))+maxval(gm_ion/U_h(:,:,:,1)) 
     if (flag_IR_type .eq. 0) then
         call get_Pr_MHD(U_m,pr)
     	dttemp=min(dttemp,safety/max(maxval(ion_pot/(pr/(gm-1.d0))),1.0d-5)) 
     endif
     if (flag_rad .eq. 3) then 
     	dttemp=min(dttemp,min(safety,0.1d0)/max(maxval(gm_rec_rad),maxval(gm_ion_rad),1.0d-5)) 
     endif
   end subroutine cfl_pip_ir

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cfl_pip_ir_nexcite(U_m,U_h,dttemp)
     double precision,intent(inout)::dttemp
     double precision,intent(inout)::U_m(ix,jx,kx,nvar_m),U_h(ix,jx,kx,nvar_h)
     double precision::dneut(ix,jx,kx,5)
     call hydrogen_excitation_timestep(U_m,U_h,dneut)
     dttemp=min(dttemp,min(safety,0.3d0)/max(maxval(dneut),1.0d-5))
   end subroutine cfl_pip_ir_nexcite
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hd_fluxes(F_h,U_h)
    double precision,intent(inout):: F_h(ix,jx,kx,nvar_h,3)
    double precision,intent(inout):: U_h(ix,jx,kx,nvar_h)
    double precision :: de(ix,jx,kx),vx(ix,jx,kx)
    double precision ::vy(ix,jx,kx),vz(ix,jx,kx)
    double precision :: pr(ix,jx,kx)  
    integer, parameter :: ro=1
    integer, parameter :: mx=2
    integer, parameter :: my=3
    integer, parameter :: mz=4
    integer, parameter :: en=5    
    call cq2pv_hd(de,vx,vy,vz,pr,U_h)
    !---------------------------------------------------------------------  
    ! ideal HD part
    !---------------------------------------------------------------------  
    ! mass
    F_h(:,:,:,ro,1) = U_h(:,:,:,mx)
    F_h(:,:,:,ro,2) = U_h(:,:,:,my)
    F_h(:,:,:,ro,3) = U_h(:,:,:,mz)
    
    ! x-momentum
    F_h(:,:,:,mx,1) = U_h(:,:,:,mx)**2/U_h(:,:,:,ro) &
         + pr(:,:,:) 
    F_h(:,:,:,mx,2) = vy(:,:,:)*U_h(:,:,:,mx) 
    F_h(:,:,:,mx,3) = vz(:,:,:)*U_h(:,:,:,mx) 
    
    ! y-momentum
    F_h(:,:,:,my,1) = vx(:,:,:)*U_h(:,:,:,my)
    F_h(:,:,:,my,2) = U_h(:,:,:,my)**2/U_h(:,:,:,ro) &
         + pr(:,:,:) 
    F_h(:,:,:,my,3) = vz(:,:,:)*U_h(:,:,:,my) 

    ! z-momentum
    F_h(:,:,:,mz,1) = vx(:,:,:)*U_h(:,:,:,mz) 
    F_h(:,:,:,mz,2) = vy(:,:,:)*U_h(:,:,:,mz) 
    F_h(:,:,:,mz,3) = U_h(:,:,:,mz)**2/U_h(:,:,:,ro) &
         + pr(:,:,:) 

    ! energy
    F_h(:,:,:,en,1) = vx(:,:,:) &
         * ( U_h(:,:,:,en) + pr(:,:,:) ) 
    F_h(:,:,:,en,2) = vy(:,:,:) &
        * ( U_h(:,:,:,en) + pr(:,:,:) ) 
    F_h(:,:,:,en,3) = vz(:,:,:) &
         * ( U_h(:,:,:,en) + pr(:,:,:) ) 
   
    return
  end subroutine hd_fluxes
  subroutine mhd_fluxes(F_m,U_m)
    double precision,intent(inout):: F_m(ix,jx,kx,nvar_m,ndim)
    double precision,intent(inout):: U_m(ix,jx,kx,nvar_m)
    double precision :: de(ix,jx,kx),vx(ix,jx,kx)
    double precision ::vy(ix,jx,kx),vz(ix,jx,kx)
    double precision :: pr(ix,jx,kx)
    double precision ::bx(ix,jx,kx)
    double precision :: by(ix,jx,kx),bz(ix,jx,kx)
    double precision :: ex(ix,jx,kx),ey(ix,jx,kx),ez(ix,jx,kx)
    double precision :: ch
    integer, parameter :: ro=1
    integer, parameter :: mx=2
    integer, parameter :: my=3
    integer, parameter :: mz=4
    integer, parameter :: en=5
    integer, parameter :: b_x=6
    integer, parameter :: b_y=7
    integer, parameter :: b_z=8 
    integer, parameter :: ps=9
    integer :: n
    
    call cq2pv_mhd(de,vx,vy,vz,pr,bx,by,bz,U_m)
    !---------------------------------------------------------------------  
    ! ideal MHD part
    !---------------------------------------------------------------------  
    ! mass
    F_m(:,:,:,ro,1) = U_m(:,:,:,mx)
    ! momentum
    F_m(:,:,:,mx,1) = U_m(:,:,:,mx)**2/U_m(:,:,:,ro) &
         + pr(:,:,:) &
         + U_m(:,:,:,b_y)**2/(2.0d0) + U_m(:,:,:,b_z)**2/2.0d0 &
         - U_m(:,:,:,b_x)**2/(2.0d0)
    F_m(:,:,:,my,1) = vx(:,:,:)*U_m(:,:,:,my) &
         - U_m(:,:,:,b_x)*U_m(:,:,:,b_y)
    F_m(:,:,:,mz,1) = vx(:,:,:)*U_m(:,:,:,mz) &
         - U_m(:,:,:,b_x)*U_m(:,:,:,b_z)

    ey=-vz*bx+vx*bz
    ez=(-vx*by+vy*bx)
   
    F_m(:,:,:,b_x,1) =    0.0d0
    F_m(:,:,:,b_y,1)=-ez
    F_m(:,:,:,b_z,1)=ey    
    F_m(:,:,:,en,1) = vx(:,:,:) &
         * ( U_m(:,:,:,en) + pr(:,:,:) + 0.5d0*(U_m(:,:,:,b_x)**2 & 
         + U_m(:,:,:,b_y)**2+U_m(:,:,:,b_z)**2) ) &
         - U_m(:,:,:,b_x)*(vx(:,:,:)*U_m(:,:,:,b_x) &
               +vy(:,:,:)*U_m(:,:,:,b_y) &
               +vz(:,:,:)*U_m(:,:,:,b_z))
         
    if(ndim.ge.2) then
       ex=-vy*bz+vz*by
       F_m(:,:,:,ro,2) = U_m(:,:,:,my)
       F_m(:,:,:,mx,2) = vy(:,:,:)*U_m(:,:,:,mx) &
            - U_m(:,:,:,b_y)*U_m(:,:,:,b_x)
       F_m(:,:,:,my,2) = U_m(:,:,:,my)**2/U_m(:,:,:,ro) &
            + pr(:,:,:) &
            + U_m(:,:,:,b_x)**2/(2.0d0)+ U_m(:,:,:,b_z)**2/(2.0d0) & 
         - U_m(:,:,:,b_y)**2/(2.0d0)
       F_m(:,:,:,mz,2) = vy(:,:,:)*U_m(:,:,:,mz) &
            - U_m(:,:,:,b_y)*U_m(:,:,:,b_z)
       F_m(:,:,:,b_x,2) = ez
       F_m(:,:,:,b_y,2) =    0.0d0
       F_m(:,:,:,b_z,2) =-ex
                          
       F_m(:,:,:,en,2) = vy(:,:,:) &
            * ( U_m(:,:,:,en) + pr(:,:,:) + U_m(:,:,:,b_x)**2/(2.0d0) & 
            + U_m(:,:,:,b_y)**2/(2.0d0)+U_m(:,:,:,b_z)**2/(2.0d0) ) &
            - U_m(:,:,:,b_y)*(vx(:,:,:)*U_m(:,:,:,b_x) &
            +vy(:,:,:)*U_m(:,:,:,b_y) &
            +vz(:,:,:)*U_m(:,:,:,b_z))
       if(ndim.ge.3) then
          F_m(:,:,:,ro,3) = U_m(:,:,:,mz)
          F_m(:,:,:,mx,3) = vz(:,:,:)*U_m(:,:,:,mx) &
               - U_m(:,:,:,b_z)*U_m(:,:,:,b_x)
          F_m(:,:,:,my,3) = vz(:,:,:)*U_m(:,:,:,my) &
               - U_m(:,:,:,b_z)*U_m(:,:,:,b_y)
          F_m(:,:,:,mz,3) = U_m(:,:,:,mz)**2/U_m(:,:,:,ro) &
               + pr(:,:,:) &
               + U_m(:,:,:,b_x)**2/(2.0d0) + U_m(:,:,:,b_y)**2/(2.0d0) &
               - U_m(:,:,:,b_z)**2/(2.0d0)
          F_m(:,:,:,b_x,3) = -ey
          F_m(:,:,:,b_y,3) =  ex 
          F_m(:,:,:,b_z,3) =    0.0d0          
          F_m(:,:,:,en,3) = vz(:,:,:) &
               * ( U_m(:,:,:,en) + pr(:,:,:) + U_m(:,:,:,b_x)**2/(2.0d0)+ &
               U_m(:,:,:,b_y)**2/(2.0d0)+U_m(:,:,:,b_z)**2/(2.0d0) ) &
               - U_m(:,:,:,b_z)*(vx(:,:,:)*U_m(:,:,:,b_x) &
               +vy(:,:,:)*U_m(:,:,:,b_y) &
               +vz(:,:,:)*U_m(:,:,:,b_z))
       endif
    endif

    !! divB cleaning (original Dedner's 9-wave method)
    if(flag_divb.eq.1) then
       ch = cmax
       do n=1,ndim
          F_m(:,:,:,5+n,n) =U_m(:,:,:,9)
          F_m(:,:,:,9,n) = ch**2*U_m(:,:,:,5+n)
       enddo
    endif

    return
  end subroutine mhd_fluxes

  ! subroutine divb_flux(F_m,U_m)
  !   double precision,intent(inout):: F_m(ix,jx,kx,nvar_m,3)
  !   double precision,intent(in):: U_m(ix,jx,kx,nvar_m)
  !   double precision :: ch
  !   integer n
  !   ch = cmax
  !   if(flag_sch.eq.1) then
  !      do n=1,ndim
  !         F_m(:,:,:,5+n,n) =U_m(:,:,:,9)
  !         F_m(:,:,:,9,n) = ch**2*U_m(:,:,:,5+n)
  !      enddo
  !   else
  !      F_m(1:ix-1,:,:,6,1)=0.5d0*(U_m(1:ix-1,:,:,9)+U_m(2:ix,:,:,9)) 
  !      F_m(1:ix-1,:,:,9,1)=ch**2*0.5d0*(U_m(1:ix-1,:,:,6)+U_m(2:ix,:,:,6)) 
  !      if(ndim.ge.2) then
  !         F_m(:,1:jx-1,:,7,2)=0.5d0*(U_m(:,1:jx-1,:,9)+U_m(:,2:jx,:,9))
  !         F_m(:,1:jx-1,:,9,2)=ch**2*0.5d0*(U_m(:,1:jx-1,:,7)+U_m(:,2:jx,:,7))
  !         if(ndim.ge.3) then
  !            F_m(:,:,1:kx-1,8,3)=0.5d0*(U_m(:,:,1:kx-1,9)+U_m(:,:,2:kx,9))
  !            F_m(:,:,1:kx-1,9,3)=ch**2*0.5d0*(U_m(:,:,1:kx-1,8)+U_m(:,:,2:kx,8))
  !         endif
  !      endif
  !   endif
  ! end subroutine divb_flux


  ! subroutine divb_flux_tmp(F_m,U_m)
  !   double precision,intent(inout):: F_m(ix,jx,kx,nvar_m,3)
  !   double precision,intent(in):: U_m(ix,jx,kx,nvar_m)
  !   double precision :: ch
  !   integer n
  !   ch = cmax
  !   print*,'ch: ',ch
  !   if(flag_sch.eq.1) then
  !      do n=1,ndim
  !         F_m(:,:,:,5+n,n) =U_m(:,:,:,9)
  !         F_m(:,:,:,9,n) = ch**2*U_m(:,:,:,5+n)
  !      enddo
  !   else
  !      F_m(1:ix-1,:,:,6,1)=0.5d0*(U_m(1:ix-1,:,:,9)+U_m(2:ix,:,:,9)) 
  !      F_m(1:ix-1,:,:,9,1)=ch**2*0.5d0*(U_m(1:ix-1,:,:,6)+U_m(2:ix,:,:,6)) 
  !      if(ndim.ge.2) then
  !         F_m(:,1:jx-1,:,7,2)=0.5d0*(U_m(:,1:jx-1,:,9)+U_m(:,2:jx,:,9))
  !         F_m(:,1:jx-1,:,9,2)=ch**2*0.5d0*(U_m(:,1:jx-1,:,7)+U_m(:,2:jx,:,7))
  !         if(ndim.ge.3) then
  !            F_m(:,:,1:kx-1,8,3)=0.5d0*(U_m(:,:,1:kx-1,9)+U_m(:,:,2:kx,9))
  !            F_m(:,:,1:kx-1,9,3)=ch**2*0.5d0*(U_m(:,:,1:kx-1,8)+U_m(:,:,2:kx,8))
  !         endif
  !      endif
  !   endif
  ! end subroutine divb_flux_tmp

  
  subroutine divb_cleaning_Dedner_source(dt_sub,U_m)
    double precision,intent(in):: dt_sub
    double precision,intent(inout)::U_m(ix,jx,kx,nvar_m)
    integer, parameter :: ps=9
! divB cleaning
! source term integration for psi
     U_m(:,:,:,ps) = U_m(:,:,:,ps)*exp(-dt_sub/dt*safety*db_clean)    
    return
  end subroutine divb_cleaning_Dedner_source




  subroutine divb_cleaning_Dedner_iter(dt0,U_m)
    !Written by  S. Takasao
    ! b_cr and db_clean shoud be read from setting.txt file.
    ! b_cr: b_cr is a field strength value which satisfies b_cr > |divB*dx|
    ! Max iteration number is typically several, and less than 20.
    ! db_clean: db_clean is a ratio of the diffusion time to propagation time,
    !         and 0 < db_clean < 1. Typically db_clean = 0.5-0.7.
    double precision,intent(in):: dt0
    double precision,intent(inout)::U_m(ix,jx,kx,nvar_m)
    integer :: ncount
    double precision :: dpsdx(ix,jx,kx),dpsdy(ix,jx,kx),dpsdz(ix,jx,kx)
    double precision :: divb(ix,jx,kx),psi(ix,jx,kx),b_vec(ix,jx,kx,3)

    double precision :: eps,epsg,eps_limit,ch,cmaxg
    double precision,parameter::slow=0.9d0
    integer,parameter :: it_max = 100
    integer,parameter::i_bx=6,i_by=7,i_bz=8,i_psi=9

    ncount=0
    eps      = 1.d4
    eps_limit= 1.d-2 * b_cr ! 1% of the reliable field strength

    ! the propagation speed is slightly slowed down to avoid numerical instability.
    ! (not to let divB go to infinity)
    ch = slow*cmax

    b_vec = U_m(:,:,:,i_bx:i_bz)
    call div_cal(divb,b_vec)
    call bnd_divb(divb)
    psi = U_m(:,:,:,i_psi)

    do while(eps.ge.eps_limit .and. ncount.lt.it_max)
       !! initialization
       dpsdx = 0.d0
       dpsdy = 0.d0
       dpsdz = 0.d0

       psi = -(ch**2*dt0) * divb + (1.d0-db_clean)*psi       

       call derivative(dpsdx,psi,1)
       call derivative(dpsdy,psi,2)
       call derivative(dpsdz,psi,3)
       b_vec(:,:,:,1) = b_vec(:,:,:,1) - dt0*dpsdx
       b_vec(:,:,:,2) = b_vec(:,:,:,2) - dt0*dpsdy
       b_vec(:,:,:,3) = b_vec(:,:,:,3) - dt0*dpsdz
       call div_cal(divb,b_vec)
       call bnd_divb(divb)

       select case(ndim)
       case(2)
!          eps = sum(dsc(:,1:2))/size(dsc(:,1:2)) * maxval(abs(divb(:,:,1)))
          eps = min(minval(dx),minval(dy)) * maxval(abs(divb(:,:,1)))
       case(3)
!          eps = sum(dsc(:,1:3))/size(dsc(:,1:3)) * maxval(abs(divb))
          eps = min(minval(dx),minval(dy),minval(dz)) * maxval(abs(divb(:,:,:)))
       end select

       if(flag_mpi.eq.1) then
          eps = mpi_double_interface(eps,3)
       endif
!    if (MY_RANK.eq.0) print*,eps,eps_limit,ncount !test to check the number of iteration and convergence
       ncount = ncount+1
       if(ncount.eq.it_max) exit
    end do

    ! update
    U_m(:,:,:,i_bx:i_bz) = b_vec
    U_m(:,:,:,i_psi) = psi

    !    if(my_rank.eq.0) print*,'iteration #: ',ncount,'h*divB: ',eps

    return
  end subroutine divb_cleaning_Dedner_iter

  
  subroutine div_cal(div,vec)
    double precision,intent(in)::vec(ix,jx,kx,3)
    double precision,intent(out)::div(ix,jx,kx)
    if(s_order.eq.4) then
       if(ndim.eq.1) then
          do i=3,ix-2 
             div(i,:,:)=(-vec(i+2,:,:,1)+8.0d0*vec(i+1,:,:,1) &
                  -8.0d0*vec(i-1,:,:,1)+vec(i-2,:,:,1))/(12.0d0*dx(i))
          enddo
       else if(ndim.eq.2) then
          do j=3,jx-2
             do i=3,ix-2 
                div(i,j,:)=((-vec(i+2,j,:,1)+8.0d0*vec(i+1,j,:,1) &
                     -8.0d0*vec(i-1,j,:,1)+vec(i-2,j,:,1))/dx(i) &
                     +(-vec(i,j+2,:,2)+8.0d0*vec(i,j+1,:,2) &
                     -8.0d0*vec(i,j-1,:,2)+vec(i,j-2,:,2))/dy(j))/12.0d0
             enddo
          enddo
       else if(ndim.eq.3) then
          do k=3,kx-2
             do j=3,jx-2
                do i=3,ix-2 
                   div(i,j,k)=((-vec(i+2,j,k,1)+8.0d0*vec(i+1,j,k,1) &
                        -8.0d0*vec(i-1,j,k,1)+vec(i-2,j,k,1))/dx(i) &
                        +(-vec(i,j+2,k,2)+8.0d0*vec(i,j+1,k,2) &
                        -8.0d0*vec(i,j-1,k,2)+vec(i,j-2,k,2))/dy(j) &
                        +(-vec(i,j,k+2,3)+8.0d0*vec(i,j,k+1,3) &
                        -8.0d0*vec(i,j,k-1,3)+vec(i,j,k-2,3))/dz(k))/12.0d0
                enddo
             enddo
          enddo
       endif
    else 
       if(ndim.eq.1) then
          do i=2,ix-1 
             div(i,:,:)=(vec(i+1,:,:,1)-vec(i-1,:,:,1))/(2.0d0*dx(i))
          enddo
       else if(ndim.eq.2) then
          do j=2,jx-1
             do i=2,ix-1
                div(i,j,:)=((vec(i+1,j,:,1)-vec(i-1,j,:,1))/dx(i) &
                     +(vec(i,j+1,:,2)-vec(i,j-1,:,2))/dy(j))/2.0d0
             enddo
          enddo
       else if(ndim.eq.3) then
          do k=2,kx-1
             do j=2,jx-1
                do i=2,ix-1 
                   div(i,j,k)=((vec(i+1,j,k,1)-vec(i-1,j,k,1))/dx(i) &
                        +(vec(i,j+1,k,2)-vec(i,j-1,k,2))/dy(j) &
                        +(vec(i,j,k+1,3)-vec(i,j,k-1,3))/dz(k))/2.0d0
                enddo
             enddo
          enddo
       endif       
    endif
  end subroutine div_cal

  subroutine derivative(der,var,direct)
    integer,intent(in)::direct
    double precision,intent(in)::var(ix,jx,kx)
    double precision,intent(inout)::der(ix,jx,kx)
    
    !4th order derivative
    if (s_order .eq. 4) then
	    if(direct.eq.1) then 
	       do i=3,ix-2 
		  der(i,:,:)=(-var(i+2,:,:)+8.0d0*var(i+1,:,:) &
		       -8.0d0*var(i-1,:,:)+var(i-2,:,:))/(12.0d0*dx(i))
	       enddo
	    else if(direct.eq.2) then
	       do j=3,jx-2 
		  der(:,j,:)=(-var(:,j+2,:)+8.0d0*var(:,j+1,:) &
		       -8.0d0*var(:,j-1,:)+var(:,j-2,:))/(12.0d0*dy(j))
	       enddo
	    else if(direct.eq.3) then
	       do k=3,kx-2 
		  der(:,:,k)=(-var(:,:,k+2)+8.0d0*var(:,:,k+1) &
		       -8.0d0*var(:,:,k-1)+var(:,:,k-2))/(12.0d0*dz(k))
	       enddo
	    endif
    endif
    !1st order derivative
    if (s_order .eq. 1) then
	    if(direct.eq.1) then 
	       do i=2,ix-1 
		  der(i,:,:)=(var(i+1,:,:) &
		       -var(i-1,:,:))/(2.d0*dx(i))
	       enddo
	    else if(direct.eq.2) then
	       do j=2,jx-1 
		  der(:,j,:)=(var(:,j+1,:) &
		       -var(:,j-1,:))/(2.d0*dy(j))
	       enddo
	    else if(direct.eq.3) then
	       do k=2,kx-1 
		  der(:,:,k)=(var(:,:,k+1) &
		       -var(:,:,k-1))/(2.d0*dz(k))
	       enddo
	    endif
    endif
    
  end subroutine derivative

!  subroutine source_divb(S_m,U_m)
!    double precision,intent(inout)::S_m(ix,jx,kx,nvar_m)
!    double precision,intent(inout)::U_m(ix,jx,kx,nvar_m)
!  end subroutine source_divb

 subroutine artvis(dir,U,dt,nvar,mhd,n_target)    
    ! n: loop number
    integer, intent(in) :: dir,nvar,mhd,n_target
    double precision,intent(inout)::U(ix,jx,kx,nvar)
    double precision,intent(in)::dt
    integer :: i,j,k,var0
    integer :: is,js,ks
    double precision :: F(ix,jx,kx,n_target,ndim)
    double precision :: bb(ix,jx,kx),rovv(ix,jx,kx),pr(ix,jx,kx)
    double precision :: cc(ix,jx,kx),phi(ix,jx,kx,n_target)
    double precision :: a(ix,jx,kx,n_target) &
                       ,b(ix,jx,kx,n_target) &
                       ,c(ix,jx,kx,n_target)
    double precision :: sss(ix,jx,kx,n_target),de(ix,jx,kx)
    double precision :: UL(ix,jx,kx,n_target),UR(ix,jx,kx,n_target)
    double precision :: dU(ix,jx,kx,n_target)
    double precision :: dv(ix,jx,kx)
    double precision :: vx(ix,jx,kx),vy(ix,jx,kx),vz(ix,jx,kx)
    double precision :: bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
    double precision :: ww = 1.d-8
    integer xs,xe,ys,ye,zs,ze

    if (dir.eq.1) then
       is = 1
       js = 0
       ks = 0
       xs=2
       xe=ix-1
       ys=min(2,jx)
       ye=max(jx-1,1)
       zs=min(2,kx)
       ze=max(kx-1,1)
    elseif (dir.eq.2) then
       is = 0
       js = 1
       ks = 0
       xs=2
       xe=ix-1
       ys=2
       ye=jx-1
       zs=min(2,kx)
       ze=max(kx-1,1)
!    elseif (dir.eq.3) then
    else 
       is = 0
       js = 0
       ks = 1
       xs=2
       xe=ix-1
       ys=2
       ye=jx-1
       zs=2
       ze=kx-1
    endif
    
    !coefficient (velocity)
    if(mhd.eq.1)then
       bb(:,:,:) = U(:,:,:,6)**2 + U(:,:,:,7)**2 + U(:,:,:,8)**2
       call cq2pv_mhd(de,vx,vy,vz,pr,bx,by,bz,U)
       rovv=max(de,ro_lim)*(vx*vx+vy*vy+vz*vz)
       cc(:,:,:) = dsqrt( gm*max(pr(:,:,:),pr_lim)/max(U(:,:,:,1),ro_lim) ) &
            + dsqrt( rovv(:,:,:)/max(U(:,:,:,1),ro_lim) ) &
            + dsqrt( bb(:,:,:)/max(U(:,:,:,1),ro_lim) )
    else 
       call cq2pv_hd(de,vx,vy,vz,pr,U)
       cc=dsqrt(gm*max(pr(:,:,:),pr_lim)/max(de,ro_lim))+dsqrt(vx*vx+vy*vy+vz*vz)       
    endif

    cc(xs:xe,ys:ye,zs:ze) = &
         0.5d0*( cc(xs:xe,ys:ye,zs:ze) &
         +cc(xs+is:xe+is,ys+js:ye+js,zs+ks:ze+ks) )

    ! obtain the UR and UL 
    ! by using the generalized minmod function    
    a(xs:xe,ys:ye,zs:ze,:)  = &
         theta*(U(xs:xe,ys:ye,zs:ze,1:n_target) &
         -U(xs-is:xe-is,ys-js:ye-js,zs-ks:ze-ks,1:n_target))

    b(xs:xe,ys:ye,zs:ze,:)  = &
         0.5*(U(xs+is:xe+is,ys+js:ye+js,zs+ks:ze+ks,1:n_target) &
         -U(xs-is:xe-is,ys-js:ye-js,zs-ks:ze-ks,1:n_target))

    c(xs:xe,ys:ye,zs:ze,:) = &
         theta*(U(xs+is:xe+is,ys+js:ye+js,zs+ks:ze+ks,1:n_target) &
         -U(xs:xe,ys:ye,zs:ze,1:n_target))

    sss(:,:,:,:) = dble(floor( &
          sign(1.d0,a(:,:,:,:)) &
         +sign(1.d0,b(:,:,:,:)) &
         +sign(1.d0,c(:,:,:,:)) )/3)
    dU(:,:,:,:) = abs(sss(:,:,:,:)) &
         * ( 0.5*( sss(:,:,:,:)+1.d0) &
            *min(a(:,:,:,:),b(:,:,:,:),c(:,:,:,:)) &
            +0.5*(-sss(:,:,:,:)+1.d0) &
            *max(a(:,:,:,:),b(:,:,:,:),c(:,:,:,:)) )

    UL(xs:xe,ys:ye,zs:ze,:)=U(xs:xe,ys:ye,zs:ze,1:n_target) &
         +0.5d0*du(xs:xe,ys:ye,zs:ze,:)

    UR(xs:xe,ys:ye,zs:ze,:)=U(xs+is:xe+is,ys+js:ye+js,zs+ks:ze+ks,1:n_target) &
         -0.5d0*du(xs+is:xe+is,ys+js:ye+js,zs+ks:ze+ks,:)
    phi(:,:,:,:)=0.d0
    
    where((UR(xs:xe,ys:ye,zs:ze,:)-UL(xs:xe,ys:ye,zs:ze,:)) &
         *(U(xs+is:xe+is,ys+js:ye+js,zs+ks:ze+ks,1:n_target)- &
         U(xs:xe,ys:ye,zs:ze,1:n_target))>0.0d0) &
         phi(xs:xe,ys:ye,zs:ze,:)= &
         ((UL(xs:xe,ys:ye,zs:ze,:)-UR(xs:xe,ys:ye,zs:ze,:)) /&
         (U(xs+is:xe+is,ys+js:ye+js,zs+ks:ze+ks,1:n_target)- &
         U(xs:xe,ys:ye,zs:ze,1:n_target)))**2
!    phi(xs:xe,ys:ye,zs:ze,:)=-0.5d0*(sign(1.0d0, &
!         -(UR(xs:xe,ys:ye,zs:ze,:)-UL(xs:xe,ys:ye,zs:ze,:)) &
!         *(U(xs+is:xe+is,ys+js:ye+js,zs+ks:ze+ks,:)- &
!         U(xs:xe,ys:ye,zs:ze,:)))-1.0d0)*&
!         ((UL(xs:xe,ys:ye,zs:ze,:)-UR(xs:xe,ys:ye,zs:ze,:)) /&
!         (U(xs+is:xe+is,ys+js:ye+js,zs+ks:ze+ks,:)- &
!         U(xs:xe,ys:ye,zs:ze,:)))**2


    do var0=1,n_target
       F(:,:,:,var0,dir) = 0.5d0*cc(:,:,:)*phi(:,:,:,var0) &
            *( UL(:,:,:,var0)-UR(:,:,:,var0) )
    enddo
    if (dir.eq.1) then       
       dv(xs:xe,ys:ye,zs:ze)= &
            (abs(vx(xs+1:xe+1,ys:ye,zs:ze)-vx(xs:xe,ys:ye,zs:ze)))/ww
       do k=zs,ze
          do j=ys,ye
             do i=xs,xe       
                U(i,j,k,1:n_target) = U(i,j,k,1:n_target) &
                     - dt*( F(i,j,k,:,dir) &
                     -F(i-1,j,k,:,dir) )/dx(i) &
                     *tanh(dv(i,j,k))
             enddo
          enddo
       enddo
    elseif(dir.eq.2) then
       dv(xs:xe,ys:ye,zs:ze)= &
            (abs(vy(xs:xe,ys+1:ye+1,zs:ze)-vy(xs:xe,ys:ye,zs:ze)))/ww  
       do k=zs,ze
          do j=ys,ye
             do i=xs,xe       
                U(i,j,k,1:n_target) = U(i,j,k,1:n_target) &
                     - dt*( F(i,j,k,:,dir) &
                     -F(i,j-1,k,:,dir) )/dy(j) &
                     *tanh(dv(i,j,k))
             enddo
          enddo
       enddo       
    elseif (dir.eq.3) then
       dv(xs:xe,ys:ye,zs:ze)= &
            (abs(vz(xs:xe,ys:ye,zs+1:ze+1)-vz(xs:xe,ys:ye,zs:ze)))/ww
       do k=zs,ze
          do j=ys,ye
             do i=xs,xe       
                U(i,j,k,1:n_target) = U(i,j,k,1:n_target) &
                     - dt*( F(i,j,k,:,dir) &
                     -F(i,j,k-1,:,dir) )/dz(k) &
                     *tanh(dv(i,j,k))
             enddo
          enddo
       enddo
    end if
    return
  end subroutine artvis
  
  subroutine set_heat_conductivity(U_h,U_m)
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)
  end subroutine set_heat_conductivity

  subroutine set_viscosity(U_h,U_m)
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)
  end subroutine set_viscosity

  subroutine vel_damp(U_h,U_m)
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)
    double precision :: damp_time1(ix,jx,kx)
    double precision :: mach,rcom,rpres,f_n

!Damping based on sucsessive kinetic energy 
if ((flag_damp.eq.1).or.(flag_damp.eq.2)) then

    if (flag_damp.eq.2) then
	damp_time=1.0d0
	damp_time1=spread(spread(spread(1.d0*damp_time ,1,ix),2,jx),3,kx)
	do while (maxval(U_h(:,:,:,5)-damp_time1*dt*(U_h(:,:,:,2)**2+U_h(:,:,:,2)**2+U_h(:,:,:,2)**2)/U_h(:,:,:,1)/2.0d0) .GT. oldke_damp)
		damp_time=damp_time*2.0d0
		damp_time1=spread(spread(spread(1.d0*damp_time ,1,ix),2,jx),3,kx)
	enddo
!	print*,damp_time,maxval(U_h(:,:,:,5)-damp_time1*dt*(U_h(:,:,:,2)**2+U_h(:,:,:,2)**2+U_h(:,:,:,2)**2)/U_h(:,:,:,1)/2.0d0),oldke_damp
	damp_time=min(damp_time,1.0e6)
    endif

    damp_time1=spread(spread(spread(1.d0*damp_time ,1,ix),2,jx),3,kx)

    if (flag_mhd.eq.1) then
    U_m(:,:,:,2:4)=U_m(:,:,:,2:4)-spread(damp_time1,4,3)*dt*U_m(:,:,:,2:4)
    U_m(:,:,:,5)=U_m(:,:,:,5)-damp_time1*dt*(U_m(:,:,:,2)**2+U_m(:,:,:,3)**2 &
      +U_m(:,:,:,4)**2)/U_m(:,:,:,1)/2.d0
    endif
    if(flag_mhd.eq.0.or.flag_pip.eq.1) then
    U_h(:,:,:,2:4)=U_h(:,:,:,2:4)-spread(damp_time1,4,3)*dt*U_h(:,:,:,2:4)
    U_h(:,:,:,5)=U_h(:,:,:,5)-damp_time1*dt*(U_h(:,:,:,2)**2+U_h(:,:,:,3)**2 &
      +U_h(:,:,:,4)**2)/U_h(:,:,:,1)/2.d0
    endif

    oldke_damp=0.9d0*maxval(U_h(:,:,:,5)-(U_h(:,:,:,2)**2+U_h(:,:,:,2)**2+U_h(:,:,:,2)**2)/U_h(:,:,:,1)/2.0d0)


endif

if (flag_damp.eq.3) then
!Damping based on the shock frame parameters
	mach=2.d0
	rcom=(gm+1.d0)*mach**2/(2.d0+(gm-1.d0)*mach**2)
	rpres=1.d0+gm*mach**2*(1.d0-1.d0/rcom)

	damp_time1=spread(spread((tanh((x+70.d0)/10.0d0)+1.0d0)/2.0d0&
+(1.0d0-tanh((x-70.0d0)/10.0d0))/2.0d0-1.0d0,2,jx),3,kx)

!print*,damp_time1(:,1,1)
!print*,damp_time1(:,2,1)
!stop
f_n=0.0d0

if (flag_pip .eq. 1) then
	f_n=n_fraction
	U_h(:,:,:,5)=U_h(:,:,:,5)-0.5d0*(U_h(:,:,:,2)**2+U_h(:,:,:,3)**2+U_h(:,:,:,4)**2)/U_h(:,:,:,1)
	U_h(:,:,:,2)=(U_h(:,:,:,2)+f_n*mach)*damp_time1-spread(spread(spread(f_n*mach,1,ix),2,jx),3,kx)
	U_h(:,:,:,3)=(U_h(:,:,:,3))*damp_time1
	U_h(:,:,:,4)=(U_h(:,:,:,4))*damp_time1
	U_h(:,:,:,5)=U_h(:,:,:,5)+0.5d0*(U_h(:,:,:,2)**2+U_h(:,:,:,3)**2+U_h(:,:,:,4)**2)/U_h(:,:,:,1)
endif

	U_m(:,:,:,5)=U_m(:,:,:,5)-0.5d0*(U_m(:,:,:,2)**2+U_m(:,:,:,3)**2+U_m(:,:,:,4)**2)/U_m(:,:,:,1)
	U_m(:,:,:,2)=(U_m(:,:,:,2)+(1.0d0-f_n)*mach)*damp_time1-spread(spread(spread((1.0d0-f_n)*mach,1,ix),2,jx),3,kx)
	U_m(:,:,:,3)=(U_m(:,:,:,3))*damp_time1
	U_m(:,:,:,4)=(U_m(:,:,:,4))*damp_time1
	U_m(:,:,:,5)=U_m(:,:,:,5)+0.5d0*(U_m(:,:,:,2)**2+U_m(:,:,:,3)**2+U_m(:,:,:,4)**2)/U_m(:,:,:,1)	

!		print*,U_m(i,1,1,:)
endif

  end subroutine vel_damp

  subroutine get_vel_diff(vd,U_h,U_m)
  double precision,intent(inout)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)
  double precision,intent(out)::vd(ix,jx,kx,3)  
    integer :: i

     if (flag_pip.eq.1) then 
     do i=1,3
       vd(:,:,:,i)=(U_h(:,:,:,i+1)/U_h(:,:,:,1)-U_m(:,:,:,i+1)/U_m(:,:,:,1))
     enddo
     endif

  end subroutine get_vel_diff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine hydrogen_excitation_timestep(U_m,U_h,dnexcite)
! calculate the hydrogen excite timestep
  double precision,intent(in)::U_m(ix,jx,kx,nvar_m),U_h(ix,jx,kx,nvar_h)
  double precision,intent(out)::dnexcite(ix,jx,kx,5)
  double precision::dneut(6),rom(ix,jx,kx),roh(ix,jx,kx)
  double precision::dntot
  double precision::dneutv(ix,jx,kx,6)
  double precision::conv(ix,jx,kx,6),conv_temp(ix,jx,kx)
  double precision::dntotv(ix,jx,kx)
  integer::i,j,k,ii,jj,nmaxloc
  rom(:,:,:)=U_m(:,:,:,1)
  roh(:,:,:)=U_h(:,:,:,1)
!Vector form
!The new number of electrons/protons is:
            ! Calculate the change in each neutral species
            dntotv(:,:,:)=0.d0
            dneutv(:,:,:,:)=0.d0
            do ii=1,6
                do jj=1,6 
                    dneutv(:,:,:,ii)=dneutv(:,:,:,ii)+Nexcite(:,:,:,jj)*colrat(:,:,:,jj,ii)/Gm_rec_ref*t_ir - &
                    			Nexcite(:,:,:,ii)*colrat(:,:,:,ii,jj)/Gm_rec_ref*t_ir
                    if (flag_rad .ge. 2) then
                        dneutv(:,:,:,ii)=dneutv(:,:,:,ii)+Nexcite(:,:,:,jj)*radrat(:,:,:,jj,ii)/Gm_rec_ref*t_ir - &
                        			Nexcite(:,:,:,ii)*radrat(:,:,:,ii,jj)/Gm_rec_ref*t_ir
                    endif
                enddo
                
                !Convective term (neutrals)
!                if (s_order .ne. 4) then
!                	print*,'Convective term only works in 4th order at the moment'
!                	stop
!            	endif
                if (ii .le. 5) then
		            call derivative(conv_temp,Nexcite(:,:,:,ii),1)
		            conv(:,:,:,ii)=conv_temp*U_h(:,:,:,2)/U_h(:,:,:,1)
		            if (ndim .gt. 1) then
		            	call derivative(conv_temp,Nexcite(:,:,:,ii),2)
		            	conv(:,:,:,ii)=conv(:,:,:,ii)+conv_temp*U_h(:,:,:,3)/U_h(:,:,:,1)
		           		if (ndim .gt. 2) then
		           			call derivative(conv_temp,Nexcite(:,:,:,ii),3)
		           			conv(:,:,:,ii)=conv(:,:,:,ii)+conv_temp*U_h(:,:,:,4)/U_h(:,:,:,1)
	           			endif
           			endif
       			endif
            enddo

	    dnexcite(:,:,:,1:5)=dneutv(:,:,:,1:5)-conv(:,:,:,1:5)
  end subroutine hydrogen_excitation_timestep
  
end module scheme_rot
