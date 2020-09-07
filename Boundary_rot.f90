module boundary_rot
!====================================================================
! This is the boundary module for the PIP code.
! Module for calling the different boundaries
! Author A.Hillier
! First version - 2013/05/31 (Fixed routines not completed)
! Second version - 2013/06/27 (New labelling for boundaries)
!====================================================================
  use globalvar,only:ix,jx,kx,ndim,flag_pip,flag_mhd,nvar_h,nvar_m,&
       flag_mpi,flag_bnd,margin,neighbor,flag_divb,time,x,gm,dt,nt,t_order, beta, n_fraction
  use mpi_rot,only:mpi_bnd,mpi_bnd_onevar
  use parameters,only:pi !for the periodic drivers

  implicit none
  integer,allocatable,save ::sym_mhd(:,:)  
  integer,allocatable,save ::sym_hd(:,:)  
contains 
!------------------------------------------------------------------------
  subroutine initialize_bnd
    if (flag_pip.eq.1) then       
       allocate(sym_mhd(nvar_m,6),sym_hd(nvar_h,6))
       call set_sym_mhd      
       call set_sym_hd 
    elseif (flag_mhd.eq.1) then       
       allocate(sym_mhd(nvar_m,6))
       call set_sym_mhd   
    else       
       allocate(sym_hd(nvar_h,6))
       call set_sym_hd 
    endif    
  end subroutine initialize_bnd

  subroutine PIPbnd(U_h,U_m,istep)
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)
    integer, intent(in)::istep

    if(flag_mpi.eq.1) then
       call mpi_bnd(U_h,U_m)
    endif    
    if (flag_pip.eq.1) then       
       call MHDbnd(U_m,istep)   
       call HDbnd(U_h,istep)   
    elseif (flag_mhd.eq.1) then       
       call MHDbnd(U_m,istep)          
    else       
       call HDbnd(U_h,istep)          
    endif    



end subroutine PIPbnd

  subroutine MHDbnd(U_m,istep)   
    double precision,intent(inout)::U_m(ix,jx,kx,nvar_m)
    integer,intent(in)::istep
    integer dir,upper_lower,type,i
    do i=1,2*ndim
       dir=(i+1)/2
       upper_lower=mod(i+1,2)
!print*,upper_lower       
       if(flag_bnd(i).eq.3) then
          type=2
       else if(flag_bnd(i).eq.4) then
          type=2
       else if(flag_bnd(i).eq.5) then
          type=2
       else if(flag_bnd(i).eq.6) then
          type=2
       else if(flag_bnd(i).eq.11) then
          type=2
!       else if(flag_bnd(i).eq.20) then
!          type=2
       else
          type=flag_bnd(i)
       endif

       if(neighbor(i).eq.-1) then
	if(flag_bnd(i) .EQ. 20) then
          call boundary_control_custom_mhd(U_m,ix,jx,kx,nvar_m,upper_lower,istep)
	else
          call boundary_control(U_m,ix,jx,kx,nvar_m,margin,dir,upper_lower, &
               sym_mhd(:,i),type)
	endif
       endif
    enddo
  end subroutine MHDbnd

  subroutine HDbnd(U_h,istep)   
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h)
    integer,intent(in)::istep
    integer dir,upper_lower,type,i
    do i=1,2*ndim
       dir=(i+1)/2
       upper_lower=mod(i+1,2)       
       if(flag_bnd(i).eq.3) then
          type=2
       else if(flag_bnd(i).eq.4) then
          type=2
       else if(flag_bnd(i).eq.5) then
          type=2
       else if(flag_bnd(i).eq.6) then
          type=2          
       else if(flag_bnd(i).eq.11) then
          type=2
       else
          type=flag_bnd(i)
       endif
       if(neighbor(i).eq.-1) then
	if(flag_bnd(i) .EQ. 20) then
          call boundary_control_custom_hd(U_h,ix,jx,kx,nvar_h,upper_lower,istep)
	else
          call boundary_control(U_h,&
               ix,jx,kx,nvar_h,margin,dir,upper_lower, &
               sym_hd(:,i),type)
	endif
       endif
    enddo
  end subroutine HDbnd

  subroutine bnd_energy(en1)   
    double precision,intent(inout)::en1(ix,jx,kx)
    integer dir,upper_lower,type,i
    integer :: sym0(6)
    do i=1,2*ndim
       dir=(i+1)/2
       upper_lower=mod(i+1,2)       
       if(flag_bnd(i).eq.3) then
          type=2
       else if(flag_bnd(i).eq.4) then
          type=2
       else if(flag_bnd(i).eq.5) then
          type=2
       else if(flag_bnd(i).eq.6) then
          type=2
       else if(flag_bnd(i).eq.11) then
          type=2
       else
          type=flag_bnd(i)
       endif

       select case(flag_mhd)
       case(0)
          sym0(:) = sym_hd(5,:)
       case(1)
          sym0(:) = sym_mhd(5,:)
       end select

       if(flag_mpi.eq.1) call mpi_bnd_onevar(en1)
       if(neighbor(i).eq.-1) then
          call boundary_control_onevar(en1,ix,jx,kx,margin,dir,upper_lower, &
               sym0(i),type)
       endif
    enddo
  end subroutine Bnd_Energy


  subroutine bnd_divb(divb1)   
    double precision,intent(inout)::divb1(ix,jx,kx)
    integer dir,upper_lower,type,i
    integer :: sym0(6)
    do i=1,2*ndim
       dir=(i+1)/2
       upper_lower=mod(i+1,2)       
       if(flag_bnd(i).eq.3) then
          type=2
       else if(flag_bnd(i).eq.4) then
          type=2
       else if(flag_bnd(i).eq.5) then
          type=2
       else if(flag_bnd(i).eq.6) then
          type=2          
       else if(flag_bnd(i).eq.11) then
          type=2
!       else if(flag_bnd(i).eq.20) then
!          type=2
       else
          type=flag_bnd(i)
       endif


       sym0(:) = sym_mhd(9,:) !! Dedner's psi

       if(flag_mpi.eq.1) call mpi_bnd_onevar(divb1)       
       if(neighbor(i).eq.-1) then
          call boundary_control_onevar(divb1,ix,jx,kx,margin,dir,upper_lower, &
               sym0(i),type)
       endif
    enddo
  end subroutine Bnd_Divb

  ! subroutine bnd_gravity(gra1)   
  !   double precision,intent(inout)::gra1(ix,jx,kx,3)
  !   integer dir,upper_lower,type,i
  !   integer :: sym0(6)
  !   do i=1,2*ndim
  !      dir=(i+1)/2
  !      upper_lower=mod(i+1,2)       
  !      if(flag_bnd(i).eq.3) then
  !         type=2
  !      else if(flag_bnd(i).eq.4) then
  !         type=2
  !      else if(flag_bnd(i).eq.5) then
  !         type=2
  !      else if(flag_bnd(i).eq.6) then
  !         type=2
  !      else if(flag_bnd(i).eq.11) then
  !         type=2
  !      else
  !         type=flag_bnd(i)
  !      endif

  !      select case(flag_mhd)
  !      case(0)
  !         sym0(:) = sym_hd(5,:)
  !      case(1)
  !         sym0(:) = sym_mhd(5,:)
  !      end select

  !      if(flag_mpi.eq.1) call mpi_bnd_onevar(en1)
  !      if(neighbor(i).eq.-1) then
  !         call boundary_control_onevar(en1,ix,jx,kx,margin,dir,upper_lower, &
  !              sym0(i),type)
  !      endif
  !   enddo
  ! end subroutine Bnd_Energy
  
  
  subroutine set_sym_mhd
    integer dir,i
    sym_mhd(:,:)=1
    do i=1,2*ndim
       dir=(i+1)/2
       if(flag_bnd(i).eq.2) then
          sym_mhd(1+dir,i)=-1            !! vn
          sym_mhd(5+dir,i)=-1            !! bn
       else if(flag_bnd(i).eq.3) then
          sym_mhd(1+dir,i)=-1            !! vn
          sym_mhd(6+mod(dir,3),i)=-1     !! bn
          sym_mhd(2+mod(dir+1,3),i)=-1   !! vt2  
       else if(flag_bnd(i).eq.4) then
          sym_mhd(1+dir,i)=-1            !! vn
          sym_mhd(6+mod(dir+1,3),i)=-1   !! bt2
          sym_mhd(2+mod(dir,3),i)=-1     !! vt1
!          sym_mhd(6+mod(dir+1,3),i)=-1          
       else if(flag_bnd(i).eq.5) then
          sym_mhd(1+dir,i)=-1            !! vn
          sym_mhd(5+dir,i)=-1            !! bn
          sym_mhd(2+mod(dir+1,3),i)=-1   !! vt2
          sym_mhd(6+mod(dir+1,3),i)=-1   !! bt2
       else if(flag_bnd(i).eq.6) then
          sym_mhd(1+dir,i)=-1            !! vn
          sym_mhd(6+mod(dir,3),i)=-1     !! bt1
          sym_mhd(6+mod(dir+1,3),i)=-1   !! bt2       
       endif
    enddo

    if(flag_divb.eq.1) then
       do i=1,2*ndim
          dir=(i+1)/2
          if(flag_bnd(i).eq.2) then
             sym_mhd(9,i) = 1  ! Bn is odd, so psi is even
          else if(flag_bnd(i).eq.20) then
             sym_mhd(9,i) = 1  ! Bn is odd, so psi is even
          else if(flag_bnd(i).eq.3) then
             sym_mhd(9,i) = -1 ! Bn is even, so psi is odd
          else if(flag_bnd(i).eq.4) then
             sym_mhd(9,i) = -1 ! Bn is even, so psi is odd
          else if(flag_bnd(i).eq.5) then
             sym_mhd(9,i) = 1 ! Bn is odd, so psi is even
          else if(flag_bnd(i).eq.6) then
             sym_mhd(9,i) = -1 ! Bn is even, so psi is odd
          endif
       enddo
    endif

  end subroutine set_sym_mhd

  subroutine set_sym_hd

    integer dir,i
    sym_hd(:,:)=1
    do i=1,2*ndim
       dir=(i+1)/2
       if(flag_bnd(i).eq.2) then
          sym_hd(1+dir,i)=-1
       else if(flag_bnd(i).eq.20) then
          sym_hd(1+dir,i)=-1
       else if(flag_bnd(i).eq.3) then
          sym_hd(1+dir,i)=-1
          sym_hd(2+mod(dir+1,3),i)=-1          
       else if(flag_bnd(i).eq.4) then
          sym_hd(1+dir,i)=-1
!          sym_hd(2+mod(dir+1,3),i)=-1       
          sym_hd(2+mod(dir,3),i)=-1
       else if(flag_bnd(i).eq.5) then
          sym_hd(1+dir,i)=-1            !! vn
          sym_hd(2+mod(dir+1,3),i)=-1   !! vt2
       else if(flag_bnd(i).eq.6) then          
          sym_hd(1+dir,i)=-1            !! vn
       endif
    enddo
  end subroutine set_sym_hd
  

  subroutine boundary_control(U,ix,jx,kx,nvar,margin,dir,upper_lower,sym,type)
    implicit none
    integer,intent(in)::ix,jx,kx,nvar
    integer,intent(in)::dir,upper_lower,sym(nvar),type,margin(3)
    double precision,intent(inout)::U(ix,jx,kx,nvar) 
    integer ig(3),sl(3),el(3),dl(3),sr(3),er(3),dr(3)
    integer i,n
    ig(1)=ix
    ig(2)=jx
    ig(3)=kx
    
    !caluculate region for each direction
    do i=1,3
       if(dir.eq.i) then 
          !lower
          if (upper_lower.eq.0) then
             sl(i)=1 ; el(i)=margin(i) ; dl(i)=1
             select case(type)
                !periodic
             case(1)
                sr(i)=ig(i)-2*margin(i)+1 ; er(i)=ig(i)-margin(i) ; dr(i)=1     
                !symmetric or free2
             case(2)
                sr(i)=2*margin(i) ; er(i)=margin(i)+1 ; dr(i)=-1             
                !free-1
             case(10)
                sr(i)=margin(i)+1 ; er(i)=margin(i)+1 ; dr(i)=1            
             end select
             !upper
          else if(upper_lower.eq.1) then
             sl(i)=ig(i)-margin(i)+1 ; el(i)=ig(i) ; dl(i)=1
             select case(type)
                !periodic
             case(1)
                sr(i)=margin(i)+1 ; er(i)=2*margin(i) ; dr(i)=1     
                !symmetric or !free-2
             case(2)
                sr(i)=ig(i)-margin(i) ; er(i)=ig(i)-margin(i)*2+1 ; dr(i)=-1  
                !free-1
             case(10)
                sr(i)=ig(i)-margin(i)   ; er(i)=ig(i)-margin(i) ; dr(i)=1  
             end select
          endif
          !all
       else
          sl(i)=1 ; el(i)=ig(i) ; dl(i)=1
          sr(i)=1 ; er(i)=ig(i) ; dr(i)=1
       endif
    enddo
    
    
    if ((el(1)-sl(1).eq.abs(er(1)-sr(1))).and. &
         (el(2)-sl(2).eq.abs(er(2)-sr(2))).and. &
         (el(3)-sl(3).eq.abs(er(3)-sr(3)))) then

       do n=1,nvar
          u(sl(1):el(1):dl(1),sl(2):el(2):dl(2),sl(3):el(3):dl(3),n) &
               =sym(n)*u(sr(1):er(1):dr(1),sr(2):er(2):dr(2),sr(3):er(3):dr(3),n)
       enddo
    else
       !     allocate(buf(el(1)-sl(1)+1,el(2)-sl(2)+1,el(3)-sl(3)+1,nvar))
       !not same size 
       !free -boundary 1
       if(type.eq.10) then 
          if(dir.eq.1) then 
             do n=1,margin(1)
                u(sl(1)+n-1,:,:,:)=u(sr(1),:,:,:)
             enddo
          else if(dir.eq.2) then   
             do n=1,margin(2)
                u(:,sl(2)+n-1,:,:)=u(:,sr(2),:,:)
             enddo
          else if(dir.eq.3) then   
             do n=1,margin(3)
                u(:,:,sl(3)+n-1,:)=u(:,:,sr(3),:)
             enddo
          endif
       endif
    endif
  end subroutine boundary_control


  subroutine boundary_control_onevar(u1,ix,jx,kx,margin,dir,upper_lower,sym1,type)
    implicit none
    integer,intent(in)::ix,jx,kx
    integer,intent(in)::dir,upper_lower,sym1,type,margin(3)
    double precision,intent(inout)::u1(ix,jx,kx) 
    integer ig(3),sl(3),el(3),dl(3),sr(3),er(3),dr(3)
    integer i,n
    ig(1)=ix
    ig(2)=jx
    ig(3)=kx
    
    !caluculate region for each direction
    do i=1,3
       if(dir.eq.i) then 
          !lower
          if (upper_lower.eq.0) then
             sl(i)=1 ; el(i)=margin(i) ; dl(i)=1
             select case(type)
                !periodic
             case(1)
                sr(i)=ig(i)-2*margin(i)+1 ; er(i)=ig(i)-margin(i) ; dr(i)=1     
                !symmetric or free2
             case(2)
                sr(i)=2*margin(i) ; er(i)=margin(i)+1 ; dr(i)=-1             
                !free-1
             case(10)
                sr(i)=margin(i)+1 ; er(i)=margin(i)+1 ; dr(i)=1            
             end select
             !upper
          else if(upper_lower.eq.1) then
             sl(i)=ig(i)-margin(i)+1 ; el(i)=ig(i) ; dl(i)=1
             select case(type)
                !periodic
             case(1)
                sr(i)=margin(i)+1 ; er(i)=2*margin(i) ; dr(i)=1     
                !symmetric or !free-2
             case(2)
                sr(i)=ig(i)-margin(i) ; er(i)=ig(i)-margin(i)*2+1 ; dr(i)=-1  
                !free-1
             case(10)
                sr(i)=ig(i)-margin(i)   ; er(i)=ig(i)-margin(i) ; dr(i)=1  
             end select
          endif
          !all
       else
          sl(i)=1 ; el(i)=ig(i) ; dl(i)=1
          sr(i)=1 ; er(i)=ig(i) ; dr(i)=1
       endif
    enddo
    
    
    if ((el(1)-sl(1).eq.abs(er(1)-sr(1))).and. &
         (el(2)-sl(2).eq.abs(er(2)-sr(2))).and. &
         (el(3)-sl(3).eq.abs(er(3)-sr(3)))) then

       u1(sl(1):el(1):dl(1),sl(2):el(2):dl(2),sl(3):el(3):dl(3)) &
            =sym1*u1(sr(1):er(1):dr(1),sr(2):er(2):dr(2),sr(3):er(3):dr(3))
    else
       !     allocate(buf(el(1)-sl(1)+1,el(2)-sl(2)+1,el(3)-sl(3)+1,nvar))
       !not same size 
       !free -boundary 1
       if(type.eq.10) then 
          if(dir.eq.1) then 
             do n=1,margin(1)
                u1(sl(1)+n-1,:,:)=u1(sr(1),:,:)
             enddo
          else if(dir.eq.2) then   
             do n=1,margin(2)
                u1(:,sl(2)+n-1,:)=u1(:,sr(2),:)
             enddo
          else if(dir.eq.3) then   
             do n=1,margin(3)
                u1(:,:,sl(3)+n-1)=u1(:,:,sr(3))
             enddo
          endif
       endif
    endif
  end subroutine boundary_control_onevar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine boundary_control_custom_hd(U,ix,jx,kx,nvar,dir,istep)
!Custom boundary conditions (bnd_flag=20)
!Velocity driver
implicit none
integer,intent(in)::dir,ix,jx,kx,nvar,istep
double precision,intent(inout)::U(ix,jx,kx,nvar) 
double precision :: vxpert,vdxpert,period,ppert,vamp
double precision :: rodxpert,ropert,cs,pdxpert,pr0,ro0
double precision :: f_n,f_p,f_p_p,f_p_n
double precision :: mach,rcom,rpres

period=0.5d0 
vamp=2.0e-1
ro0=276.17855d0 !xi_n=0.9
pr0=82.853565d0 !xi_n=0.9
!ro0=3037.9641d0 !xi_n=0.99
!pr0=911.38922d0 !xi_n=0.99

!if(dir.eq.1) then 
!	if ((time .le. 0.5d0*period) .and. (istep .eq. 0)) then
!		cs=sqrt(gm*pr0/ro0)
!		vxpert=vamp*dsin(time*2.0d0*pi/period)
!		vdxpert=vamp*2.0d0*pi/period*dcos(time*2.0d0*pi/period)
!		ppert=pr0*gm/cs*vxpert 
!		pdxpert=pr0*gm/cs*vdxpert
!		ropert=ro0/cs*vxpert
!		rodxpert=ro0/cs*vdxpert
!	else
!		vxpert=0.0d0
!		vdxpert=0.0d0
!		ppert=0.0d0
!		pdxpert=0.0d0
!		ropert=0.0d0
!		rodxpert=0.0d0
!	endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Drive velocity, density and pressure
!	u(1,:,:,1)=u(1,:,:,1)+rodxpert*dt
!	u(1,:,:,2)=u(1,:,:,2)+(rodxpert*vxpert+(ropert+ro0)*vdxpert)*dt
!	u(1,:,:,3)=u(1,:,:,3) !vy momentum 
!	u(1,:,:,4)=u(1,:,:,4) !vz momentum 
!	u(1,:,:,5)=u(1,:,:,5) +(0.5d0*vxpert**2.0d0*rodxpert+!(ropert+ro0)*vxpert*vdxpert+pdxpert/(gm-1.0d0))*dt
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!endif

  !set ionization fraction-----------------
  if(flag_pip.eq.0) then
     f_n=1.0d0
     f_p=1.0d0
     f_p_n=1.0d0
     f_p_p=1.0d0
  else
     f_n=n_fraction
     f_p=1.0d0-n_fraction     
     f_p_n=f_n/(f_n+2.0d0*f_p)
     f_p_p=2.0d0*f_p/(f_n+2.0d0*f_p)
  endif

if(dir.eq.0) then 
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Fix boundary values - parallel shock
	mach=2.d0
	rcom=(gm+1.d0)*mach**2/(2.d0+(gm-1.d0)*mach**2)
	rpres=1.d0+gm*mach**2*(1.d0-1.d0/rcom)

	u(1,:,:,1)=f_n*rcom
	u(1,:,:,2)=-mach*f_n
	u(1,:,:,3)=0.0d0!-0.1d0*u(1,:,:,1) 
	u(1,:,:,4)=0.0d0 !vz momentum 
	u(1,:,:,5)=(f_p_n*rpres/gm)/(gm-1.d0)+0.5d0*(f_n*mach**2)/rcom! +0.5d0*(2.d0/gm/beta)

	u(2,:,:,1)=f_n*rcom
	u(2,:,:,2)=-mach*f_n
	u(2,:,:,3)=0.0d0!-0.1d0*u(1,:,:,1) 
	u(2,:,:,4)=0.0d0 !vz momentum 
	u(2,:,:,5)=(f_p_n*rpres/gm)/(gm-1.d0)+0.5d0*(f_n*mach**2)/rcom !+0.5d0*(2.d0/gm/beta)

	u(3,:,:,1)=f_n*rcom
	u(3,:,:,2)=-f_n*mach
	u(3,:,:,3)=0.0d0!-0.1d0*u(1,:,:,1) 
	u(3,:,:,4)=0.0d0 !vz momentum 
	u(3,:,:,5)=(f_p_n*rpres/gm)/(gm-1.d0)+0.5d0*(f_n*mach**2)/rcom !+0.5d0*(2.d0/gm/beta)

	u(4,:,:,1)=f_n*rcom
	u(4,:,:,2)=-f_n*mach
	u(4,:,:,3)=0.0d0!-0.1d0*u(1,:,:,1) 
	u(4,:,:,4)=0.0d0 !vz momentum 
	u(4,:,:,5)=(f_p_n*rpres/gm)/(gm-1.d0)+0.5d0*(f_n*mach**2)/rcom !+0.5d0*(2.d0/gm/beta)

endif

if (dir .eq. 1) then 
	mach=2.d0
	rcom=(gm+1.d0)*mach**2/(2.d0+(gm-1.d0)*mach**2)
	rpres=1.d0+gm*mach**2*(1.d0-1.d0/rcom)

	u(ix,:,:,1)=f_n*1.0d0
	u(ix,:,:,2)=-f_n*mach
	u(ix,:,:,3)=0.0d0!-0.1d0*u(1,:,:,1) 
	u(ix,:,:,4)=0.0d0 !vz momentum 
	u(ix,:,:,5)=(f_p_n*1.d0/gm)/(gm-1.d0)+0.5d0*f_n*mach**2 !+0.5d0*(2.d0/gm/beta)

	u(ix-1,:,:,1)=f_n*1.0d0
	u(ix-1,:,:,2)=-f_n*mach
	u(ix-1,:,:,3)=0.0d0!-0.1d0*u(1,:,:,1) 
	u(ix-1,:,:,4)=0.0d0 !vz momentum 
	u(ix-1,:,:,5)=(f_p_n*1.d0/gm)/(gm-1.d0)+0.5d0*f_n*mach**2!+0.5d0*(2.d0/gm/beta)

	u(ix-2,:,:,1)=f_n*1.0d0
	u(ix-2,:,:,2)=-f_n*mach
	u(ix-2,:,:,3)=0.0d0!-0.1d0*u(1,:,:,1) 
	u(ix-2,:,:,4)=0.0d0 !vz momentum 
	u(ix-2,:,:,5)=(f_p_n*1.d0/gm)/(gm-1.d0)+0.5d0*f_n*mach**2!+0.5d0*(2.d0/gm/beta)

	u(ix-3,:,:,1)=f_n*1.0d0
	u(ix-3,:,:,2)=-f_n*mach
	u(ix-3,:,:,3)=0.0d0!-0.1d0*u(1,:,:,1) 
	u(ix-3,:,:,4)=0.0d0 !vz momentum 
	u(ix-3,:,:,5)=(f_p_n*1.d0/gm)/(gm-1.d0)+0.5d0*f_n*mach**2!+0.5d0*(2.d0/gm/beta)
!!!!!!!!!!!!!!!!!!!!!!!!!!!
endif
end subroutine boundary_control_custom_hd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine boundary_control_custom_mhd(U,ix,jx,kx,nvar,dir,istep)
!Custom boundary conditions (bnd_flag=20)
!Velocity driver
implicit none
integer,intent(in)::dir,ix,jx,kx,nvar,istep
double precision,intent(inout)::U(ix,jx,kx,nvar) 
double precision :: vxpert,vdxpert,period,ppert,vamp
double precision :: rodxpert,ropert,cs,pdxpert,pr0,ro0
double precision mach,rcom,rpres,alf,ang, byrat,vyrat, vu,bxu,byu,rou,pru,vxu
double precision :: f_n,f_p,f_p_p,f_p_n

period=0.5d0 
vamp=1.0e-2
!ro0=1.0d0
ro0=30.686506d0
!pr0=0.6d0
pr0=18.411901d0

  !set ionization fraction-----------------
  if(flag_pip.eq.0) then
     f_n=1.0d0
     f_p=1.0d0
     f_p_n=1.0d0
     f_p_p=1.0d0
  else
     f_n=n_fraction
     f_p=1.0d0-n_fraction     
     f_p_n=f_n/(f_n+2.0d0*f_p)
     f_p_p=2.0d0*f_p/(f_n+2.0d0*f_p)
  endif

if(dir.eq.0) then 
!	if ((time .le. 0.5d0*period) .and. (istep .eq. 0)) then
!		cs=sqrt(gm*pr0/ro0)
!		vxpert=vamp*dsin(time*2.0d0*pi/period)
!		vdxpert=vamp*2.0d0*pi/period*dcos(time*2.0d0*pi/period)
!		ppert=pr0*gm/cs*vxpert 
!		pdxpert=pr0*gm/cs*vdxpert
!		ropert=ro0/cs*vxpert
!		rodxpert=ro0/cs*vdxpert
!	else
!		vxpert=0.0d0
!		vdxpert=0.0d0
!		ppert=0.0d0
!		pdxpert=0.0d0
!		ropert=0.0d0
!		rodxpert=0.0d0
!	endif
!	rodxperp=30.686506d0*vdxpert*ppert/gm
!	roperp=30.686506d0*vxpert*ppert/gm
!	rodxperp=vdxpert*ppert/gm
!	roperp=vxpert*ppert/gm
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Drive velocity only.
!	u(1,:,:,1)=u(1,:,:,1)
!	u(1,:,:,2)=u(1,:,:,2)+u(1,:,:,1)*vxpert 
!	u(1,:,:,3)=u(1,:,:,3) !vy momentum 
!	u(1,:,:,4)=u(1,:,:,4) !vz momentum 
!	u(1,:,:,5)=u(1,:,:,5) +0.5d0*u(1,:,:,1)* vxpert**2.0d0
!	u(1,:,:,6)=u(1,:,:,6) !bx
!	u(1,:,:,7)=u(1,:,:,7) !by
!	u(1,:,:,8)=u(1,:,:,8) !bz
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Drive velocity, density and pressure
!	u(1,:,:,1)=u(1,:,:,1)+rodxpert*dt
!	u(1,:,:,2)=u(1,:,:,2)+(rodxpert*vxpert+(ropert+ro0)*vdxpert)*dt
!	u(1,:,:,3)=u(1,:,:,3)+(rodxpert*vxpert+(ropert+ro0)*vdxpert)*dt !vy momentum 
!	u(1,:,:,4)=u(1,:,:,4) !vz momentum 
!	u(1,:,:,5)=u(1,:,:,5) +(0.5d0*vxpert**2.0d0*rodxpert+(ropert+ro0)*vxpert*vdxpert+pdxpert/(gm-1.0d0))*dt
!	u(1,:,:,6)=u(1,:,:,6) !bx
!	u(1,:,:,7)=u(1,:,:,7) !by
!	u(1,:,:,8)=u(1,:,:,8) !bz
!	print*,u(1,1,1,2)/u(1,1,1,1),vxpert
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Fix boundary values.
!	u(1,:,:,1)=1.368d0!u(1,:,:,1)
!	u(1,:,:,2)=u(1,:,:,1)*0.269d0 
!	u(1,:,:,3)=u(1,:,:,1)*1.d0 
!	u(1,:,:,4)=u(1,:,:,1)*0.0d0 !vz momentum 
!	u(1,:,:,5)=1.769d0/(gm-1)+0.5d0*u(1,:,:,1)* (1.d0+0.269**2.0d0)
!	u(1,:,:,6)=1.d0 !bx
!	u(1,:,:,7)=0.d0!u(1,:,:,7) !by
!	u(1,:,:,8)=0.d0!u(1,:,:,8) !bz
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Fix boundary values - parallel shock
	mach=2.d0
	rcom=(gm+1.d0)*mach**2/(2.d0+(gm-1.d0)*mach**2)
	rpres=1.d0+gm*mach**2*(1.d0-1.d0/rcom)
!     v_l=(/rcom,rpres/gm,-mach/rcom,0.0d0,0.0d0,dsqrt(2.d0/gm/beta),0.0d0,0.0d0/) !Use this one!
!     v_r=(/1.0d0,1.0d0/gm,mach,0.0d0,0.0d0,dsqrt(2.d0/gm/beta),0.0d0,0.0d0/)  

	u(1,:,:,1)=f_p*rcom
	u(1,:,:,2)=-f_p*mach
	u(1,:,:,3)=0.0d0!-0.1d0*u(1,:,:,1) 
	u(1,:,:,4)=0.0d0 !vz momentum 
	u(1,:,:,5)=(f_p_p*rpres/gm)/(gm-1.d0)+0.5d0*(f_p*mach**2)/rcom +0.5d0*(2.d0/gm/beta)
	u(1,:,:,6)=dsqrt(2.d0/gm/beta) !bx
	u(1,:,:,7)=0.d0!u(1,:,:,7) !by
	u(1,:,:,8)=0.d0!u(1,:,:,8) !bz

	u(2,:,:,1)=f_p*rcom
	u(2,:,:,2)=-f_p*mach
	u(2,:,:,3)=0.0d0!-0.1d0*u(1,:,:,1) 
	u(2,:,:,4)=0.0d0 !vz momentum 
	u(2,:,:,5)=(f_p_p*rpres/gm)/(gm-1.d0)+0.5d0*(f_p*mach**2)/rcom +0.5d0*(2.d0/gm/beta)
	u(2,:,:,6)=dsqrt(2.d0/gm/beta) !bx
	u(2,:,:,7)=0.d0!u(1,:,:,7) !by
	u(2,:,:,8)=0.d0!u(1,:,:,8) !bz

	u(3,:,:,1)=f_p*rcom
	u(3,:,:,2)=-f_p*mach
	u(3,:,:,3)=0.0d0!-0.1d0*u(1,:,:,1) 
	u(3,:,:,4)=0.0d0 !vz momentum 
	u(3,:,:,5)=(f_p_p*rpres/gm)/(gm-1.d0)+0.5d0*(f_p*mach**2)/rcom +0.5d0*(2.d0/gm/beta)
	u(3,:,:,6)=dsqrt(2.d0/gm/beta) !bx
	u(3,:,:,7)=0.d0!u(1,:,:,7) !by
	u(3,:,:,8)=0.d0!u(1,:,:,8) !bz

	u(4,:,:,1)=f_p*rcom
	u(4,:,:,2)=-f_p*mach
	u(4,:,:,3)=0.0d0!-0.1d0*u(1,:,:,1) 
	u(4,:,:,4)=0.0d0 !vz momentum 
	u(4,:,:,5)=(f_p_p*rpres/gm)/(gm-1.d0)+0.5d0*(f_p*mach**2)/rcom +0.5d0*(2.d0/gm/beta)
	u(4,:,:,6)=dsqrt(2.d0/gm/beta) !bx
	u(4,:,:,7)=0.d0!u(1,:,:,7) !by
	u(4,:,:,8)=0.d0!u(1,:,:,8) !bz

endif

if (dir .eq. 1) then 
	mach=2.d0
	rcom=(gm+1.d0)*mach**2/(2.d0+(gm-1.d0)*mach**2)
	rpres=1.d0+gm*mach**2*(1.d0-1.d0/rcom)
!     v_l=(/rcom,rpres/gm,-mach/rcom,0.0d0,0.0d0,dsqrt(2.d0/gm/beta),0.0d0,0.0d0/) !Use this one!
!     v_r=(/1.0d0,1.0d0/gm,mach,0.0d0,0.0d0,dsqrt(2.d0/gm/beta),0.0d0,0.0d0/)  
!print*,rpres/gm
	u(ix,:,:,1)=1.0d0*f_p
	u(ix,:,:,2)=-f_p*mach
	u(ix,:,:,3)=0.0d0!-0.1d0*u(1,:,:,1) 
	u(ix,:,:,4)=0.0d0 !vz momentum 
	u(ix,:,:,5)=(f_p_p*1.d0/gm)/(gm-1.d0)+0.5d0*f_p*mach**2 +0.5d0*(2.d0/gm/beta)
	u(ix,:,:,6)=dsqrt(2.d0/gm/beta) !bx
	u(ix,:,:,7)=0.d0!u(1,:,:,7) !by
	u(ix,:,:,8)=0.d0!u(1,:,:,8) !bz

	u(ix-1,:,:,1)=f_p*1.0d0
	u(ix-1,:,:,2)=-f_p*mach
	u(ix-1,:,:,3)=0.0d0!-0.1d0*u(1,:,:,1) 
	u(ix-1,:,:,4)=0.0d0 !vz momentum 
	u(ix-1,:,:,5)=(f_p_p*1.d0/gm)/(gm-1.d0)+0.5d0*f_p*mach**2+0.5d0*(2.d0/gm/beta)
	u(ix-1,:,:,6)=dsqrt(2.d0/gm/beta) !bx
	u(ix-1,:,:,7)=0.d0!u(1,:,:,7) !by
	u(ix-1,:,:,8)=0.d0!u(1,:,:,8) !bz

	u(ix-2,:,:,1)=f_p*1.0d0
	u(ix-2,:,:,2)=-f_p*mach
	u(ix-2,:,:,3)=0.0d0!-0.1d0*u(1,:,:,1) 
	u(ix-2,:,:,4)=0.0d0 !vz momentum 
	u(ix-2,:,:,5)=(f_p_p*1.d0/gm)/(gm-1.d0)+0.5d0*f_p*mach**2+0.5d0*(2.d0/gm/beta)
	u(ix-2,:,:,6)=dsqrt(2.d0/gm/beta) !bx
	u(ix-2,:,:,7)=0.d0!u(1,:,:,7) !by
	u(ix-2,:,:,8)=0.d0!u(1,:,:,8) !bz

	u(ix-3,:,:,1)=f_p*1.0d0
	u(ix-3,:,:,2)=-f_p*mach
	u(ix-3,:,:,3)=0.0d0!-0.1d0*u(1,:,:,1) 
	u(ix-3,:,:,4)=0.0d0 !vz momentum 
	u(ix-3,:,:,5)=(f_p_p*1.d0/gm)/(gm-1.d0)+0.5d0*f_p*mach**2+0.5d0*(2.d0/gm/beta)
	u(ix-3,:,:,6)=dsqrt(2.d0/gm/beta) !bx
	u(ix-3,:,:,7)=0.d0!u(1,:,:,7) !by
	u(ix-3,:,:,8)=0.d0!u(1,:,:,8) !bz
!!!!!!!!!!!!!!!!!!!!!!!!!!!
endif
end subroutine boundary_control_custom_mhd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module boundary_rot
