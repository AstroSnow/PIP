module boundary_rot
!====================================================================
! This is the boundary module for the PIP code.
! Module for calling the different boundaries
! Author A.Hillier
! First version - 2013/05/31 (Fixed routines not completed)
! Second version - 2013/06/27 (New labelling for boundaries)
!====================================================================
  use globalvar,only:ix,jx,kx,ndim,flag_pip,flag_mhd,nvar_h,nvar_m,&
       flag_mpi,flag_bnd,margin,neighbor,flag_divb,time
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

  subroutine PIPbnd(U_h,U_m)
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)

    if(flag_mpi.eq.1) then
!       print *,"MY_RANK",my_rank
!       if(my_rank.eq.0) then
!          print *,my_rank,U_m(4,303:306,1,6),"BEF"
!       else
!          print *,my_rank,U_m(4,3:6,1,6),"BEF"
!       endif
       call mpi_bnd(U_h,U_m)
!       if(my_rank.eq.0) then
!          print *,my_rank,U_m(4,303:306,1,6),"AFT"
!       else
!          print *,my_rank,U_m(4,3:6,1,6),"AFT"
!       endif
!       stop
    endif    
    if (flag_pip.eq.1) then       
       call MHDbnd(U_m)   
       call HDbnd(U_h)   
    elseif (flag_mhd.eq.1) then       
       call MHDbnd(U_m)          
    else       
       call HDbnd(U_h)          
    endif    



end subroutine PIPbnd

  subroutine MHDbnd(U_m)   
    double precision,intent(inout)::U_m(ix,jx,kx,nvar_m)
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
!       else if(flag_bnd(i).eq.20) then
!          type=2
       else
          type=flag_bnd(i)
       endif

       if(neighbor(i).eq.-1) then
          call boundary_control(U_m,ix,jx,kx,nvar_m,margin,dir,upper_lower, &
               sym_mhd(:,i),type)
       endif
    enddo
  end subroutine MHDbnd

  subroutine HDbnd(U_h)   
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h)
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
          call boundary_control(U_h,&
               ix,jx,kx,nvar_h,margin,dir,upper_lower, &
               sym_hd(:,i),type)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Velocity driver
       if (type .EQ. 20) then 
          if(dir.eq.1) then 
             do n=1,margin(1)
                u(1,:,:,1)=u(2,:,:,1) !Density
!                u(2,:,:,2)=u(2,:,:,2)+u(2,:,:,1)*0.1d0*dsin(time*2.0d0*pi/30.0d0) !vx momentum
                u(2,:,:,2)=u(2,:,:,2)+u(2,:,:,1)*1.0e-7*dcos(time*2.0d0*pi/1.0d0) !vx momentum
                u(1,:,:,2)=u(2,:,:,2)
                u(1,:,:,3)=u(2,:,:,3) !vy momentum 
                u(1,:,:,4)=u(2,:,:,4) !vz momentum 
!                u(2,:,:,5)=u(2,:,:,5) +0.5d0*u(2,:,:,1)* (0.1d0*dsin(time*2.0d0*pi/30.0d0))**2.0d0 !energy
                u(2,:,:,5)=u(2,:,:,5) +0.5d0*u(2,:,:,1)* (1.0e-7*dcos(time*2.0d0*pi/1.0d0))**2.0d0 !energy
                 u(1,:,:,5)=u(2,:,:,5) 
		!print*,sl(1)+n-1,sr(1),u(sl(1)+n-1,:,:,2),u(sl(1)+n-1,:,:,5),0.1d0*dsin(time*2.0d0*pi/30.0d0),time
             enddo
	  endif
       endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Velocity driver
       if (type .EQ. 20) then 
	print*,'onevar type=',type
	stop
          if(dir.eq.1) then 
             do n=1,margin(1)
                u1(sl(1)+n-1,:,:)=u1(sr(1),:,:)+0.5d0*u1(sr(1),:,:)*(100.0d0*sin(time*1000.0d0))**2.d0 !KINETIC ENERGY?!!!
		print*,u1(sl(1)+n-1,:,:)
             enddo
          endif
       endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif
  end subroutine boundary_control_onevar

  
end module boundary_rot
