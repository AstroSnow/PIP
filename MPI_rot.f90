module mpi_rot
  use globalvar,only:total_prc,my_rank,nvar_h,nvar_m,mpi_siz,mpi_pos,margin,&
       ix,jx,kx,ndim,flag_mpi_split,neighbor,U_h,U_m,flag_mhd,flag_pip,&
       cno,vmpi,flag_mpi,flag_bnd, &
       dimsFile, ig
  implicit none
  include "mpif.h"
  integer ierr
  integer status(6)
contains
  subroutine set_mpi
    integer i
    !charm for MPI
    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD,total_prc,ierr)
    call mpi_comm_rank(MPI_COMM_WORLD,my_rank,ierr)

    if(flag_mpi.eq.0 .or.my_rank.eq.0) print *,"total cpu",total_prc
    vmpi=nvar_h+nvar_m
    write(cno,"(i4.4)")my_rank
    !   print *,"TOT",total_prc,my_rank,cno
    if(total_prc.eq.1) then
      mpi_siz(:)=1
      mpi_pos(:)=0
    else
      call set_mpi_pos
    endif

    ! store full ranges + margins for use in write-IO
    dimsFile(1:3) = [ix, jx, kx]
    ! compute # of grid-points for individual process
    do i=1,3
      ig(i) = (dimsFile(i)-2*margin(i))/mpi_siz(i)+2*margin(i)
      ! Add excess grid-points from uneven splitting to last block along dimension
      if (mpi_pos(i) .eq. (mpi_siz(i)-1)) then
        ig(i) = ig(i) + mod(dimsFile(i)-2*margin(i), mpi_siz(i))
      end if
    end do
    ! replace _x range values with those for the process
    ix=ig(1); jx=ig(2); kx=ig(3)
  end subroutine set_mpi

  subroutine set_mpi_pos
    integer n0,tmp!,tmp2,n2,n3
    integer n!ixy,iyz,izx
    if(mpi_siz(1)*mpi_siz(2)*mpi_siz(3).eq.0  &
         .or. mpi_siz(1)*mpi_siz(2)*mpi_siz(3).ne.total_prc) then
       select case(flag_mpi_split)
       !1D domain decomposition
       case(1)
          mpi_siz(1)=total_prc
          mpi_siz(2)=1
          mpi_siz(3)=1
       !n-dimensional decomposition
       case(2)
          call set_domain_equally(mpi_siz(1),mpi_siz(2),mpi_siz(3),3,total_prc)
       !adjust
       case(3)
          n0=0
          do n=1,3
             if(mpi_siz(n).eq.0) n0=n0+1
          enddo
          select case(n0)
          case(1)
             do n=1,3
                if(mpi_siz(n).eq.0) mpi_siz(n)=total_prc/ &
                     (mpi_siz(mod(n,3)+1)*mpi_siz(mod(n+1,3)+1))
             enddo
          case(2)
             do n=1,3
                call set_domain_equally(mpi_siz(mod(n,3)+1), &
                     mpi_siz(mod(n+1,3)+1),tmp,2,total_prc/mpi_siz(n))
             enddo
          case(3)
             call set_domain_equally(mpi_siz(1),mpi_siz(2),mpi_siz(3), &
                  3,total_prc)
          end select                    
       end select
    endif    
    mpi_pos(1)=mod(my_rank,mpi_siz(1))
    mpi_pos(2)=mod(my_rank/mpi_siz(1),mpi_siz(2))
    mpi_pos(3)=my_rank/(mpi_siz(1)*mpi_siz(2))

    if(flag_mpi.eq.1) then
       if(my_rank.eq.0) print *,"siz",mpi_siz    
    endif
   
  end subroutine set_mpi_pos

  subroutine set_mpi_neighbor
    integer nn,offsets(3),n
    neighbor(1)=my_rank-1
    neighbor(2)=my_rank+1
    neighbor(3)=my_rank-mpi_siz(1)
    neighbor(4)=my_rank+mpi_siz(1)
    neighbor(5)=my_rank-mpi_siz(1)*mpi_siz(2)
    neighbor(6)=my_rank+mpi_siz(1)*mpi_siz(2)
    
    !connect the other side of domain if periodic boundary is set
    !else set neighbor -1 (numerical boundary)
    offsets(1)=1
    offsets(2)=mpi_siz(1)
    offsets(3)=mpi_siz(1)*mpi_siz(2)
    do n=1,3
       if(mpi_pos(n).eq.0) then
          nn=2*n-1
          if(flag_bnd(nn).eq.1) then
             neighbor(nn)=my_rank+(mpi_siz(n)-1)*offsets(n)
          else
             neighbor(nn)=-1
          endif
       endif
       if(mpi_pos(n).eq.(mpi_siz(n)-1)) then
          nn=2*n
          if(flag_bnd(nn).eq.1) then
             neighbor(nn)=my_rank-(mpi_siz(n)-1)*offsets(n)
          else
             neighbor(nn)=-1
          endif
       endif       
    enddo
    do n=1,3
       if(mpi_siz(n).eq.1) neighbor(2*n-1:2*n)=-1
    enddo

  end subroutine set_mpi_neighbor

  subroutine mpi_bnd(U_h,U_m)
    integer n0,nsend,nrecv,ncount
!    integer xcount,ycount,zcount
    double precision,intent(inout)::U_h(ix,jx,kx,nvar_h),U_m(ix,jx,kx,nvar_m)
    double precision sendx(margin(1),jx,kx,vmpi),recvx(margin(1),jx,kx,vmpi)   
    double precision sendy(ix,margin(2),kx,vmpi),recvy(ix,margin(2),kx,vmpi)   
    double precision sendz(ix,jx,margin(3),vmpi),recvz(ix,jx,margin(3),vmpi)   
    
!for x-direction

    if(mpi_siz(1).gt.1) then
       !=============================================================
       !left-boundary of x-direction
       n0=1
       if(flag_pip.eq.1 .or. flag_mhd.eq.0) then
          sendx(:,:,:,1:nvar_h)=U_h(ix-2*margin(1)+1:ix-margin(1),:,:,:)
          n0=n0+nvar_h
       endif
       if(flag_mhd.eq.1) then
          sendx(:,:,:,n0:n0+nvar_m-1)=U_m(ix-2*margin(1)+1:ix-margin(1),:,:,:)
          n0=n0+nvar_m
       endif
       ncount=kx*jx*margin(1)*(n0-1)

       nsend=neighbor(2)
       nrecv=neighbor(1)       
       if(nsend.eq.-1) nsend=mpi_proc_null
       if(nrecv.eq.-1) nrecv=mpi_proc_null
 
       call mpi_SendRecv(sendx,ncount,mpi_double_precision,nsend,0, &
            recvx,ncount,mpi_double_precision,nrecv,0,mpi_comm_world,status,ierr)       
       if(neighbor(1).ne.-1) then
          n0=1
          if(flag_pip.eq.1.or.flag_mhd.eq.0) then
             U_h(1:margin(1),:,:,:)=recvx(:,:,:,n0:n0+nvar_h-1)
             n0=n0+nvar_h
          endif
         
          if(flag_mhd.eq.1) then
             U_m(1:margin(1),:,:,:)=recvx(:,:,:,n0:n0+nvar_m-1)
             n0=n0+nvar_m
          endif          
       endif
       !=============================================================
       !right-boundary of x-direction
       n0=1
       if(flag_pip.eq.1 .or. flag_mhd.eq.0) then
          sendx(:,:,:,1:nvar_h)=U_h(margin(1)+1:2*margin(1),:,:,:)
          n0=n0+nvar_h
       endif
       if(flag_mhd.eq.1) then
          sendx(:,:,:,n0:n0+nvar_m-1)=U_m(margin(1)+1:2*margin(1),:,:,:)
          n0=n0+nvar_m
       endif
       ncount=kx*jx*margin(1)*(n0-1)

       nsend=neighbor(1)
       nrecv=neighbor(2)       
       if(nsend.eq.-1) nsend=mpi_proc_null
       if(nrecv.eq.-1) nrecv=mpi_proc_null

       call mpi_SendRecv(sendx,ncount,mpi_double_precision,nsend,0, &
            recvx,ncount,mpi_double_precision,nrecv,0,mpi_comm_world,status,ierr)

       if(neighbor(2).ne.-1) then
          n0=1
          if(flag_pip.eq.1.or.flag_mhd.eq.0) then
             U_h(ix-margin(1)+1:ix,:,:,:)=recvx(:,:,:,n0:n0+nvar_h-1)
             n0=n0+nvar_h
          endif
          if(flag_mhd.eq.1) then
             U_m(ix-margin(1)+1:ix,:,:,:)=recvx(:,:,:,n0:n0+nvar_m-1)
             n0=n0+nvar_m
          endif          
       endif
    endif


!for y-direction
    if(mpi_siz(2).gt.1) then
       !=============================================================
       !left-boundary of y-direction
       n0=1
       if(flag_pip.eq.1 .or. flag_mhd.eq.0) then
          sendy(:,:,:,1:nvar_h)      =U_h(:,jx-2*margin(2)+1:jx-margin(2),:,:)
          n0=n0+nvar_h
       endif
       if(flag_mhd.eq.1) then
          sendy(:,:,:,n0:n0+nvar_m-1)=U_m(:,jx-2*margin(2)+1:jx-margin(2),:,:)
          n0=n0+nvar_m
       endif
       ncount=kx*ix*margin(2)*(n0-1)

       nsend=neighbor(4)
       nrecv=neighbor(3)       
       if(nsend.eq.-1) nsend=mpi_proc_null
       if(nrecv.eq.-1) nrecv=mpi_proc_null
       call mpi_SendRecv(sendy,ncount,mpi_double_precision,nsend,0, &
            recvy,ncount,mpi_double_precision,nrecv,0,mpi_comm_world,status,ierr)
       
       if(neighbor(3).ne.-1) then
          n0=1
          if(flag_pip.eq.1.or.flag_mhd.eq.0) then
             U_h(:,1:margin(2),:,:)=recvy(:,:,:,n0:n0+nvar_h-1)
             n0=n0+nvar_h
          endif
          if(flag_mhd.eq.1) then
             U_m(:,1:margin(2),:,:)=recvy(:,:,:,n0:n0+nvar_m-1)
             n0=n0+nvar_m
          endif          
       endif


       !=============================================================
       !right-boundary of y-direction
       n0=1
       if(flag_pip.eq.1 .or. flag_mhd.eq.0) then
          sendy(:,:,:,1:nvar_h)=U_h(:,margin(2)+1:2*margin(2),:,:)
          n0=n0+nvar_h
       endif
       if(flag_mhd.eq.1) then
          sendy(:,:,:,n0:n0+nvar_m-1)=U_m(:,margin(2)+1:2*margin(2),:,:)
          n0=n0+nvar_m
       endif
       ncount=kx*ix*margin(2)*(n0-1)

       nsend=neighbor(3)
       nrecv=neighbor(4)       
       if(nsend.eq.-1) nsend=mpi_proc_null
       if(nrecv.eq.-1) nrecv=mpi_proc_null
       call mpi_SendRecv(sendy,ncount,mpi_double_precision,nsend,0, &
            recvy,ncount,mpi_double_precision,nrecv,0,mpi_comm_world,status,ierr)       
       if(neighbor(4).ne.-1) then
          n0=1
          if(flag_pip.eq.1.or.flag_mhd.eq.0) then
             U_h(:,jx-margin(2)+1:jx,:,:)=recvy(:,:,:,n0:n0+nvar_h-1)
             n0=n0+nvar_h
          endif
          if(flag_mhd.eq.1) then
             U_m(:,jx-margin(2)+1:jx,:,:)=recvy(:,:,:,n0:n0+nvar_m-1)
             n0=n0+nvar_m
          endif          

       endif

    endif
    
    !for z-direction
    if(mpi_siz(3).gt.1) then
       !=============================================================
       !left-boundary of z-direction
       n0=1
       if(flag_pip.eq.1 .or. flag_mhd.eq.0) then
          sendz(:,:,:,1:nvar_h)=U_h(:,:,kx-2*margin(3)+1:kx-margin(3),:)
          n0=n0+nvar_h
       endif
       if(flag_mhd.eq.1) then
          sendz(:,:,:,n0:n0+nvar_m-1)=U_m(:,:,kx-2*margin(3)+1:kx-margin(3),:)
          n0=n0+nvar_m
       endif
       ncount=jx*ix*margin(3)*(n0-1)

       nsend=neighbor(6)
       nrecv=neighbor(5)       
       if(nsend.eq.-1) nsend=mpi_proc_null
       if(nrecv.eq.-1) nrecv=mpi_proc_null
       call mpi_SendRecv(sendz,ncount,mpi_double_precision,nsend,0, &
            recvz,ncount,mpi_double_precision,nrecv,0,mpi_comm_world,status,ierr)       
       if(neighbor(5).ne.-1) then
          n0=1
          if(flag_pip.eq.1.or.flag_mhd.eq.0) then
             U_h(:,:,1:margin(3),:)=recvz(:,:,:,n0:n0+nvar_h-1)
             n0=n0+nvar_h
          endif
          if(flag_mhd.eq.1) then
             U_m(:,:,1:margin(3),:)=recvz(:,:,:,n0:n0+nvar_m-1)
             n0=n0+nvar_m
          endif          
       endif

       !=============================================================
       !right-boundary of z-direction
       n0=1
       if(flag_pip.eq.1 .or. flag_mhd.eq.0) then
          sendz(:,:,:,1:nvar_h)=U_h(:,:,margin(3)+1:2*margin(3),:)
          n0=n0+nvar_h
       endif
       if(flag_mhd.eq.1) then
          sendz(:,:,:,n0:n0+nvar_m-1)=U_m(:,:,margin(3)+1:2*margin(3),:)
          n0=n0+nvar_m
       endif
       ncount=jx*ix*margin(3)*(n0-1)

       nsend=neighbor(5)
       nrecv=neighbor(6)       
       if(nsend.eq.-1) nsend=mpi_proc_null
       if(nrecv.eq.-1) nrecv=mpi_proc_null
       call mpi_SendRecv(sendz,ncount,mpi_double_precision,nsend,0, &
            recvz,ncount,mpi_double_precision,nrecv,0,mpi_comm_world,status,ierr)       
       if(neighbor(6).ne.-1) then
          n0=1
          if(flag_pip.eq.1.or.flag_mhd.eq.0) then
             U_h(:,:,kx-margin(3)+1:kx,:)=recvz(:,:,:,n0:n0+nvar_h-1)
             n0=n0+nvar_h
          endif
          if(flag_mhd.eq.1) then
             U_m(:,:,kx-margin(3)+1:kx,:)=recvz(:,:,:,n0:n0+nvar_m-1)
             n0=n0+nvar_m
          endif          
       endif
    endif
  end subroutine mpi_bnd



  subroutine mpi_bnd_onevar(u1)
    integer n0,nsend,nrecv,ncount
!    integer xcount,ycount,zcount
    double precision,intent(inout)::u1(ix,jx,kx)
    double precision sendx(margin(1),jx,kx),recvx(margin(1),jx,kx)   
    double precision sendy(ix,margin(2),kx),recvy(ix,margin(2),kx)   
    double precision sendz(ix,jx,margin(3)),recvz(ix,jx,margin(3))   
    
!for x-direction

    if(mpi_siz(1).gt.1) then
       !=============================================================
       !left-boundary of x-direction
       sendx(:,:,:)=u1(ix-2*margin(1)+1:ix-margin(1),:,:)
       ncount=kx*jx*margin(1)

       nsend=neighbor(2)
       nrecv=neighbor(1)       
       if(nsend.eq.-1) nsend=mpi_proc_null
       if(nrecv.eq.-1) nrecv=mpi_proc_null
 
       call mpi_SendRecv(sendx,ncount,mpi_double_precision,nsend,0, &
            recvx,ncount,mpi_double_precision,nrecv,0,mpi_comm_world,status,ierr)       
       if(neighbor(1).ne.-1) then
          n0=1
          u1(1:margin(1),:,:)=recvx(:,:,:)
       endif
       !=============================================================
       !right-boundary of x-direction
       sendx(:,:,:)=u1(margin(1)+1:2*margin(1),:,:)
       ncount=kx*jx*margin(1)

       nsend=neighbor(1)
       nrecv=neighbor(2)       
       if(nsend.eq.-1) nsend=mpi_proc_null
       if(nrecv.eq.-1) nrecv=mpi_proc_null

       call mpi_SendRecv(sendx,ncount,mpi_double_precision,nsend,0, &
            recvx,ncount,mpi_double_precision,nrecv,0,mpi_comm_world,status,ierr)

       if(neighbor(2).ne.-1) then
          u1(ix-margin(1)+1:ix,:,:)=recvx(:,:,:)
       endif
    endif


!for y-direction
    if(mpi_siz(2).gt.1) then
       !=============================================================
       !left-boundary of y-direction
       sendy(:,:,:)      =u1(:,jx-2*margin(2)+1:jx-margin(2),:)
       ncount=kx*ix*margin(2)

       nsend=neighbor(4)
       nrecv=neighbor(3)       
       if(nsend.eq.-1) nsend=mpi_proc_null
       if(nrecv.eq.-1) nrecv=mpi_proc_null
       call mpi_SendRecv(sendy,ncount,mpi_double_precision,nsend,0, &
            recvy,ncount,mpi_double_precision,nrecv,0,mpi_comm_world,status,ierr)
       
       if(neighbor(3).ne.-1) then
          u1(:,1:margin(2),:)=recvy(:,:,:)
       endif


       !=============================================================
       !right-boundary of y-direction
       sendy(:,:,:)=u1(:,margin(2)+1:2*margin(2),:)
       ncount=kx*ix*margin(2)

       nsend=neighbor(3)
       nrecv=neighbor(4)       
       if(nsend.eq.-1) nsend=mpi_proc_null
       if(nrecv.eq.-1) nrecv=mpi_proc_null
       call mpi_SendRecv(sendy,ncount,mpi_double_precision,nsend,0, &
            recvy,ncount,mpi_double_precision,nrecv,0,mpi_comm_world,status,ierr)       
       if(neighbor(4).ne.-1) then
          u1(:,jx-margin(2)+1:jx,:)=recvy(:,:,:)
       endif

    endif
    
    !for z-direction
    if(mpi_siz(3).gt.1) then
       !=============================================================
       !left-boundary of z-direction
       sendz(:,:,:)=u1(:,:,kx-2*margin(3)+1:kx-margin(3))
       ncount=jx*ix*margin(3)

       nsend=neighbor(6)
       nrecv=neighbor(5)       
       if(nsend.eq.-1) nsend=mpi_proc_null
       if(nrecv.eq.-1) nrecv=mpi_proc_null
       call mpi_SendRecv(sendz,ncount,mpi_double_precision,nsend,0, &
            recvz,ncount,mpi_double_precision,nrecv,0,mpi_comm_world,status,ierr)       
       if(neighbor(5).ne.-1) then
          u1(:,:,1:margin(3))=recvz(:,:,:)
       endif

       !=============================================================
       !right-boundary of z-direction
       sendz(:,:,:)=u1(:,:,margin(3)+1:2*margin(3))
       ncount=jx*ix*margin(3)

       nsend=neighbor(5)
       nrecv=neighbor(6)       
       if(nsend.eq.-1) nsend=mpi_proc_null
       if(nrecv.eq.-1) nrecv=mpi_proc_null
       call mpi_SendRecv(sendz,ncount,mpi_double_precision,nsend,0, &
            recvz,ncount,mpi_double_precision,nrecv,0,mpi_comm_world,status,ierr)       
       if(neighbor(6).ne.-1) then
          u1(:,:,kx-margin(3)+1:kx)=recvz(:,:,:)
       endif
    endif
  end subroutine mpi_bnd_onevar

  
  subroutine set_domain_equally(mpi1,mpi2,mpi3,ndim,total)
    integer,intent(in)::total,ndim
    integer,intent(inout)::mpi1,mpi2,mpi3
    integer n2,n3,tmp
    tmp=total
    !prime decomposition for 2 and 3
    n2=0
    do while(mod(tmp,2).eq.0.and.tmp.gt.1)
       tmp=int(tmp/2)
       n2=n2+1
    end do
    n3=0
    do while(mod(tmp,3).eq.0.and.tmp.gt.1)
       tmp=int(tmp/3)
       n3=n3+1
    end do
    mpi1=2**(n2/ndim)*3**(n3/ndim)
    if (ndim.eq.2) then
       mpi2=total_prc/mpi1
    else if(ndim.eq.3) then
       mpi2=2**(n2/ndim)*3**(n3/ndim)
       mpi3=total_prc/(mpi1*mpi2)               
    endif
  end subroutine set_domain_equally

  function mpi_double_interface(data,type)
    integer,intent(in)::type
    integer MPI_TYPE
    double precision,intent(in)::data
    double precision mpi_double_interface
    select case(type)
    case(1)
       MPI_TYPE=MPI_MIN
    case(2)
       MPI_TYPE=MPI_SUM
    case(3)
       MPI_TYPE=MPI_MAX
    end select
    call mpi_allreduce(data,mpi_double_interface,1,mpi_double_precision,MPI_TYPE, &
         mpi_comm_world,ierr)
  end function mpi_double_interface

  subroutine end_mpi
    call mpi_finalize(ierr)
  end subroutine end_mpi

  subroutine my_mpi_barrier
    call MPI_barrier(MPI_COMM_WORLD,ierr)
  end subroutine my_mpi_barrier  

  subroutine send_Gm_rec_ref(Gm_rec_ref)
    double precision, intent(in):: Gm_rec_ref
    !Send the normalised value to all the processors
    call MPI_Bcast(Gm_rec_ref,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!    call MPI_Bcast(ierr,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
  end subroutine send_Gm_rec_ref

end module mpi_rot
