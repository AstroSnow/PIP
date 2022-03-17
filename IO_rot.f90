module io_rot
!====================================================================
! This is the in_out module for the PIP code.
! Module for input and output
! Author A.Hillier
! First version - 2013/04/26 NN
! Modification history
! 2013/05/29  S.Takasao
!====================================================================
  use globalvar,only:ix,jx,kx,ndim,dt,dt_cnd,U_h,U_m,time,cno,flag_stop,&
       flag_mhd,flag_pip,flag_mpi,my_rank,indir,outdir,&
       ac,gm_ion,gm_rec,gra,eta,dxc,dyc,dzc,dx,dy,dz,x,y,z,&
       flag_ir,nvar_h,nvar_m,flag_resi,nt,nout,margin,gm,flag_restart,&
       flag_bnd,flag_col,flag_grav,tend,mpi_pos,xi_n,mu,flag_visc,&
       total_iter,flag_amb,dtout,mpi_siz,nt,nmax,output_type,flag_ps,flag_divb,&
       flag_damp,damp_time,flag_rad,flag_ir_type,arb_heat,visc,esav,emsavtime,&
       ac_sav, xi_sav, ion_sav, rec_sav, col_sav, gr_sav, vs_sav, heat_sav, et_sav, ps_sav,&
       Nexcite, &
       file_id, plist_id, hdf5_error, filespace_id, memspace_id,&
       dimsFile, dimsMem, start_stop, hdf5_offset, neighbor, ig
  use mpi_rot,only:end_mpi
  use IOT_rot,only:initialize_IOT,get_next_output
  use Util_rot,only:get_word,get_value_integer
  use HDF5
  implicit none
  include "mpif.h"
  integer ios
  integer,allocatable,save::mf_m(:,:),mf_h(:,:)
!  integer,save::mf_m(8,2),mf_h(5,2)
  character*15,allocatable::file_m(:),file_h(:)
  character*4 tno
  integer, parameter :: mf_t=10 &
                       ,mf_x=11, mf_y=12, mf_z=13 &
                       ,mf_dx=14, mf_dy=15, mf_dz=16
  ! version number (date)
  integer, parameter :: mf_info=9, version=20140708,restart_unit=77
  double precision start_time,end_time
contains
!output subroutine called from main
!  subroutine output(nout,time)
  subroutine output(out)
    integer,intent(in)::out
    integer i,j,k,outesav
!    integer j
    double precision total_divB,cx,cy,max_C,divb
!    double precision,intent(in)::time
    !if nout = 0 initial setting for output is done^^^^^^^^^^^^^^^^^^^^^^
    nt=nt+1
    if(nout.eq.0) then
       call set_initial_out
       call save_coordinates
       call def_varfiles(0)
       start_time=MPI_Wtime()
       if(flag_mpi.eq.0 .or. my_rank.eq.0) then
          call mk_config
       endif
       call initialize_IOT(dtout,tend,output_type)
    endif
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!Check for physical-time save
end_time=MPI_Wtime()
outesav=-1
if (((end_time-start_time) .ge. emsavtime) .and. (esav .eq. 1)) then
	outesav=1
	esav=2
     if(flag_mpi.eq.0 .or.my_rank.eq.0) &
          write(6,*) 'calling emergency save'
endif
    !output variables
  if((time.ge.get_next_output(nout,time,esav)) .or. (out.eq.1) .or. (outesav.eq.1)) then
     !     if(flag_mpi.eq.0 .or.my_rank.eq.0) write(6,*) 'Time,dt,nout,nt,total_iter: ',time,dt,nout,nt,total_iter
    end_time=MPI_Wtime()
     if(flag_mpi.eq.0 .or.my_rank.eq.0) &
!          write(6,*) 'Time,dt,dt_cnd,nout,nt,total_iter: ',time,dt,dt_cnd,nout,nt,total_iter
          write(6,*) 'Time,dt,nout,nt,elapsed time: ',time,dt,nout,nt,end_time-start_time
       total_divb=0.0d0
       max_C=0.0d0
       if(ndim.eq.100) then
          do j=margin(2)+1,jx-margin(2);do i=margin(1)+1,ix-margin(1)
             divb=abs((-U_m(i+2,j,1,6)+8.0d0*U_m(i+1,j,1,6)&
                  -8.0d0*U_m(i-1,j,1,6)+U_m(i-2,j,1,6))/dx(i) +&
                  (-U_m(i,j+2,1,7)+8.0d0*U_m(i,j+1,1,7)&
                  -8.0d0*U_m(i,j-1,1,7)+U_m(i,j-2,1,7))/dy(j))

             total_divb=total_divb +divb
!             if(divb.gt.1.0d-13) then
!                print *,x(i),y(j),divb
!             endif
             max_C=max(max_C,((-U_m(i+2,j,1,7)+8.0d0*U_m(i+1,j,1,7) &
                  -8.0d0*U_m(i-1,j,1,7)+U_m(i-2,j,1,7))/(12.0d0*dx(i)))**2+ &
                  ((-U_m(i,j+2,1,6)+8.0d0*U_m(i,j+1,1,6) &
                  -8.0d0*U_m(i,j-1,1,6)+U_m(i,j-2,1,6))/(12.0d0*dy(j)))**2)
          enddo;enddo
          print *,"NT,TOTAL_DIVB, maxJ =",nt,total_divb,sqrt(max_C)
       endif
       call save_varfiles(nout)
       nout=nout+1
    endif

  end subroutine output

  subroutine save_coordinates  
    character*4 tmp_id

    if(flag_mpi.eq.0 .or.(mpi_pos(2).eq.0.and.mpi_pos(3).eq.0)) then
       write(tmp_id,"(i4.4)")mpi_pos(1)
       call dacdef1d(mf_x,trim(outdir) // 'x.dac.'//tmp_id,6,ix)
       write(mf_x) x
       call dacdef1d(mf_dx,trim(outdir) // 'dx.dac.'//tmp_id,6,ix)
       write(mf_dx) dx
       close(mf_x)
       close(mf_dx)
    endif

    if(ndim.ge.2) then
       if(flag_mpi.eq.0 .or.(mpi_pos(1).eq.0.and.mpi_pos(3).eq.0)) then
          write(tmp_id,"(i4.4)")mpi_pos(2)
          call dacdef1d(mf_y,trim(outdir) // 'y.dac.'//tmp_id,6,jx)
          write(mf_y) y
          call dacdef1d(mf_dy,trim(outdir) // 'dy.dac.'//tmp_id,6,jx)
          write(mf_dy) dy
          close(mf_y)
          close(mf_dy)

       endif
       if(ndim.ge.3) then
          if(flag_mpi.eq.0 .or.(mpi_pos(1).eq.0.and.mpi_pos(2).eq.0)) then
             write(tmp_id,"(i4.4)")mpi_pos(3)
             call dacdef1d(mf_z,trim(outdir) // 'z.dac.'//tmp_id,6,kx)
             write(mf_z) z
             call dacdef1d(mf_dz,trim(outdir) // 'dz.dac.'//tmp_id,6,kx)
             write(mf_dz) dz
             close(mf_z)
             close(mf_dz)
          endif
       endif
    endif

  end subroutine save_coordinates

  subroutine write_1D_array(n, varname, data_array)
    integer(HID_T) :: dset_id
    integer(HSSIZE_T) :: offset(1)
    integer(HSIZE_T) :: dims_1D(1)
    integer :: n
    character(*) :: varname
    double precision, dimension(1) :: data_array(ig(n))

    ! Creating dataset
    CALL h5dcreate_f(file_id, trim(varname), H5T_NATIVE_DOUBLE, filespace_id(n), &
                     dset_id,  hdf5_error)
    ! Select hyperslab in the file.
    CALL h5dget_space_f(dset_id, filespace_id(n), hdf5_error)
    offset(1) = hdf5_offset(n)
    dims_1D(1) = dimsMem(n)
    CALL h5sselect_hyperslab_f(filespace_id(n), H5S_SELECT_SET_F, offset, dims_1D, hdf5_error)
    ! write data to file
    dims_1D(1) = dimsFile(n)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data_array(start_stop(n,1):start_stop(n,2)), &
                    dims_1D, hdf5_error, &
                    file_space_id=filespace_id(n), mem_space_id=memspace_id(n), &
                    xfer_prp=plist_id)
    ! Closing dataset connections
    CALL h5dclose_f(dset_id, hdf5_error)
    end subroutine write_1D_array

  subroutine def_varfiles(append)
    integer,intent(in)::append
    integer i

    write(tno,"(i4.4)")nout


    if(flag_pip.eq.1.or.flag_amb.eq.1) then
!print*,'ac_sav',ac_sav
!print*,'xi_sav',xi_sav
       if(ac_sav.eq.0) call save1param(ac,tno//'ac.dac.',1)
       if(xi_sav.eq.0) call save1param(xi_n,tno//'xi.dac.',1)
    endif
    if(flag_pip.eq.1.and.flag_ir.eq.1) then
       if(ion_sav.eq.0) call save1param(Gm_ion,tno//'ion.dac.',1)
       if(rec_sav.eq.0) call save1param(Gm_rec,tno//'rec.dac.',1)
    endif
    if(flag_mhd.eq.1.and.flag_resi.eq.1) then
       if(et_sav.eq.0) call save1param(eta,tno//'et.dac.',1)
    endif
    if(flag_col.eq.1) then
       if(col_sav.eq.0) call save1param(ac,tno//'col.dac.',1)
    endif

    if(flag_grav.eq.1) then
       if(gr_sav.eq.0) call save1param(gra,tno//'gr.dac.',3)
    endif
    if(flag_visc.eq.1) then
       if(vs_sav.eq.0) call save1param(mu,tno//'vs.dac.',1)
    endif
    if(flag_pip.eq.1.and.flag_ir_type.eq.0.and.flag_IR.ne.0) then
       if(heat_sav.eq.0) call save1param(arb_heat,tno//'aheat.dac.',1)
    endif


    if(flag_mpi.eq.0 .or.my_rank.eq.0)      &
         call dacdef0s(mf_t,trim(outdir) // 't.dac.'//cno,6,append)

  end subroutine def_varfiles


!   subroutine def_varfiles(append)
!     integer,intent(in)::append
!     integer i

!     if(flag_mpi.eq.0 .or.my_rank.eq.0)      &
!          call dacdef0s(mf_t,trim(outdir) // '/t.dac.'//cno,6,append)
!     if(flag_mhd.eq.1) then
! !       do i=1,nvar_m
! !       print *,"OK",my_rank,file_m(1)//cno
!        do i=1,8
!           call dacdef3s(mf_m(i,1),trim(outdir)//trim(file_m(i))//cno,6,append)
!        enddo
!        if(flag_resi.ge.2) then
!           call dacdef3s(77,trim(outdir)//'/et.dac.'//cno,6,append)
!        endif
!        if(flag_ir.ge.2) then
!           call dacdef3s(78,trim(outdir)//'/ion.dac.'//cno,6,append)
!           call dacdef3s(79,trim(outdir)//'/rec.dac.'//cno,6,append)
! !          call dacdef3s(Gm_ion,'ion.dac.',1)
! !          call dacdef3s(Gm_rec,'rec.dac.',1)
!        endif

!     endif
!     if(flag_pip.eq.1.or.flag_mhd.eq.0) then
!        do i=1,nvar_h
!           call dacdef3s(mf_h(i,1),trim(outdir) // trim(file_h(i))//cno,6,append)
!        enddo
!     endif

!   end subroutine def_varfiles

 !close file units
  subroutine epilogue
    integer i
    if(flag_mpi.eq.0 .or.my_rank.eq.0) then
       print *,"END of simulation total loop ",nt
       open(99,file=trim(outdir) // "/config.txt",status="old",form="formatted",position="append")
!       open(99,file=trim(outdir) // "/config.txt",status="replace",form="formatted")
       write(99,*)"nout:",nout
       close(99)

    endif
    end_time=MPI_Wtime()
    if(my_rank.eq.0) then
       write(*,*)"CPU TIME FOR CALCULATION IS :",end_time-start_time
    endif
!    call mk_result
    if(flag_mpi.eq.1) then
       call end_mpi
    endif
    close(mf_t)
  end subroutine epilogue

!  subroutine save_varfiles(t)
  subroutine save_varfiles(n_out)
    integer n_out
    integer i
!    double precision, intent(in) :: t


    if(n_out.ne.0) then
       call def_varfiles(1)
    endif
    write(mf_t) time
    close(mf_t)
    if(flag_mhd.eq.1) then
       do i=1,nvar_m
          call save1param(U_m(:,:,:,i),tno//trim(file_m(i)),1)
       enddo
       if(flag_resi.ge.2) then
          if(et_sav.eq.0) call save1param(eta,tno//"et.dac.",1)
       endif
       if(flag_ir.ge.1) then
!	print*,gm_ion
          if(ion_sav.eq.0) call save1param(Gm_ion,tno//'ion.dac.',1)
          if(rec_sav.eq.0) call save1param(Gm_rec,tno//'rec.dac.',1)
       endif
      if(flag_ir.eq.4) then
        !print*,Nexcite(1,1,1,:)
        call save1param(Nexcite(:,:,:,1),tno//'nexcite1.dac.',1)
        call save1param(Nexcite(:,:,:,2),tno//'nexcite2.dac.',1)
        call save1param(Nexcite(:,:,:,3),tno//'nexcite3.dac.',1)
        call save1param(Nexcite(:,:,:,4),tno//'nexcite4.dac.',1)
        call save1param(Nexcite(:,:,:,5),tno//'nexcite5.dac.',1)
        call save1param(Nexcite(:,:,:,6),tno//'nexcite6.dac.',1)
      endif
       if((flag_visc.ge.1).and.(vs_sav.eq.0)) then
          call save1param(visc(:,:,:,1),tno//"viscx.dac.",1)
          call save1param(visc(:,:,:,2),tno//"viscy.dac.",1)
          call save1param(visc(:,:,:,3),tno//"viscz.dac.",1)
       endif
    endif
    if(flag_pip.eq.1 .or.flag_mhd.eq.0) then
       do i=1,nvar_h
          call save1param(U_h(:,:,:,i),tno//trim(file_h(i)),1)
       enddo
    endif
    if(flag_divb.eq.1 .and. flag_mhd.eq.1 .and. ps_sav .eq.0) then
       call save1param(U_m(:,:,:,9),tno//trim(file_m(9)),1)
    endif

  end subroutine save_varfiles

  subroutine write_3D_array(varname, data_array)
    integer(HID_T) :: dset_id         ! Dataset identifier
    integer :: per, n
    character(*) :: varname
    double precision :: data_array(ix,jx,kx)

    ! strip off '.dac.' suffix on certain variable names
    per = index(varname, '.')
    if(per /= 0) varname = varname(1:per-1)

    ! Creating dataset
    CALL h5dcreate_f(file_id, trim(varname), H5T_NATIVE_DOUBLE, filespace_id(4), &
                     dset_id, hdf5_error)
    ! Select hyperslab in the file.
    CALL h5dget_space_f(dset_id, filespace_id(4), hdf5_error)
    CALL h5sselect_hyperslab_f(filespace_id(4), H5S_SELECT_SET_F, hdf5_offset, dimsMem, hdf5_error)
    ! write data to file
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
                    data_array(start_stop(1,1):start_stop(1,2), start_stop(2,1):start_stop(2,2), &
                               start_stop(3,1):start_stop(3,2)), &
                    dimsFile, hdf5_error, &
                    file_space_id=filespace_id(4), mem_space_id=memspace_id(4), &
                    xfer_prp=plist_id)
    ! Closing dataset connections
    CALL h5dclose_f(dset_id, hdf5_error)
  end subroutine write_3D_array

  subroutine save1param(q,name,nvar)
    integer,intent(in)::nvar
    double precision, intent(in) :: q(ix,jx,kx,nvar)
    character(*), intent(in) :: name
    integer, parameter :: mf_q = 999
    integer nn

    call dacdef3s(mf_q,trim(outdir) // '/' // name // cno,6,0)
    do nn=1,nvar
       write(mf_q) q(:,:,:,nn)
    enddo
!print*,mf_q,name,q(1,1,1,:)
    close(mf_q)
  end subroutine save1param


  subroutine set_initial_out
    integer i
    if(flag_mhd.eq.1) then
          allocate(mf_m(nvar_m,2),file_m(nvar_m))
          do i=1,nvar_m
             mf_m(i,1)=19+i
             mf_m(i,2)=i
          enddo
          file_m(1)='ro_p.dac.'
          file_m(2)='mx_p.dac.'
          file_m(3)='my_p.dac.'
          file_m(4)='mz_p.dac.'
          file_m(5)='en_p.dac.'
          file_m(6)='bx.dac.'
          file_m(7)='by.dac.'
          file_m(8)='bz.dac.'
          if(flag_divb.eq.1.or.flag_divb.eq.2) &
               file_m(9)='ps.dac.'
       endif
       if(flag_pip.eq.1.or.flag_mhd.eq.0) then
          allocate(mf_h(nvar_h,2),file_h(nvar_h))
!          allocate(file_h(nvar_h))
          do i=1,nvar_h
             mf_h(i,1)=30+i
             mf_h(i,2)=i
          enddo

          file_h(1)='ro_n.dac.'
          file_h(2)='mx_n.dac.'
          file_h(3)='my_n.dac.'
          file_h(4)='mz_n.dac.'
          file_h(5)='en_n.dac.'
       endif
  end subroutine set_initial_out
  subroutine reset_out
    if(flag_mhd.eq.1) then
       deallocate(file_m,mf_m)
    endif
    if(flag_pip.eq.1.or.flag_mhd.eq.0) then
       deallocate(file_h,mf_h)
    endif
  end subroutine reset_out




  !Modification for restart setting NN 2015/08/27
  subroutine mk_config
    character*15 key
    character*200  tmp
    integer ind_e
    open(11,file='setting.txt',form='formatted',status='old')
    open(99,file=trim(outdir) // "/config.txt",status="replace",form="formatted")
    do
       read(11,"(A)",end=999)tmp
       call get_word(tmp,key,ind_e)
       if(ind_e.gt.1) write(99,"(A)")key//":"//tmp(1:ind_e-1)
    enddo
999 continue
    close(11)
    write(99,"(A)")"ENDSETTING"
    write(99,*)flag_mhd,flag_pip, " #mhd and pip flag"
    write(99,*)nvar_h,nvar_m, " #number of variables"
    write(99,*)ix,jx,kx, " # used grid numbers"
    write(99,*)margin, " # used margin grid numbers"
    write(99,*)gm," #Abiabatic constant"
    write(99,*)flag_bnd
!    write(99,*)flag_damp,damp_time, "velocity damping"
    if(flag_mpi.eq.1) then
       write(99,*)mpi_siz," #mpi domain size"
    endif
    close(99)
  end subroutine mk_config

  !! restart routine should be modified
  subroutine restart
    integer tmp,out_tmp
    character*200 line
    character*15 key
!    open(99,file=trim(indir)//"result.txt",status="old",form="formatted"

    !Modification restart setting 2015/08/27 NN======================

    open(restart_unit,file=trim(indir)//"config.txt",status="old",form="formatted")
    do
       read(restart_unit,*,end=111)line
       if(trim(line)=="ENDSETTING") exit
    enddo
111 continue
    read(restart_unit,*)flag_mhd,flag_pip
    read(restart_unit,*)nvar_h,nvar_m
    read(restart_unit,*)ix,jx,kx
    read(restart_unit,*)margin
!    read(99,*)tmp
!    read(99,*)tmp
!    tmp=tmp+1
    read(restart_unit,*)gm
    read(restart_unit,*,end=777)flag_bnd
!    read(restart_unit,*)flag_damp,damp_time
    if(flag_mpi.eq.1) then
       read(restart_unit,*,end=777)mpi_siz
    endif
    if(flag_restart.eq.0) then
       key="nout"
       do
          read(restart_unit,"(A)",end=888)line
          call get_value_integer(line,key,out_tmp)
       enddo
888    continue
       flag_restart=out_tmp-1
    endif
!    if(flag_restart.eq.-1) then
!       key="nout"
!       do
!          read(restart_unit,"(A)",end=889)line
!          call get_value_integer(line,key,out_tmp)
!       enddo
!889    continue
!       flag_restart=out_tmp-1
!    endif
777 continue


    close(restart_unit)

    !=========================Modification end

    call reread_coordinate
    if (flag_mpi.eq.0 .or. my_rank.eq.0) then
       print *,"Now reading data from [",trim(indir),"] ..."
       print *,"start step is : ",flag_restart
    endif
    call reread_variables


    if (flag_mpi.eq.0 .or. my_rank.eq.0) print *,"reading Finish."
!    if (flag_mpi.eq.0 .or. my_rank.eq.1) print *,"dtout=",dtout

    call reconf_grid_space

    nout = flag_restart+1

    start_time=MPI_Wtime()
    tend=tend+time
!    if (flag_mpi.eq.0 .or. my_rank.eq.1) print *,"time",start_time,tend,dtout
    call initialize_IOT(dtout,tend,output_type)
  end subroutine restart

  subroutine reconf_grid_space

    dxc(2:ix-1)=0.5*(x(3:ix)-x(1:ix-2))
    dxc(1)=dxc(2) ; dxc(ix)=dxc(ix-1)
    if(ndim.ge.2) then
       dyc(2:jx-1)=0.5*(y(3:jx)-y(1:jx-2))
       dyc(1)=dyc(2) ; dyc(jx)=dyc(jx-1)
       if(ndim.ge.3) then
          dzc(2:kx-1)=0.5*(z(3:kx)-z(1:kx-2))
          dzc(1)=dzc(2) ; dzc(kx)=dzc(kx-1)
       else
          dzc(1)=1.0d2
       endif
    else
       dyc(1)=1.0d2;dzc(1)=1.0d2
    endif
  end subroutine reconf_grid_space

  subroutine reread_coordinate
    character*4 tmp_id
    write(tmp_id,"(i4.4)")mpi_pos(1)

    call dacget(51,trim(indir) // 'x.dac.'//tmp_id,ix,x)
    call dacget(51,trim(indir) // 'dx.dac.'//tmp_id,ix,dx)
    if(ndim.ge.2)then
       write(tmp_id,"(i4.4)")mpi_pos(2)
       call dacget(51,trim(indir) // 'y.dac.'//tmp_id,jx,y)
       call dacget(51,trim(indir) // 'dy.dac.'//tmp_id,jx,dy)
       if(ndim.ge.3)then
          write(tmp_id,"(i4.4)")mpi_pos(3)
          call dacget(51,trim(indir) // 'z.dac.'//tmp_id,kx,z)
          call dacget(51,trim(indir) // 'dz.dac.'//tmp_id,kx,dz)
       endif
    endif
  end subroutine reread_coordinate

  subroutine reread_variables
    integer nvar,i
    character*4 step_char
    write(step_char,"(i4.4)")flag_restart
    nvar=ix
    if(ndim.ge.2)nvar=nvar*jx
    if(ndim.ge.3)nvar=nvar*kx

    call set_initial_out
    if(flag_mhd.eq.1) then
       do i=1,nvar_m
          call dacget(mf_m(i,1),trim(indir)//step_char//trim(file_m(i))//cno,nvar, &
               U_m(:,:,:,mf_m(i,2)))
       enddo
    endif

    if(flag_pip.eq.1.or.flag_mhd.eq.0) then
       do i=1,nvar_h
          call dacget(mf_h(i,1),trim(indir)//step_char//trim(file_h(i))//cno,nvar,&
               U_h(:,:,:,mf_h(i,2)))
       enddo
    endif
!    call reset_out

    if(flag_resi.ge.1) then
       call dacget(11,trim(indir)//step_char//'et.dac.'//cno,nvar,eta(:,:,:))
    endif


    if(flag_pip.eq.1.and.flag_col.ge.1) then
       call dacget(11,trim(indir)//step_char//'ac.dac.'//cno,nvar,ac(:,:,:))
    endif

    if(flag_pip.eq.1.and.flag_ir.ge.1) then
       call dacget(11,trim(indir)//step_char//'ion.dac.'//cno,nvar,gm_ion(:,:,:))
       call dacget(11,trim(indir)//step_char//'rec.dac.'//cno,nvar,gm_rec(:,:,:))
	if (flag_ir_type .eq.0) then
		allocate(arb_heat(ix,jx,kx))
		call dacget(11,trim(indir)//step_char//'aheat.dac.'//cno,nvar,arb_heat(:,:,:))
	endif
!       print *,"GM_ION",maxval(gm_ion),minval(gm_ion)
!       print *,"GM_REC",maxval(gm_rec),minval(gm_rec)
    endif

    if(flag_grav.ge.1) then
!      do i=1,3
          call dacget(11,trim(indir)//step_char//'gr.dac.'//cno,nvar*3,gra,3)
!       enddo
    endif

    !! read time (by Tak)
    call get_time(mf_t,trim(indir)//'t.dac.0000',flag_restart,time)
    if (flag_mpi.eq.0 .or. my_rank.eq.0) print *, 'Read time: ',time
    if(my_rank.eq.0) then
       call copy_time(mf_t,trim(indir)//'t.dac.0000',trim(outdir)//'t.dac.0000',flag_restart+1)
    endif

  end subroutine reread_variables


!  subroutine dacget(idf,file,nvar,restart,var)
  subroutine dacget(idf,file,nvar,var,read_num)
    integer,intent(in)::idf
    integer,intent(in)::nvar
    character*(*),intent(in)::file
    double precision,intent(out)::var(nvar)
    integer,optional::read_num
    integer tmp,i
    integer n_read,read_size
    if(present(read_num)) then
       n_read=read_num
    else
       n_read=1
    endif
    read_size=nvar/n_read

    open(idf,file=file,form="unformatted",status="old")
    !remove .dac. header----------------------------
    do i=1,5
       read(idf)tmp
    enddo
    i=tmp
    !-----------------------------------------------
!    if(present(restart)) then
    do i=1,n_read
       read(idf)var(1+(i-1)*read_size:read_size*i)
    enddo
    close(idf)

  end subroutine dacget

  subroutine copy_time(idf,in_file,out_file,start_step)
    integer,intent(in)::idf,start_step
    integer i,tmp
    character*(*),intent(in)::in_file
    character*(*),intent(in)::out_file
    double precision time_tmp
    double precision times(start_step)
    open(idf,file=in_file,form="unformatted",status="old")
    !remove .dac. header----------------------------
    do i=1,5
       read(idf)tmp
    enddo
    !-----------------------------------------------
    do i=1,start_step
       read(idf)time_tmp
       times(i)=time_tmp
    enddo
    close(idf)

    call dacdef0s(idf,out_file,6,0)
    do i=1,start_step
       write(idf)times(i)
    end do
    close(idf)

  end subroutine copy_time

  !! Get time when the calculation is restarted (by Tak)
  subroutine get_time(idf,file,restart,time_tmp)
    integer,intent(in)::idf
    integer,intent(in)::restart
    character*(*),intent(in)::file
    double precision,intent(out)::time_tmp
    integer tmp,i
    open(idf,file=file,form="unformatted",status="old")
    !remove .dac. header----------------------------
    do i=1,5
       read(idf)tmp
    enddo
    !-----------------------------------------------
    do i=1,restart+1
       read(idf)time_tmp
    enddo
    close(idf)
  end subroutine get_time


  subroutine dacdef1d(idf,file,mtype,in)
    integer,intent(in)::idf,mtype,in
    character*(*) file
    open(idf,file=file,form='unformatted')
    write(idf)1
    write(idf)0
    write(idf)mtype
    write(idf)1
    write(idf)in
  end subroutine dacdef1d

  subroutine dacdef0s(idf,file,mtype,append)
    integer,intent(in)::idf,mtype,append
    character*(*) file
    if(append.ne.1) then
       open(idf,file=file,form='unformatted')
       write(idf)1
       write(idf)0
       write(idf)mtype
       write(idf)1
       write(idf)-1
    else
       open(idf,file=file,form='unformatted',position='append')
    endif
  end subroutine dacdef0s

  subroutine dacdef3s(idf,file,mtype,append)
    integer,intent(in)::idf,mtype,append
    character*(*) file
    if(append.ne.1) then
       open(idf,file=file,form='unformatted')
       write(idf)1
       write(idf)0
       write(idf)mtype
       write(idf)4
       write(idf)ix,jx,kx,-1
    else

       open(idf,file=file,form='unformatted',position='append')
    endif
  end subroutine dacdef3s



  subroutine stop_sim
    !----------------------------------------------------------------------|
    ! Checks on code to see if it has run for too many steps
    ! or has a too small dt
    !----------------------------------------------------------------------|
    integer i,j,k,o_count,tmp_stop,ierr
    flag_stop=0
    if (nt.eq.nmax) then
       write(6,*) 'max loop number reached'
       flag_stop=1
    endif
    if (dt .le. 1.d-9) then
       write(6,*) 't=',time,'dt=',dt
       write(6,*) 'dt too small'
       flag_stop=1
    endif
    if (flag_mhd.eq.1) then
       if(sum(U_m*0).ne.0)then
          o_count=0
          print *,"NAN Appear"
          do k=1,kx;do j=1,jx;do i=1,ix
             if(sum(U_m(i,j,k,:)*0) .ne.0.and.o_count<10) then
                print *,'NAN at (i,j,k) =  ',i,j,k
                print *,'NAN at (x,y,z) =  ',x(i),y(j),z(k)
                print *,'U_m ',U_m(i,j,k,:)
                o_count=o_count+1
             endif
          enddo;enddo;enddo
          flag_stop=1
       endif
    else
       if(sum(U_h*0).ne.0)then
          print *,"NAN Appear"
          do k=1,kx;do j=1,jx;do i=1,ix
             if(sum(U_h(i,j,k,:)*0) .ne.0.and.o_count<10) then
                print *,'NAN at (i,j,k) =  ',i,j,k
                print *,'NAN at (x,y,z) =  ',x(i),y(j),z(k)
                print *,'U_h ',U_h(i,j,k,:)
                o_count=o_count+1
             endif
          enddo;enddo;enddo
          flag_stop=1
       endif
    endif
    if(flag_mpi.eq.1) then
       call mpi_allreduce(flag_stop,tmp_stop,1,mpi_integer,MPI_MAX, &
            mpi_comm_world,ierr)
       flag_stop=tmp_stop
    endif
  end subroutine stop_sim
end module io_rot
