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
  subroutine output(out)
    integer,intent(in)::out
    integer i,j,k,outesav
    double precision total_divB,cx,cy,max_C,divb

    ! if nout = 0 initial setting for output is done
    nt=nt+1
    if(nout.eq.0) then
       call set_initial_out
       start_time=MPI_Wtime()
       if(flag_mpi.eq.0 .or. my_rank.eq.0) then
          call mk_config
       endif
       call initialize_IOT(dtout,tend,output_type)
    endif

    !Check for physical-time save
    end_time=MPI_Wtime()
    outesav=-1
    if (((end_time-start_time) .ge. emsavtime) .and. (esav .eq. 1)) then
      outesav=1
      esav=2
      if(flag_mpi.eq.0 .or.my_rank.eq.0) &
        write(6,*) 'calling emergency save'
    endif

    ! output variables
    if((time.ge.get_next_output(nout,time,esav)) .or. (out.eq.1) .or. (outesav.eq.1)) then
      end_time=MPI_Wtime()

      if(flag_mpi.eq.0 .or.my_rank.eq.0) &
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
          max_C=max(max_C,((-U_m(i+2,j,1,7)+8.0d0*U_m(i+1,j,1,7) &
                -8.0d0*U_m(i-1,j,1,7)+U_m(i-2,j,1,7))/(12.0d0*dx(i)))**2+ &
                ((-U_m(i,j+2,1,6)+8.0d0*U_m(i,j+1,1,6) &
                -8.0d0*U_m(i,j-1,1,6)+U_m(i,j-2,1,6))/(12.0d0*dy(j)))**2)
        enddo;enddo
        print *,"NT,TOTAL_DIVB, maxJ =",nt,total_divb,sqrt(max_C)
      endif

      ! create output parallel HDF5 file for all data in time-step
      call create_hdf5_outfile
      ! create file & local-memory dataspaces
      call create_write_dataspaces
      ! save spatial grid coordinates
      call save_coordinates
      ! save initial variables to output HDF5 file
      if(nout .eq. 0) call def_varfiles(0)
      ! Save simulation variables to output HDF5 file
      call save_varfiles(nout)
      call close_write_outfile

      nout=nout+1
    endif
  end subroutine output

  subroutine create_hdf5_outfile
    character(40) :: file_path

    ! Initialize FORTRAN predefined datatypes
    CALL h5open_f(hdf5_error)

    ! Setup file access property list with parallel I/O access.
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf5_error)
    CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, hdf5_error)

    ! Create the output file
    write(tno, "(i4.4)") nout
    file_path = trim(outdir) // 't' // tno // '.h5'
    CALL h5fcreate_f(trim(file_path), H5F_ACC_TRUNC_F, file_id, hdf5_error, access_prp = plist_id)

    ! Create property list for collective dataset write
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf5_error)
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf5_error)
  end subroutine create_hdf5_outfile

  subroutine create_write_dataspaces
    integer(HSIZE_T) :: dims_1D(1)
    integer(HSIZE_T) :: unif_cell_size
    integer :: i

    ! iterate over grid coordinates (X, Y, Z)
    do i=1,3
      unif_cell_size = (dimsFile(i)-2*margin(i))/mpi_siz(i)
      ! Define local process start index & global offset
      if(neighbor(2*i-1).eq.-1 .or. mpi_pos(i).eq.0) then
        start_stop(i,1) = 1
        hdf5_offset(i) = 0
      else
        start_stop(i,1) = 1 + margin(i)
        hdf5_offset(i) = margin(i) + unif_cell_size*mpi_pos(i)
      end if
      ! Define local process stop index
      if(neighbor(2*i).eq.-1 .or. mpi_pos(i).eq.(mpi_siz(i)-1)) then
        start_stop(i,2) = ig(i)
      else
        start_stop(i,2) = ig(i) - margin(i)
      end if

      dimsMem(i) = start_stop(i,2) - start_stop(i,1) + 1
      ! Create 1D dataspaces for spatial coordinates
      dims_1D(1) = dimsMem(i)
      CALL h5screate_simple_f(1, dims_1D, memspace_id(i), hdf5_error)
      dims_1D(1) = dimsFile(i)
      CALL h5screate_simple_f(1, dims_1D, filespace_id(i), hdf5_error)
    end do

    ! Create 3D dataspaces
    CALL h5screate_simple_f(3, dimsFile, filespace_id(4), hdf5_error)
    CALL h5screate_simple_f(3, dimsMem, memspace_id(4), hdf5_error)
  end subroutine create_write_dataspaces

  subroutine close_write_outfile
    integer :: i

    ! Close the dataspaces
    do i=1,4
      CALL h5sclose_f(filespace_id(i), hdf5_error)
      CALL h5sclose_f(memspace_id(i), hdf5_error)
    end do
    ! Close the file and property list.
    CALL h5pclose_f(plist_id, hdf5_error)
    CALL h5fclose_f(file_id, hdf5_error)
  end subroutine close_write_outfile

  subroutine save_coordinates
    ! Save the coordinate grids into the parallel HDF5 file
    call write_1D_array(1, "xgrid", x)
    call write_1D_array(2, "ygrid", y)
    call write_1D_array(3, "zgrid", z)
    ! Also include the grid spacings
    call write_1D_array(1, "dx", dx)
    call write_1D_array(2, "dy", dy)
    call write_1D_array(3, "dz", dz)
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

    if(flag_pip.eq.1.or.flag_amb.eq.1) then
      if(ac_sav.eq.0) call write_3D_array("ac", ac)
      if(xi_sav.eq.0) call write_3D_array("xi", xi_n)
    endif
    if(flag_pip.eq.1.and.flag_ir.eq.1) then
      if(ion_sav.eq.0) call write_3D_array("ion", Gm_ion)
      if(rec_sav.eq.0) call write_3D_array("rec", Gm_rec)
    endif
    if(flag_mhd.eq.1.and.flag_resi.eq.1) then
      if(et_sav.eq.0) call write_3D_array("et", eta)
    endif
    if(flag_col.eq.1) then
      if(col_sav.eq.0) call write_3D_array("col", ac)
    endif

    if(flag_grav.eq.1) then
      if(gr_sav.eq.0) call write_3D_array("gr", gra)
    endif
    if(flag_visc.eq.1) then
      if(vs_sav.eq.0) call write_3D_array("vs", mu)
    endif
    if(flag_pip.eq.1.and.flag_ir_type.eq.0.and.flag_IR.ne.0) then
      if(heat_sav.eq.0) call write_3D_array("aheat", arb_heat)
    endif

    if(flag_mpi.eq.0 .or.my_rank.eq.0) then
      call dacdef0s(mf_t,trim(outdir) // 't.dac.'//cno,6,append)
    end if
  end subroutine def_varfiles

  !close file units
  subroutine epilogue
    integer i
    if(flag_mpi.eq.0 .or.my_rank.eq.0) then
      print *,"END of simulation total loop ",nt
      open(99,file=trim(outdir) // "/config.txt",status="old",form="formatted",position="append")
      write(99,*)"nout:",nout
      close(99)
    endif
    end_time=MPI_Wtime()
    if(my_rank.eq.0) then
      write(*,*)"CPU TIME FOR CALCULATION IS :",end_time-start_time
    endif

    if(flag_mpi.eq.1) then
      call end_mpi
      ! Close FORTRAN interface
      if(my_rank.eq.0) call h5close_f(hdf5_error)
    endif
    close(mf_t)
  end subroutine epilogue

  subroutine save_varfiles(n_out)
    integer n_out, i
    character(8) :: Nexc_name

    if(n_out.ne.0) call def_varfiles(1)
    write(mf_t) time
    close(mf_t)

    if(flag_mhd.eq.1) then
      do i=1,nvar_m
        if((i.lt.9) .or. (flag_divb.eq.1 .and. ps_sav .eq.0)) then
          call write_3D_array(trim(file_m(i)), U_m(:,:,:,i))
        end if
      enddo
      if(flag_resi.ge.2) then
        if(et_sav.eq.0) call write_3D_array("et", eta)
      endif
      if(flag_ir.ge.1) then
        if(ion_sav.eq.0) call write_3D_array("ion", Gm_ion)
        if(rec_sav.eq.0) call write_3D_array("rec", Gm_rec)
      endif
      if(flag_ir.eq.4) then
        do i=1,6
          write(Nexc_name, '(a, i0)') 'nexcite', i
          call write_3D_array(Nexc_name, Nexcite(:,:,:,i))
        end do
      endif
      if((flag_visc.ge.1).and.(vs_sav.eq.0)) then
        call write_3D_array("viscx", visc(:,:,:,1))
        call write_3D_array("viscy", visc(:,:,:,2))
        call write_3D_array("viscz", visc(:,:,:,3))
      endif
    endif
    if(flag_pip.eq.1 .or.flag_mhd.eq.0) then
      do i=1,nvar_h
        call write_3D_array(trim(file_h(i)), U_h(:,:,:,i))
      enddo
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

    !Modification restart setting 2015/08/27 NN======================
    open(restart_unit,file=trim(indir)//"config.txt",status="old",form="formatted")
    do
      read(restart_unit,*,end=111)line
      if(trim(line)=="ENDSETTING") exit
    enddo
111 continue
    read(restart_unit,*)flag_mhd,flag_pip
    read(restart_unit,*)nvar_h,nvar_m
    read(restart_unit,*)
    read(restart_unit,*)margin
    read(restart_unit,*)gm
    read(restart_unit,*,end=777)flag_bnd
    if(flag_mpi.eq.1) then
      read(restart_unit,*,end=777)
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
777 continue
    close(restart_unit)

    !=========================Modification end
    if (flag_mpi.eq.0 .or. my_rank.eq.0) then
       print *,"Now reading data from [",trim(indir),"] ..."
       print *,"start step is : ",flag_restart
    endif
    ! open file and initialize local-memory dataspaces
    call read_restart_hdf5
    ! load spatial grid coordinates
    call reread_coordinate
    ! load simulation variables
    call reread_variables
    call close_restart_hdf5
    if (flag_mpi.eq.0 .or. my_rank.eq.0) print *,"reading Finish."

    call reconf_grid_space

    nout = flag_restart+1

    start_time=MPI_Wtime()
    tend=tend+time
    call initialize_IOT(dtout,tend,output_type)
  end subroutine restart

  subroutine read_restart_hdf5()
    character(40) :: file_path
    character*4 step_char
    integer(HSIZE_T) :: dims_1D(1)
    integer(HSIZE_T) :: proc_block
    integer :: i

    CALL h5open_f(hdf5_error)
    ! open HDF5 file for requested time-step
    write(step_char, "(i4.4)") flag_restart
    file_path = trim(indir) // 't' // step_char // '.h5'
    CALL h5fopen_f(trim(file_path), H5F_ACC_RDONLY_F, file_id, hdf5_error)

    ! iterate over grid coordinates (X, Y, Z)
    do i=1,3
      dimsMem(i) = ig(i)

      proc_block = (dimsFile(i)-2*margin(i))/mpi_siz(i)
      hdf5_offset(i) = proc_block*mpi_pos(i)
      ! Set 1D memory dataspace
      dims_1D(1) = dimsMem(i)
      CALL h5screate_simple_f(1, dims_1D, memspace_id(i), hdf5_error)
    end do
    ! Set 3D memory dataspace
    CALL h5screate_simple_f(3, dimsMem, memspace_id(4), hdf5_error)
  end subroutine read_restart_hdf5

  subroutine close_restart_hdf5()
    integer :: i

    ! Close the memory dataspaces
    do i=1,4
      CALL h5sclose_f(memspace_id(i), hdf5_error)
    end do
    ! Close the file.
    CALL h5fclose_f(file_id, hdf5_error)
  end subroutine close_restart_hdf5

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
    ! Save the coordinate grids into the parallel HDF5 file
    call read_1D_array(1, "xgrid", x)
    call read_1D_array(2, "ygrid", y)
    call read_1D_array(3, "zgrid", z)
    ! Also include the grid spacings
    call read_1D_array(1, "dx", dx)
    call read_1D_array(2, "dy", dy)
    call read_1D_array(3, "dz", dz)
  end subroutine reread_coordinate

  subroutine read_1D_array(n, varname, data_out)
    integer(HID_T) :: dset_id, dataspace_id
    integer(HSIZE_T) :: dims_1D(1)
    integer(HSSIZE_T) :: offsetFile(1)
    integer :: n
    character(*) :: varname
    double precision, dimension(1) :: data_out(dimsMem(n))

    CALL h5dopen_f(file_id, varname, dset_id, hdf5_error)
    ! Get file dataspace
    dims_1D(1) = dimsMem(n)
    offsetFile(1) = hdf5_offset(n)
    CALL h5dget_space_f(dset_id, dataspace_id, hdf5_error)
    CALL h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, offsetFile, dims_1D, hdf5_error)
    ! read file
    dims_1D(1) = dimsFile(n)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, dims_1D, hdf5_error, &
                   mem_space_id=memspace_id(n), file_space_id=dataspace_id)
    ! close connections
    CALL h5sclose_f(dataspace_id, hdf5_error)
    CALL h5dclose_f(dset_id, hdf5_error)
  end subroutine read_1D_array

  subroutine reread_variables
    integer i

    call set_initial_out

    if(flag_mhd.eq.1) then
      do i=1,nvar_m
        call read_3D_array(trim(file_m(i)), U_m(:,:,:,mf_m(i,2)))
      enddo
    endif
    if(flag_pip.eq.1.or.flag_mhd.eq.0) then
      do i=1,nvar_h
        call read_3D_array(trim(file_h(i)), U_h(:,:,:,mf_h(i,2)))
      enddo
    endif

    if(flag_resi.ge.1) then
      call read_3D_array('et', eta(:,:,:))
    endif

    if(flag_pip.eq.1.and.flag_col.ge.1) then
      call read_3D_array('ac', ac(:,:,:))
    endif
    if(flag_pip.eq.1.and.flag_ir.ge.1) then
      call read_3D_array('ion', gm_ion(:,:,:))
      call read_3D_array('rec', gm_rec(:,:,:))
      if (flag_ir_type .eq.0) then
        allocate(arb_heat(ix,jx,kx))
        call read_3D_array('aheat', arb_heat(:,:,:))
  	  endif
    endif

    if(flag_grav.ge.1) then
      call read_3D_array('gr', gra)
    endif

    !! read time (by Tak)
    call get_time(mf_t,trim(indir)//'t.dac.0000',flag_restart,time)
    if (flag_mpi.eq.0 .or. my_rank.eq.0) print *, 'Read time: ',time
    if(my_rank.eq.0) then
       call copy_time(mf_t,trim(indir)//'t.dac.0000',trim(outdir)//'t.dac.0000',flag_restart+1)
    endif
  end subroutine reread_variables

  subroutine read_3D_array(varname, data_out)
    integer(HID_T) :: dset_id, dataspace_id     ! dataset and file-dataspace IDs
    integer(HSIZE_T) :: dimsFile(3), dimsMem(3)
    integer(HSSIZE_T) :: offsetFile(3)
    integer :: per
    character(*) :: varname
    double precision, dimension(3) :: data_out(ix,jx,kx)

    ! strip off '.dac.' suffix on certain variable names
    per = index(varname, '.')
    if(per /= 0) varname = varname(1:per-1)

    CALL h5dopen_f(file_id, varname, dset_id, hdf5_error)
    ! Get file dataspace
    CALL h5dget_space_f(dset_id, dataspace_id, hdf5_error)
    CALL h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, hdf5_offset, dimsMem, hdf5_error)
    ! read file
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, dimsFile, hdf5_error, &
                   mem_space_id=memspace_id(4), file_space_id=dataspace_id)
    ! close connections
    CALL h5sclose_f(dataspace_id, hdf5_error)
    CALL h5dclose_f(dset_id, hdf5_error)
  end subroutine read_3D_array

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
