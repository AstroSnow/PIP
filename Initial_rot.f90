module initial_rot 
!====================================================================
! This is the initial module for the PIP code.
! Module for a personal library of your initial conditions
! First version - 2013/04/26 NN  
! Second version - 2013/07/17 NN Move shocktube setting to shock_tube.f90  
!====================================================================
  use model_rot,only:get_parameters,allocate_vars,set_dsc
  use globalvar,only:ix,jx,kx,time,dt,flag_stop,&
       flag_mpi,cno,neighbor,nt,x,y,z,dx,dy,dz,gra,&
       flag_restart,flag_ini,dsc,dxc,dyc,dzc,nout,U_m,U_h,tend,my_rank,start_t
  use mpi_rot,only:set_mpi,set_mpi_neighbor,my_mpi_barrier
  use IO_rot,only:restart,output
  use solver_rot,only:set_coefficients
  use boundary_rot,only:initialize_bnd,PIPbnd
  use Scheme_rot,only:cfl  
  implicit none
contains
  !============================================================================
  !start main part of initial setting
  !============================================================================
  subroutine prologue    
    integer i
    time=0.d0
    nt=0
    dt=0.d0
    flag_stop=0
    !get parameters from setting file
    call get_parameters
    if(flag_mpi.eq.1)then 
       call set_mpi
    else
       write(cno,"(i4.4)")0
       neighbor(:)=-1
    endif    
    !allocate variables array
    call allocate_vars

    !set initial coordinates and physical variables
    if(flag_restart.ge.0) then
       call restart       
       call set_dsc
       start_t=time
       if (my_rank.eq.0) print *, start_t,time
    else         
       select case(flag_ini)
       case('linear_wave')          
          call linear_wave
       case('shock_tube') 
          call shock_tube
       case('explosion')          
          call explosion
       case('currentsheet')          
          call currentsheet
       case('RT')          
          call RT
       case('KH')          
          call KH
       case('Orszag_Tang')
          call Orszag_Tang
       case('FieldLoop')          
          call FieldLoop
       case('ambipolar')          
          call ambipolar
       case('PNcoupling')          
          call PNcoupling
       case('HCtest')          
          call HCtest
       case('HCtest_Tonly')
          call HCtest_Tonly
       case('Flare')          
          call Flare
       case('Coalescence')          
          call Coalescence
       case('Sedov')          
          call Sedov
       case('tearing')          
          call tearing
       case('NUwave')          
          call NUwave
       case('Reconnection')          
          call Reconnection
       case('BMP_density')          
          call BMP_density
       case('Alfven_damping')          
          call Alfven_damping
       case('blob')          
          call blob
       case('field_diffusion')          
          call field_diffusion
       case('CS_collapse')          
          call CS_collapse
       case('asym_currentsheet')          
          call asym_currentsheet
       case('CSC')          
          call CSC
       case('cnd_tube') 
          call cnd_tube
       case('MRI') 
          call MRI
       case('disk_flare') 
          call disk_flare
       case('load_prom')
          call mass_load_prom
       case('relax_prom')
          call relax_prom
       case('relax_prom2')
          call relax_prom2
       case('hsstatic')
          call hsstatic
       case('resonator')
          call resonator
       case('ionrectest')
          call ionrectest
       case('shock_tube_ion')
          call shock_tube_ion
       case('shock_tube_stab')
          call shock_tube_stab
       case('shock_tube_stab2')
          call shock_tube_stab2
       case('com_spec')
          call Complete_spectrum
       end select
       call set_coefficients(U_h,U_m,0)              
       start_t=0.d0
    endif

    !initialize boundary
    call initialize_bnd
    if(flag_mpi.eq.1) call set_mpi_neighbor
    dsc(1:ix,1)=dxc
    dsc(1:jx,2)=dyc
    dsc(1:kx,3)=dzc
    call my_mpi_barrier
    !set dt    
    call cfl(U_h,U_m)    
    if(flag_restart.lt.0) then
       nout=0
       !first output 
       call output(0)
    endif
!    start_time=0.d0
    !Boundary conditions to IC
    call PIPbnd(U_h,U_m,0)
!    if(flag_grav.eq.1) call bnd_grav
    
    
  end subroutine prologue
end module initial_rot
