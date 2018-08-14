program main
  !====================================================================
  ! This is the main program for the PIP2 code. (2014/07/08)
  ! Author A.Hillier
  ! First version - 2013/04/23
  ! Second version - 2013/05/29 AH - Main changes to put the do loop on the top level
  ! Third version -2013/06/04 NN bug fix
  ! Fourth version - 2013/06/28 AH Changed dinmension switch to inside solver routines
  ! Fifth version - 2013/07/10 Moved scheme switch into solver routines
  ! PIP2-----------------------
  ! First version - 2014/07/08 NN
  !====================================================================
  !  use globalvar,only:U_h,U_m,nt,time,dt,flag_stop,tend,nmax
  use globalvar,only:U_h,U_m,nt,time,dt,flag_stop,tend,nmax,nout,my_rank  
  use initial_rot,only:prologue
  use io_rot,only:output,stop_sim,epilogue
  use solver_rot,only:run_solver
  use scheme_rot,only:cfl
  implicit none
  !----------------------------------------------------------------------|
  ! initial setting
  !----------------------------------------------------------------------|
  call prologue
  !----------------------------------------------------------------------|
  ! Time intergration
  !----------------------------------------------------------------------|
  do while(time .lt. tend)
     ! calculate cfl condition  -----------------------------------------
     call cfl(U_h,U_m)          
     ! Time intergration  -----------------------------------------------
     call run_solver(U_h,U_m) 
     ! Output of simulation data  ---------------------------------------
     if (dt.ne.dt) print *, my_rank,dt
     time=time+dt !! NB: dt can be modified by STS routines
     call output(0)
     ! Error check  -----------------------------------------------------
     call stop_sim
     if (flag_stop .eq.1) then
        call output(1)
        exit
     endif
  enddo
  !----------------------------------------------------------------------|
  ! End of simulation
  !----------------------------------------------------------------------|
  print *, my_rank, tend,time
  call epilogue  
end program main
