module IOT_rot
    use globalvar,only:start_t,flag_restart,my_rank
  implicit none
  double precision,save:: dt,te,tnext
  double precision,parameter::factor=1.0d1
  integer,save:: type
contains
  subroutine initialize_IOT(dtout,tend,output_type)
    integer,intent(in)::output_type
    double precision,intent(in)::dtout,tend
    dt=dtout
    te=tend
    type=output_type
    tnext=0.0d0
  end subroutine initialize_IOT


  function get_next_output(nout,time,esav)
    integer,intent(in)::nout,esav    
    double precision,intent(in)::time
    double precision get_next_output
    if(type.eq.0) then
       if (flag_restart.ge.0) then
         get_next_output=start_t+(nout-flag_restart)*dt
!       if (my_rank.eq.0) print *,get_next_output,start_t
       else
         get_next_output=nout*dt
	if (esav .eq. 2) get_next_output=(nout-1)*dt
       endif
    else if(type.eq.1) then
       if(nout.eq.0) then
          get_next_output=0.0d0
          tnext=dt
       else
          get_next_output=tnext
          if(time.ge.tnext) tnext=tnext*factor
       endif
    else
       get_next_output=get_free_output(nout,time)
    endif
  end function get_next_output
  
  function get_free_output(nout,time)
    integer,intent(in)::nout    
    double precision,intent(in)::time
    double precision get_free_output
    if(nout.eq.0) then
       get_free_output=0.0d0
       tnext=dt
    else
       get_free_output=tnext
       if(time.ge.tnext)  then
          if(time.lt.0.55d0) then
             tnext=tnext+0.05d0
          else if(time.le.3.0d0) then
             tnext=0.5d0+(tnext-0.5d0)*2.0d0
          endif
       endif
    endif
  end function get_free_output
end module IOT_rot
