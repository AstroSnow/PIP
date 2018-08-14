subroutine set_coordinate_NUG(is_all,is,s,ds,dsc,s_order, &
     s0,s1,ds0,s_fact,ds_max,i_rank,margin,sym_flag,flag_err)
  !set non-unifrom-grid with mpi
  implicit none
  integer,intent(in)::is,is_all,i_rank,margin,sym_flag,s_order
  integer is0,i,i_top,i_total,l_flag,u_flag,i_start
  double precision,intent(out)::s(is),ds(is),dsc(is)
  double precision,intent(in):: s0,s1,ds0,s_fact,ds_max
  integer,intent(out):: flag_err
  double precision s_all(is_all),ds_all(is_all)

  flag_err=0
  i_total=(is_all-2*margin)/(is-2*margin)

  !! initialization
  s(:) = 0.d0 ; ds(:) = 0.d0 ; dsc(:) = 0.d0
  s_all(:) = 0.d0 ; ds_all(:) = 0.d0
  
  if(sym_flag.eq.0) then
     is0=margin+1
  else
     is0=is_all/2+1
  endif
  s_all(is0)=ds0*0.5d0
  ds_all(:)=ds0    


  do i=is0,is_all-margin-1
     s_all(i+1)=s_all(i)+ds_all(i)
     if(s_all(i+1).ge.s1) then
        if (ds_all(i).le.ds_max) then
           ds_all(i+1)=ds_all(i)*s_fact
        else
           ds_all(i+1)=ds_all(i)
        endif
     endif
  enddo

  
  !! Return error when the uniform grid region is smaller than s1
  if(maxval(s_all).lt.s1) then
     flag_err = 1
  endif
  
  if(i_rank.eq.0) then
     print*,'Domain size: ',minval(s_all),maxval(s_all)
  endif  

  i_top=is_all-margin
  ds_all(i_top)=ds_all(i_top-1)
  do i=1,margin        
     if(i.ne.1) then
        ds_all(i_top+i-1)=ds_all(i_top-i+1)
     endif
     s_all(i_top+i)=s_all(i_top+i-1)+ds_all(i_top+i-1)
  enddo
  
  if(sym_flag.eq.0) then
     ds_all(margin)=ds0
     do i=1,margin                
        if(i.ne.1) then
           ds_all(margin+1-i)=ds_all(margin+i-1)
        endif
        s_all(margin+1-i)=s_all(margin+2-i)-ds_all(margin+1-i)        
     enddo         
  else 
     do i=1,is0
        s_all(i)=-s_all(is_all-i+1)        
        ds_all(i)=ds_all(is_all-i)
     enddo     
  endif

  i_start=i_rank*(is-2*margin)
  do i=1,is
     s(i)=s_all(i+i_start)
     ds(i)=ds_all(i+i_start)
  enddo
  if(s_order.eq.4) then
     if(i_rank.eq.0) then
        l_flag=0
     else
        l_flag=2
     endif
     if(i_rank.eq.i_total-1) then
        u_flag=0
     else
        u_flag=2
     endif        
     do i=3-l_flag,is-2+u_flag
        dsc(i)=0.25d0*(ds_all(i+i_start+1)+ds_all(i+i_start) &
             +ds_all(i+i_start-1)+ds_all(i+i_start-2))
     enddo
     if(l_flag.eq.0) dsc(1:2)=dsc(3)
     if(u_flag.eq.0) dsc(is-1:is)=dsc(is-2)    
  else
     if(i_rank.eq.0) then
        l_flag=0
     else
        l_flag=1
     endif
     if(i_rank.eq.i_total-1) then
        u_flag=0
     else
        u_flag=1
     endif        
     do i=2-l_flag,is-1+u_flag
!        dsc(i)=0.5d0*(ds_all(i+i_rank*(is-margin))+ds_all(i+i_rank*(is-margin)-1))
        dsc(i)=0.5d0*(ds_all(i+i_start)+ds_all(i+i_start-1))
     enddo
     
     if(l_flag.eq.2) dsc(1)=dsc(2)
     if(u_flag.eq.2) dsc(is)=dsc(is-1)             
  endif  

end subroutine set_coordinate_NUG
