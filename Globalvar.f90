module globalvar
!====================================================================
! This is the global variables for the PIP code. 
! Author A.Hillier
! First version - 2013/04/23
! Second version - 2013/06/02
! Third version - 2013/06/04 NN
!====================================================================

!Define global constant------------------------------------------------------

!Coordinate system (only cartesian coordinates used)
  double precision,allocatable,save ::x(:),dx(:),y(:),dy(:),z(:),dz(:),dsc(:,:)
  double precision,allocatable,save:: dxc(:),dyc(:),dzc(:)

!Conserved variables
  double precision,allocatable,save ::U_m(:,:,:,:),U_h(:,:,:,:)

!Magnetic field at cell-face
  
!other variables
  double precision,allocatable,save ::kap(:,:,:),mu(:,:,:),eta(:,:,:)
  double precision,allocatable,save::ac(:,:,:),xi_n(:,:,:),GM_rec(:,:,:),Gm_ion(:,:,:)
  double precision,allocatable,save::gra(:,:,:,:),visc(:,:,:,:)

  double precision,allocatable,save::arb_heat(:,:,:)

!6 level collisional rates
  double precision,allocatable,save::Colrat(:,:,:,:,:),Nexcite(:,:,:,:)

!Table for the exponential integral
  double precision,allocatable,save::expinttab(:,:) !columns for input and i=0,1,2

!These are the variables for heat condution, viscosity, resistivity and gravity

! for numerical calculation
  integer,save::flag_mhd,flag_pip
  integer,save::flag_sch,flag_hll,flag_artvis,no_advection  
  integer,save::flag_resi,flag_visc,flag_heat,flag_grav,flag_amb
  integer,save::flag_b_stg,flag_divb,flag_pip_imp
  integer,save::flag_cc,flag_ex,flag_restart,flag_debug,flag_ps
  character*200,save :: flag_ini
  integer,save::flag_damp, flag_ir, flag_ir_type
  integer,save::flag_cyl
  integer,save::debug_option,debug_direction,flag_col
  double precision,save::safety,gm,beta,col,t_ir,nu_0
  !for heat conduction
  integer,save::hc_split,hc_sch,hc_max,hc_integ,hc_type,nsub_max
  double precision,save::kapper,sor_omega
  double precision, save::theta,db_clean,scl_height,b_cr
  double precision, save:: safety_cnd
  !for ambipolar diffusion
  double precision,save::pip_imp_factor,eta_0,n_fraction,damp_time
  double precision,save::vd_cri,j_cri
  double precision,save::debug_parameter
! alpha is the collision parameter and flag_sch is the flag for the scheme you are using
  double precision,save :: cmax

! Radiative losses
  integer,save :: flag_rad
  double precision,save :: rad_temp
  double precision,allocatable,save::GM_rec_rad(:,:,:),Gm_ion_rad(:,:,:),radrat(:,:,:,:,:)

!for coordinate system
  integer,save::ix,jx,kx,nvar_h,nvar_m,margin(3),ndim,ig(3)
! marginz is set to 0 for 2d to let all the 3d codes that use margin work

!timestepping variables
  double precision,save ::time,dt,dtout,dt_cnd,start_t
  integer,save:: flag_stop,flag_time,t_order,s_order,total_iter
  integer,save:: flag_hc_test,nmax=500000000
  
!for output
  double precision,save:: tend ! time at which the data should be output and end time of simulation
  character*200,save:: outdir,indir
  integer, save :: nout,nt,output_type ! number of times data has been output
!  integer,allocatable,save::mf_m(:,:),mf_h(:,:)

! For the boundary conditions
  integer,save::flag_bnd(6)


!for mpi 
  integer,save::flag_mpi,flag_mpi_split,neighbor(6),my_rank,total_prc 
  integer,save::mpi_siz(3),mpi_pos(3),vmpi
  character*4 cno  

!for limit
  double precision,save :: ro_lim,pr_lim
  double precision,parameter :: tiny=1d-20

!For damping across time MAKE ALLOCATABLE??, CHOOSE SENSIBLE VALUE?
  double precision,save :: oldke_damp=1.0d0
  
!Emergency save procedure
integer,save :: esav
double precision,save :: emsavtime

!Normalisation values
  double precision, save :: T0, n0, L0

!Other useful values
  double precision, save :: f_p_ini,f_p_p_ini,n0fac,Gm_rec_ref !Constants for ionisation and recombination

!optional save parameters
  integer, save :: ac_sav, xi_sav, ion_sav, rec_sav, col_sav, gr_sav, vs_sav, heat_sav, et_sav, ps_sav
end module globalvar
