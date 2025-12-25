module rt_parameters
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    use kind_module
    implicit none
!   physical parameters 
    real(wp) :: freq
    !! Freq,     RF frequency, GHz
    real(wp) :: xmi1
    !!  Mi1/Mp,  relative mass of ions 1
    real(wp) :: zi1
    !! charge of ions 1
    real(wp) :: xmi2
    !! Mi2/Mp,  relative mass of ions 2
    real(wp) :: zi2
    !! charge of ions 2
    real(wp) :: dni2 
    !!  0.03   Ni2/Ni1, relative density of ions 2
    real(wp) :: xmi3
    !!  Mi3/Mp,  relative mass of ions 3
    real(wp) :: zi3
    !!  charge of ions 3
    real(wp) :: dni3
    !!  Ni3/Ni1, relative density of ions 3

!!!!!!!!!!!!!  parameters for alphas calculations !!!
    integer  ::  itend0
    !! itend0,   if = 0, no alphas
    real(wp) ::  energy
    !! energy,   max. perp. energy of alphas (MeV)
    real(wp) ::  factor   
    !! factor,   factor in alpha source
    real(wp) ::  dra   
    !! dra,      relative alpha source broadening (dr/a)
    integer  ::  kv    
    !! kv,       V_perp  greed number    

!!!!!!!!!!!!! numerical parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer  ::  nr     
    !! nr,  radial grid number  <= 505 (number of mag surface for FRTC)
    real(wp) ::  hmin1
    !! hmin1, rel.(hr) min. step in the Fast comp. mode, <1.d0
    real(wp) ::  rrange
    !! rrange,   rel.(hr) size of a 'turning' point region, <1.d0
    real(wp) ::  eps, eps_const
    !! eps,      accuracy
    real(wp) ::  hdrob
    !! hdrob,    h4 correction,
    real(wp) ::  cleft
    !! cleft,    left Vz plato border shift (<1)
    real(wp) ::  cright
    !! cright,   right Vz plato border shift (>1)
    real(wp) ::  cdel
    !! cdel,     (left part)/(Vz plato size)
    real(wp) ::  rbord 
    !! rbord,    relative radius of reflection, <1.
    real(wp) ::  pchm
    !! pchm,     threshold between 'strong' and weak' absorption, <1.
    real(wp) ::  pabs0
    !! pabs,     part of remaining power interp. as absorption
    real(wp) ::  pgiter
    !! pgiter,   relative accuracy to stop iterations
    integer ::   ni1     
    !! ni1,      grid number in the left part of Vz plato
    integer ::   ni2     
    !! ni2,      grid number in the right part of Vz plato
    integer ::   niterat     
    !! niterat,  maximal number of iterations
    integer ::   nmaxm(4)
    !! nmaxm(1), permitted reflections at 0 iteration
    !! nmaxm(2), permitted reflections at 1 iteration
    !! nmaxm(3), permitted reflections at 2 iteration
    !! nmaxm(4), permitted reflections at 3 iteration
    integer ::   maxstep2  
    !! maxstep2, maximal steps' number in Fast comp. mode
    integer ::   maxstep4    
    !! maxstep4, maximal steps' number in Slow comp. mode    

!!!!!!!!!!!!!  options !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer ::  ipri
    !! ipri, printing output monitoring: 0,1,2,3,4
    integer ::  iw, iw_start
    !! iw, initial mode (slow=1, fast=-1)
    integer ::  ismth
    !! ismth, if=0, no smoothing in Ne(rho),Te(rho),Ti(rho)
    integer ::  ismthalf
    !! ismthalf,  if=0, no smoothing in D_alpha(vperp)
    integer ::  ismthout
    !! ismthout,  if=0, no smoothing in output profiles
    integer ::  inew
    !! inew=0 for usual tokamak&Ntor_grill; 1 or 2 for g' in ST&Npol_grill
    integer ::  itor
    !! itor,      +-1, Btor direction in right coord{drho,dteta,dfi}
    integer ::  i_pol
    !! ipol,      +-1, Bpol direction in right coord{drho,dteta,dfi}

  !!!!!!!!!!!!!  grill parameters and input LH spectrum !!!!!!!!!!!!
    real(wp) ::  zplus 
    !! Zplus,    upper grill corner in centimeters
    real(wp) ::  zminus
    !! Zminus,   lower grill corner in centimeters
    real(wp) ::  ZGapPlus
    !! ZGapPlus -  upper grill gap corner in centimeters
    real(wp) ::  ZGapMinus
    !! ZGapMinus - lower grill gap corner in centimeters
    integer ::   ntet
    !! ntet,     theta grid number
    integer ::   nnz
    !! nnz,      N_phi grid number    

    integer ::   spectrum_type
    !! spectrum type 1 - 1D, 2 = 2D, 3, scatter
    logical ::   spectrum_PWM 
    !! PWM approximaton on/off
    integer ::   spectrum_coordinate_system
    !! coordinate system of spectrum
    logical  :: upl_fix
    !! floag for fixing upl value
    real(wp) :: upl_value
    !! used upl value 
    integer  :: fp_solver
    !! fp solfer 
    !! 0 default savelyev solver
    !! 1 next solver
    integer :: traj_len_seved
    !! traj length seved
    !! - 1 save all
    !! 0 don't save
    !! >1 length seved
    integer :: max_number_of_traj
    !! maximum number of trajectories
    integer :: max_size_of_traj
    !! maximum length of trajectories

    contains      
    subroutine show_parameters()          
      print*, "Freq = ", freq          
      print*, "xmi1 = ", xmi1     
      print*, "zi1 = ",  zi1     
      print*, "xmi2 = ", xmi2     
      print*, "zi2 = ",  zi2     
      print*, "dni2 = ", dni2     

      print*, "------- grill parameters and input LH spectrum "
      print*, "zplus = ",  zplus     
      print*, "zminus = ", zminus     
      print*, "ntet = ",  ntet     
      print*, "nnz = ", nnz      
      
      print*, "---------- options --------------"
      print*, "upl_fix = ", upl_fix
      print*, "upl_value = ", upl_value 
      print*, "fp_solver = ", fp_solver 
      print*, "max_number_of_traj = ", max_number_of_traj 
      print*, "max_size_of_traj = ", max_size_of_traj       

      print*, "---------- spectrum --------------"
      print*, "spectrum_type = ", spectrum_type
      print*, "spectrum_PWM = ", spectrum_PWM 
      print*, "spectrum_coordinate_system = ", spectrum_coordinate_system       

      
    end subroutine show_parameters
    
    subroutine read_nml_parameters(file_path)
        !! Read some parmeters,  Here we use a namelist 
        !! but if you were to change the storage format (TOML,or home-made), 
        !! this signature would not change

        character(len=*),  intent(in)  :: file_path
        !integer, intent(out) :: type_
        integer                        :: file_unit, iostat

        real(wp) ::  pabs
        !! pabs,     part of remaining power interp. as absorption
        integer ::  ipol
        integer ::   nmaxm_1, nmaxm_2, nmaxm_3, nmaxm_4

        ! Namelist definition===============================
        namelist /physical_parameters/ freq, xmi1, zi1, xmi2, zi2, dni2, xmi3, zi3, dni3
        namelist /parameters_for_alphas_calculations/ itend0, energy, factor, dra, kv
        namelist /numerical_parameters/ &
        nr, hmin1, rrange, eps, hdrob, cleft, &
        cright, cdel, rbord, pchm, pabs, pgiter, ni1, ni2, &
        niterat, nmaxm_1, nmaxm_2, nmaxm_3, nmaxm_4, &
        maxstep2, maxstep4
        namelist /options/ ipri, iw, ismth, ismthalf, ismthout, inew, itor, ipol, &
            upl_fix, upl_value, &
            fp_solver, traj_len_seved, &
            max_number_of_traj, max_size_of_traj
        namelist /grill_parameters/ Zplus, Zminus, ZGapPlus, ZGapMinus, ntet, nnz
        namelist /spectrum/ spectrum_type, spectrum_PWM, spectrum_coordinate_system 
        ! Namelist definition===============================

        call open_inputfile(file_path, file_unit, iostat)
        if (iostat /= 0) then
            print *, 'error open file'
            !! write here what to do if opening failed"
            return
        end if

        read (nml=physical_parameters, iostat=iostat, unit=file_unit)
        read (nml=parameters_for_alphas_calculations, iostat=iostat, unit=file_unit)
        read (nml=numerical_parameters, iostat=iostat, unit=file_unit)
        read (nml=options, iostat=iostat, unit=file_unit)
        read (nml=grill_parameters, iostat=iostat, unit=file_unit)
        read (nml=spectrum, iostat=iostat, unit=file_unit)

        iw_start = iw
        pabs0 = pabs
        eps_const = eps
        i_pol = ipol
        nmaxm(:) = (/nmaxm_1, nmaxm_2, nmaxm_3, nmaxm_4/)

        call close_inputfile(file_path, file_unit, iostat)
        if (iostat /= 0) then
            print *, 'error close file'
            !! write here what to do if reading failed"
            return
        end if
    end subroutine read_nml_parameters

    !! Namelist helpers

    subroutine open_inputfile(file_path, file_unit, iostat)
        !! Check whether file exists, with consitent error message
        !! return the file unit
        character(len=*),  intent(in)  :: file_path
        integer,  intent(out) :: file_unit, iostat

        inquire (file=file_path, iostat=iostat)
        if (iostat /= 0) then
            write (stderr, '(3a)') 'Error: file "', trim(file_path), '" not found!'
        end if
        open (action='read', file=file_path, iostat=iostat, newunit=file_unit)
    end subroutine open_inputfile

    subroutine close_inputfile(file_path, file_unit, iostat)
        !! Check the reading was OK
        !! return error line IF not
        !! close the unit
        character(len=*),  intent(in)  :: file_path
        character(len=1000) :: line
        integer,  intent(in) :: file_unit, iostat

        if (iostat /= 0) then
            write (stderr, '(2a)') 'Error reading file :"', trim(file_path)
            write (stderr, '(a, i0)') 'iostat was:"', iostat
            backspace(file_unit)
            read(file_unit,fmt='(A)') line
            write(stderr,'(A)') &
                'Invalid line : '//trim(line)
        end if
        close (file_unit)   
    end subroutine close_inputfile


    subroutine read_parameters(file_name)
        implicit none
        integer, parameter :: iunit = 20
        character(*) file_name
        print *, file_name

        open(iunit, file= file_name)
        !!!!!!!!!!!!!  read  physical parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               read(iunit,*)
               read(iunit,*) freq
               read(iunit,*) xmi1
               read(iunit,*) zi1
               read(iunit,*) xmi2
               read(iunit,*) zi2
               read(iunit,*) dni2
               read(iunit,*) xmi3
               read(iunit,*) zi3
               read(iunit,*) dni3
        !!!!!!!!!!!!!  read parameters for alphas calculation !!!!!!!!!!!!!!!!!!!
               read(iunit,*)
               read(iunit,*) itend0
               read(iunit,*) energy
               read(iunit,*) factor
               read(iunit,*) dra
               read(iunit,*) kv
        
        !!!!!!!!!!!!!  read  numerical parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               read(iunit,*)
               read(iunit,*) nr
               read(iunit,*) hmin1
               read(iunit,*) rrange
               read(iunit,*) eps
               eps_const = eps
               read(iunit,*) hdrob
               read(iunit,*) cleft
               read(iunit,*) cright
               read(iunit,*) cdel
               read(iunit,*) rbord
               read(iunit,*) pchm
               read(iunit,*) pabs0
               read(iunit,*) pgiter
               read(iunit,*) ni1
               read(iunit,*) ni2
               read(iunit,*) niterat
               read(iunit,*) nmaxm(1)
               read(iunit,*) nmaxm(2)
               read(iunit,*) nmaxm(3)
               read(iunit,*) nmaxm(4)
               read(iunit,*) maxstep2
               read(iunit,*) maxstep4
        
        !!!!!!!!!!!!!  read  options !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               read(iunit,*)
               read(iunit,*) ipri
               read(iunit,*) iw
               read(iunit,*) ismth
               read(iunit,*) ismthalf
               read(iunit,*) ismthout
               read(iunit,*) inew
        
               read(iunit,*) itor     !Btor direction in right-hand {drho,dteta,dfi}
               read(iunit,*) i_pol     !Bpol direction in right-hand {drho,dteta,dfi}
        
        !!!!!!!!!!!!!  read grill parameters and input LH spectrum !!!!!!!!!!!!
               read(iunit,*)
               read(iunit,*) zplus
               read(iunit,*) zminus
               read(iunit,*) ntet
               read(iunit,*) nnz
               read(iunit,*) spectrum_type
        close(iunit)        

        print *, 'checking initial parameters'
        !!!!!!!!!!!!! checking initial parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(kv.gt.50) kv=50
        if(nr.gt.100) nr=100
        if(ni1.eq.0) ni1=20
        if(ni2.eq.0) ni2=20
        if(ni1+ni2.gt.100) then
            ni1=60
            ni2=40
        end if
        if(nnz*ntet.gt.10000) then
            nnz=250
            ntet=40
            pause 'nnz and ntet changed, because nnz*ntet>10000'
        end if

        call show_parameters
    end subroutine read_parameters
end module rt_parameters