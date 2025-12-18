      subroutine lhcd2017(out_lh_power)
! ******************************************************************
!   outj(i)  = LH driven current density, MA/m^2
!   outpe(i) =LH power density (Landau+coll.) deposited into electrons, MW/m^3
!   outpe = out_lh_power
!   outpec(i) = LH power density (collisions) deposited into electrons, MW/m^3
!   outpef(i) = LH power dens. dep. into el. by fast wave, MW/m^3
!   dndt(i)  = d^2Jr1/dt^2/E, MA/m^2/sec^2/(V/m), ~runaway d(el.density)/dt/E
!   djdt(i)  = dJr2/dt, time drivative of runaway current Jr2, MA/m^2/sec
!   outpa(i)  = LH power density deposited into alphas, MW/m^3
!   outda(i)  = (Na/Ne) relative density of alpha-particles
! ******************************************************************
      use approximation
      use utils
      use plasma
      use rt_parameters
      use maxwell  
      use spectrum_mod    
      use lhcd_module  
      implicit none
      integer i
      real*8 p_in!,pe_p,pe_m,c_p,c_m
      !real*8 vint
      include 'for/parameter.inc'
      include 'for/const.inc'
      include 'for/status.inc'
      real*8 out_lh_power(NRD)
      real(wp),dimension(:),allocatable:: outpep,outpem
      type(Spectrum) spectr
      type(Timer) my_timer
! *********************************************************************
!     Co-ordinates used in ray-tracing:
!          (x-x0)/rm=r*cos(teta)-delta-gamma*sin^2(teta) !sav2008 -gamma
!          (z-z0)/rm=ell*r*sin(teta)
!     Definitions:
!     (x0,z0) - magnetic axis position, centimeters
!     rm      - minor radius in mid-plane, cenrimeters
!     r(rho_ASTRA),delta(r),gamma(r),ell(r) - dimensionless functions
!     rho_ASTRA=sqrt(Phi_tor/GP/BTOR)
!     Interval for r:  0.<= r <=1.
! *********************************************************************
      call my_timer%start('lhcd/lhcd_time.dat', time)
      print *, 'start lhcd2017 time=', time
      print *, 'NA1 NB1', NA1, NB1
      tcur=time
      out_lh_power=zero
      p_in=dble(QLH)    ! input LH power, MW

      if(p_in.eq.zero) then 
            dij(:,:,:)=zero
            return
      end if

      !call read_parameters('lhcd/ray_tracing.dat')
      call read_nml_parameters('lhcd/parameters.nml')
      call show_parameters()
      !pause
      call init_plasma(NA1,ABC,BTOR,RTOR,UPDWN,GP2,
     & AMETR,RHO,SHIF,ELON,TRIA,MU,NE,TE,TI,ZEF,UPL)

      !call write_plasma(time)
      call write_equilibrium(time)
      call write_lcms

      full_spectrum = read_spectrum('lhcd/spectrum.dat')
      full_spectrum%input_power = p_in

      pos_spectr = full_spectrum%get_positive_part()
      neg_spectr = full_spectrum%get_negative_part()
      
      call full_spectrum%normalization
      call full_spectrum%write('full_spectrum')

      allocate(outpep(ngrid),outpem(ngrid))

      !!positive spectrum:
      print *, 'positive spectrum'
      if(pos_spectr%input_power > zero) then
            call eval_lhcd(pos_spectr, 'spectrum_pos', outpep)  
      else
            dij(:,:,1)=zero
            outpep = 0
      end if      

      !!negative spectrum:
       print *, 'negative spectrum'
       if(neg_spectr%input_power > zero) then        
            call eval_lhcd(neg_spectr, 'spectrum_neg', outpem)  
       else
            dij(:,:,2)=zero
            outpem = 0
       endif     

      do i=1,ngrid
            out_lh_power(i)=outpep(i)+outpem(i)
      end do

      deallocate(outpep,outpem)
      call my_timer%stop_and_save
      print *, 'lhcd elapsed time', my_timer%elapsed_time
      !pause

      end


      subroutine eval_lhcd(in_spectrum, spectrum_name, out_power)
            use constants, only: zero
            use plasma, only: ngrid
            use rt_parameters, only: spectrum_type
            use spectrum_mod !, only: make_spline_approximation
            use lhcd_module, only: ourlhcd2017       
            implicit none
            include 'for/parameter.inc'
            include 'for/const.inc'
            !include 'for/status.inc'            
            type(Spectrum),   intent(in) :: in_spectrum
            character(len=*), intent(in) :: spectrum_name
            real(wp),         intent(inout) :: out_power(*)

            type(Spectrum) spectr
            integer :: i

            real(wp) pe, vi_pe
            real(wp) vint

            pe=0
            do i=1,ngrid
                  out_power(i)= 0
            end do

            select case (spectrum_type)
            case (0)
                  spectr = make_spline_approximation(in_spectrum)
            case (1)
                  spectr = in_spectrum 
                  call spectr%normalization
            case (2)
                  spectr = in_spectrum
                  call spectr%normalization
            case (3)
                  print *, '2D spectrum'
                  stop
            end select            

            call spectr%write(spectrum_name) !'spectrum_neg')            
            call ourlhcd2017(spectr, out_power, pe)              
             
             if(pe.ne.zero) then
                  vi_pe=vint(out_power, roc) ! define in stdfun.f
                  ! VINT:	Volume integral {0,R} of any array
                  ! Only a radially dependent array may be the 1st parameter of the function
                  ! Examples:
                  !    out\Vint(CAR3)	!Radial profile of CAR3 volume integral
                  !    out_Vint(CAR3,Ro); !Volume integral {0,Ro} of CAR3
                  !    out_Vint(CAR3B)    !Total volume integral of CAR3 (0,ROC)
                  !			(Yushmanov 26-DEC-90)
                  if(vi_pe.ne.zero) then
                        do i=1,ngrid
                              out_power(i)= pe*out_power(i)/vi_pe
                        end do
                  end if
            end if
      end subroutine