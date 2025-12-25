module FokkerPlanck_module
    use kind_module
    implicit none
    
contains
    !! calculation of distribution functions at time t1=t+dtau !!
subroutine fokkerplanck_compute(time, TAU)
    use FokkerPlanck1D_mod, only: FokkerPlanck1D
    use rt_parameters, only: nr
    use writer_module, only: write_v_array, binary_write_array
    use maxwell, only: jindex, kindex, flag_d0
    use nr_grid, only: vij, fij, fij0, dij, dfij
    use plasma, only : fvt, enorm, fst
    implicit none

    type(FokkerPlanck1D) fokker_planck
    real(wp), intent(in) :: time, TAU
    real(wp) t, dtstep, dtau
    integer, parameter :: ntau = 10
    integer n,i,j,it,nt,k
    real(wp) xend,h,dt
    real(wp) znak,alfa2,dt0,h0,r

    parameter(dt0=0.1d0,h0=0.1d0)

    dtstep=TAU/dble(ntau) !seconds 

    print *, 'fokkerplanck_compute'
    write(*,*)'time=',time,' dt=',dtstep

    do j=1, nr
        jindex=j  !common/dddql/ 
        dtau=dtstep*fst(j)
        nt=1
        if(dtau.gt.dt0) then
            nt=1+dtau/dt0
        end if
        dt=dtau/nt
        r=dble(j)/dble(nr+1)
        xend=3.d10/fvt(r)
        do k=1,2
            kindex=k
            flag_d0=.TRUE. ! d(x) enable
            znak=2.d0*dble(k)-3.d0 ! znak = -1 if k = 1 or znak = 1 if k = 2
            fokker_planck = FokkerPlanck1D(znak*enorm(j), xend, vij(:,j), fij0(:,j,k))
            call fokker_planck%init_zero_diffusion
            do i=1, ntau
                call fokker_planck%solve_time_step(dt, nt)
                !call fokkerplanck1D_iter(alfa2, h, n, dt, nt, xend, d1, d2, d3, vij(:,j), fij0(:,j,k), out_fj)
            end do
            fij0(:,j,k) = fokker_planck%f

            flag_d0=.FALSE. ! d(x) disable
            fokker_planck = FokkerPlanck1D(znak*enorm(j), xend, vij(:,j), fij(:,j,k))
            call fokker_planck%init_diffusion(dij(:,j,k))
            do i=1, ntau
                call fokker_planck%solve_time_step(dt, nt)
                !call fokkerplanck1D_iter(alfa2, h, n, dt, nt, xend, d1, d2, d3, vij(:,j), fij(:,j,k),out_fj, dfij(:,j,k))
            end do
            fij(:,j,k) = fokker_planck%f
            !call fokker_planck%eval_f_derivate(dfij(:,j,k))
        end do
    end do

    write(*,*)'fokkerplanck nr= ',nr,' ntau =',ntau, 'nt =', nt

    call binary_write_array(vij, fij0(:,1:nr,:), time, 'maxwell_fij0')
    call write_v_array(vij, fij(:,1:nr,:),  time, 'maxwell')
    call write_v_array(vij, dfij(:,1:nr,:), time, 'f_derivative')
    call write_v_array(vij,  dij(:,1:nr,:), time, 'diffusion')
    !call write_matrix(dij(1:i0,1:nr,1), time, 'diffusion')
    !time2 = sys_time() - time1
    !call my_timer%stop
    !print *, 'fokkerplanck elapsed_time: ', my_timer%elapsed_time    
   ! print *, 'fokkerplanck_new eval time: ', time2    
    !pause
    
 end


 subroutine init_diffusion(h, n, vj, dj, d1, d2, d3)
    ! инициализация диффузии для схемы савельева
    use lock_module
    implicit none
    integer, intent(in) :: n
    real(wp), intent(in) :: h
    real(wp), dimension(:), intent(in) :: vj, dj
    real(wp), dimension(:), intent(out) :: d1, d2, d3
    real(wp), dimension(:), allocatable :: xx
    integer :: i0
    integer i, klo, khi, ierr, klo1, khi1
    integer klo2, klo3, khi2, khi3, ierr1, ierr2, ierr3

    i0 = size(vj)

    allocate(xx(n+1))
    do i=1,n+1
        xx(i)=h/2.d0+h*dble(i-1) !+shift
    end do

    do i=1,n+1
        call lock(vj, i0, xx(i), klo1, khi1, ierr1)
        call lock(vj, i0, xx(i)-h/2d0, klo2, khi2, ierr2)
        call lock(vj, i0, xx(i)+h/2d0, klo3, khi3, ierr3)
        if(ierr1.eq.1) then
            write(*,*)'lock error in finction d2(x)'
            write(*,*)'j=', 123,' v=', xx(i)
            write(*,*)'klo1=', klo1, 'khi1=', khi1, 'i=', i
            write(*,*)'vj(1)=', vj(1),' vj(i0)=', vj(i0)
            pause
            stop
        end if
        if(ierr2.eq.1) then
            write(*,*)'lock error in finction d2(x)'
            write(*,*)'j=', 123, ' v=', xx(i)
            write(*,*)'klo2=', klo2, 'khi2=', khi2, 'i=',i
            write(*,*)'vj(1)=', vj(1), ' vj(i0)=', vj(i0)
            pause
            stop
        end if
        if(ierr3.eq.1) then
            write(*,*)'lock error in finction d2(x)'
            write(*,*)'j=', 123, ' v=', xx(i)
            write(*,*)'klo3=', klo3, 'khi3=', khi3, 'i=',i
            write(*,*)'vj(1)=', vj(1), ' vj(i0)=', vj(i0)
            pause
            stop
        end if
        d1(i) = dj(klo1)
        d2(i) = dj(klo2)
        d3(i) = dj(klo3)	
    end do

end subroutine
end module FokkerPlanck_module