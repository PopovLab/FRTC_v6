module nr_grid
    use kind_module
    use maxwell, only: i0
    implicit none

    integer, parameter :: MAX_NR = 200
    !! max size of nr-grids

    real(wp) :: dql(101,MAX_NR)
    !!
    real(wp) :: pdl(MAX_NR)
    real(wp) :: pdc(MAX_NR)
    real(wp) :: pda(MAX_NR)
    !!
    real(wp) :: vzmin(MAX_NR)
    !!
    real(wp) :: vzmax(MAX_NR)
    !common /a0i3/ dql(101,100),pdl(100),vzmin(100),vzmax(100)
    real(wp) :: fcoll(MAX_NR)
    real(wp) :: dens(MAX_NR) 
    real(wp) :: eta(MAX_NR)
    !common /a0i4/ fcoll(100),dens(100),eta(100)
    real(wp) :: dq1(101,MAX_NR)
    real(wp) :: dq2(101,MAX_NR)

    real(wp) :: ppv1,ppv2
    !common/vvv1/dq1(101,100),dq2(101,100),pdc(100),pda(100),ppv1,ppv2
    real(wp) :: pdfast(MAX_NR)
    !common /vvv3/ pdfast(100)
    real(wp) :: dqi0(50,MAX_NR) 
    !common /alph/ dqi0(50,100)    
    real(wp) :: dncount(101,MAX_NR)
    !common/findsigma/dncount(101,100)
    real(wp) :: pdprev1(MAX_NR), pdprev2(MAX_NR)
    !! массивы для невзязки

    real(wp) :: source(MAX_NR)
    !! это относится к альфа-частицам

    !! сетки для максвела
    real(wp) v_grid(i0,MAX_NR)
    !! сетка обычных скоростей
    real(wp) vij(i0, MAX_NR), fij0(i0, MAX_NR,2), fij(i0, MAX_NR,2)
    real(wp) dfij(i0, MAX_NR,2), dij(i0, MAX_NR,2)

contains

subroutine init_nr_grid_arrays(cltn)
    use constants, only: zero
    implicit none
    real(wp), intent(in) :: cltn
    ppv1=zero
    ppv2=zero

    pdl=zero
    pdc=zero
    pda=zero
    pdfast=zero

    dql=zero
    dq1=zero
    dq2=zero
    dqi0=zero
    dncount=zero
    vzmin=cltn
    vzmax=-cltn

end subroutine

function find_nevyazka() result(pchg)
    !! find nevyazka
    use constants, only: zero
    use rt_parameters, only : nr
    implicit none
    !real(wp), intent(inout) :: pdprev1(:), pdprev2(:)
    real(wp) :: psum1, psum2, pchg, pchg1, pchg2
    real(wp) :: dpw1, dpw2
    integer j

    psum1=zero
    psum2=zero
    pchg=zero
    pchg1=zero
    pchg2=zero
    do j=1,nr
        dpw1=pdl(j)+pdc(j)
        dpw2=pda(j)
        psum1=psum1+dpw1**2
        psum2=psum2+dpw2**2
        pchg1=pchg1+(dpw1-pdprev1(j))**2
        pchg2=pchg2+(dpw2-pdprev2(j))**2
        pdprev1(j)=dpw1
        pdprev2(j)=dpw2
    end do
    if(psum1.ne.zero) pchg=pchg1/psum1 !sav2008
    if(psum2.ne.zero) pchg=pchg+pchg2/psum2    

end function


subroutine calculate_total_current_and_power(ol, oc, oa, of)
    !! calculate total current and power
    use constants, only: zero
    use rt_parameters, only : nr
    implicit none
    real(wp), intent(inout):: ol, oc, oa, of
    real(wp) :: cppl, cppc, cppa, cppf
    integer j
    !----------------------------------------
    !     calculate total current and power
    !----------------------------------------
    cppl=zero
    cppc=zero
    cppa=zero
    cppf=zero
    do j=1,nr
        cppl=cppl+pdl(j)
        cppc=cppc+pdc(j)
        cppa=cppa+pda(j)
        cppf=cppf+pdfast(j)
    end do
    ol=cppl*1d-6
    oc=cppc*1d-6
    oa=cppa*1d-6
    of=cppf*1d-6    
end 


subroutine renormalisation_power
    !! renormalisation on xwtt 1e-7_wp
    use constants, only: xwtt
    use rt_parameters, only : nr
    implicit none
    integer j
    do j=1,nr
        pdl(j)=pdl(j)*xwtt
        pdc(j)=pdc(j)*xwtt
        pda(j)=pda(j)*xwtt
        pdfast(j)=pdfast(j)*xwtt
    end do
end


subroutine find_achieved_radial_points(nvpt)
    !!  find achieved radial points jbeg-jend
    use rt_parameters, only : nr
    implicit none
    integer, intent(in) :: nvpt
    integer i, j, jbeg, jend, nvmin, nvach

    nvmin=1 !minimum counted events at a given radius rho
    jbeg=1
    jend=0
    do j=1,nr
        nvach=0
        do i=1,nvpt
            nvach = nvach + dncount(i,j)
        end do
        if (nvach.lt.nvmin) then
            if (jend.eq.0) jbeg = jbeg + 1
        else
            jend=j
        end if
    end do
    if (jend.eq.0.or.jbeg.ge.jend) then
        write(*,*)'failure: jbeg=',jbeg,' jend=',jend 
        pause
        stop
    end if
end subroutine  

end module nr_grid