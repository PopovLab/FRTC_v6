module savelyev_solver_module
    use kind_module
    use chang_cooper_module, only: tridag
    implicit none
    
    PRIVATE :: q, k, d
contains
    subroutine savelyev_solver(alfa2, nt, h, dt, n, ybeg, yend, d1,d2,d3, y)
        !! разностная схема Савельева для уравнения Фоккера-Планка
        implicit none
        real(wp), intent(in)  :: alfa2      
        integer, intent(in)   :: nt, n
        real(wp), intent(in)  :: h, dt
        real(wp), intent(in)  :: ybeg, yend
        real(wp), intent(in)  :: d1(n+1),d2(n+1),d3(n+1)
        real(wp), intent(inout) :: y(n)

        integer i, it
        real(wp) xx(n+1), a(n),b(n),c(n),f(n)

        do i=1,n+1
            xx(i) = h/2.d0 + h*dble(i-1) !+shift
        end do
    
        do it=1,nt
            call savelyev_abccoef(alfa2, a,b,c,f,y,dt,n,ybeg,yend,xx,h,d1,d2,d3)
            call tridag(a,b,c,f,y,n)      
        end do

    end subroutine

    
    subroutine savelyev_abccoef(alfa2, a,b,c, f, y, dt, n, ybeg, yend, xx, h, d1,d2,d3)
        !! -- fill abc matrix
        implicit none
        real(wp), intent(in)    :: alfa2
        real(wp), intent(inout) :: a(n),b(n),c(n),f(n),y(n)
        real(wp), intent(in)    :: dt
        integer,  intent(in)    :: n
        real(wp), intent(in)    :: ybeg, yend, h
        real(wp), intent(in)    :: xx(n+1)
        real(wp), intent(in)    :: d1(n+1),d2(n+1),d3(n+1)

        integer i,iunit,iunit2
        real(wp) a1(n),b1(n),c1(n),f1(n),a2(n),b2(n),c2(n),f2(n)
        !real(wp) kinv,rs,rmink,rplusk,q,qf,r1,rmink2,rplusk2,kinv2
        real(wp) r,kappa,sum,bmin,bplus,sum2,sum3,sum4
        real(wp) dc,as(n+1)
        !external kinv,rs,rmink,rplusk,q,kinv2,rmink2,rplusk2,d

        sum = (kinv(xx(1) - h/2d0, d2(1)) + kinv(xx(1) + h/2d0, d3(1)))*h/2d0
        as(1) = h/sum

        sum = (kinv(xx(2)-h/2d0, d2(2))+kinv(xx(2)+h/2d0, d3(2)))*h/2d0
        as(2) = h/sum

        r = h/2d0*dabs(rs(xx(1)+h/2d0, alfa2))/k(xx(1)+h/2d0, d3(1))
        kappa = 1d0/(1d0+r)
        sum = (rmink(xx(1), d1(1), alfa2) + rmink(xx(2), d1(2), alfa2))*h/2d0
        bmin = sum/h

        sum = (rplusk(xx(1), d1(1), alfa2) + rplusk(xx(2), d1(2), alfa2))*h/2d0
        bplus = sum/h

        sum = qf(xx(2))-qf(xx(1))
        dc = sum/h

        a(1) = as(1)*(kappa/h**2 - bmin/h)
        c(1) = as(2)*(kappa/h**2 + bplus/h)
        b(1) = -(1d0/dt+a(1) + c(1) + dc)
        f(1) = -y(1)/dt-a(1)*ybeg
        do i=2,n
            sum = (kinv(xx(i+1) - h/2d0, d2(i+1)) + kinv(xx(i+1) + h/2d0, d3(i+1)))
            sum = sum*h/2d0
            as(i+1) = h/sum

            r = h/2d0*dabs(rs(xx(i) + h/2d0, alfa2))/k(xx(i) + h/2d0, d3(i))
            kappa=1d0/(1d0+r)
            sum = (rmink(xx(i), d1(i), alfa2) + rmink(xx(i+1), d1(i+1), alfa2))*h/2d0
            bmin = sum/h
            sum = (rplusk(xx(i), d1(i), alfa2) + rplusk(xx(i+1), d1(i+1), alfa2))*h/2d0
            bplus = sum/h
            sum = qf(xx(i+1)) - qf(xx(i))
            dc = sum/h
            
            a(i) = as(i)*(kappa/h**2 - bmin/h)
            c(i) = as(i+1)*(kappa/h**2 + bplus/h) 
            b(i) = -(1d0/dt + a(i) + c(i) + dc) 
            f(i) = -y(i)/dt
        end do
        f(n) = f(n)-c(n)*yend
        a(1) = 0d0
        c(n) = 0d0
    end

    real(wp) function rplusk(x, dif, alfa2)
        implicit none
        real(wp), intent(in)    :: x,dif
        real(wp), intent(in)    :: alfa2      
        rplusk=0.5d0*(rs(x, alfa2)+dabs(rs(x, alfa2)))/k(x,dif)
    end

    real(wp) function rplusk2(x, dif, alfa2)
        implicit none
        real(wp), intent(in)    :: x,dif
        real(wp), intent(in)    :: alfa2      
        rplusk2=0.5d0*(rs(x, alfa2)+dabs(rs(x, alfa2)))/k2(x,dif)
    end

    real(wp) function rmink(x, dif, alfa2)
        implicit none
        real(wp), intent(in)    :: x,dif
        real(wp), intent(in)    :: alfa2      
        rmink = 0.5d0*(rs(x, alfa2)-dabs(rs(x, alfa2)))/k(x,dif)
    end

    real(wp) function rmink2(x, dif, alfa2)
        implicit none
        real(wp), intent(in)    :: x,dif
        real(wp), intent(in)    :: alfa2      
        rmink2=0.5d0*(rs(x, alfa2)-dabs(rs(x, alfa2)))/k2(x,dif)
    end

    real(wp) function rs(x, alfa2)
        implicit none
        real(wp), intent(in) :: x
        real(wp), intent(in) :: alfa2
        !common/ef/ alfa2
        rs=1d0/x**2-alfa2
    end

    real(wp) function q(x)
        implicit none
        real(wp), intent(in) :: x
        q=2d0/x**3
    end

    real(wp) function qf(x)
        implicit none
        real(wp), intent(in) :: x
        qf=-1d0/x**2
    end

    real(wp) function k(x,dif)
        implicit none
        real(wp), intent(in) :: x, dif
        k=dif+1d0/x**3
    end

    real(wp) function k2(x,dif)
        implicit none
        real(wp), intent(in) :: x,dif
        k2=d(x)+1d0/x**3
    end

    real(wp) function kinv(x,dif)
        implicit none
        real(wp), intent(in) :: x,dif
        kinv=x**3/(dif*x**3 + 1d0)
    end

    real(wp) function kinv2(x,dif)
        implicit none
        real(wp), intent(in) :: x,dif
        kinv2=x**3/(d(x)*x**3 + 1d0)
    end





    real(wp) function d(x)
        !! возможно одна из самых замедляющих функций
        use maxwell, only: i0, flag_d0, jindex, kindex
        use nr_grid, only: vij, dij
        use constants, only : zero
        use lock_module
        implicit none
        !integer i0
        !parameter(i0=1002)
        !real*8 vij,fij0,fij,dfij,dij,enorm,fst
        !common/lh/vij(i0,100),fij0(i0,100,2),fij(i0,100,2),dfij(i0,100,2),dij(i0,100,2),enorm(100),fst(100)
        real(wp), dimension(:), allocatable:: vvj,ddj
        integer klo,khi,ierr
        real(wp)  x
        integer k,j,i
        !common/dddql/ d0,jindex,kindex

        d=zero
        if(flag_d0) return
        j=jindex
        if(x.ge.vij(i0,j)) return
        k=kindex

        allocate(vvj(i0),ddj(i0))

        do i=1,i0
            vvj(i)=vij(i,j)
            ddj(i)=dij(i,j,k)
        end do

        call lock(vvj,i0,x,klo,khi,ierr)
        if(ierr.eq.1) then
            write(*,*)'lock error in finction d2(x)'
            write(*,*)'j=',j,' v=',x
            write(*,*)'vj(1)=',vvj(1),' vj(i0)=',vvj(i0)
            pause
            stop
        end if
        d=ddj(klo)
        deallocate(vvj, ddj)
    end

end module savelyev_solver_module