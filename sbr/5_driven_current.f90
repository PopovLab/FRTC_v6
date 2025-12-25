    subroutine driven_current(outj, sigmaj, UPL, NRD, NA1, TIME, TAU, ROC, RTOR, GP2)
        !! ******************************************************************
        !!   outj(i)  = LH driven current density, MA/m^2
        !!   dndt(i)  = d^2Jr1/dt^2/E, MA/m^2/sec^2/(V/m), ~runaway d(el.density)/dt/E
        !!   djdt(i)  = dJr2/dt, time drivative of runaway current Jr2, MA/m^2/sec
        !!   outjrun(i)  = LH driven runaway current density, MA/m^2
        !!   outnerun(i) = runaway electron density/10^19 m^-3
        !! ******************************************************************
        use FokkerPlanck_module
        use driven_current_module
        implicit none
        ! variables imported from ASTRA
        integer NRD
        ! NRD     501 - Maximum size of the radial grid
        integer NA1
        ! NA1     Edge grid point number: ROC=RHO(NA1)      
        double precision TIME, TAU, RTOR, ROC, GP2
        double precision UPL(NRD)

        real*8 outj(NRD),sigmaj(NRD),afld(NRD),dtau
        integer i,inpt,ispectr

        real*8 dt, cup,cup0,cum,cum0,cp,cm,cp0,cm0,aiint
        real*8, parameter :: zero=0.d0, eps=1.d-2 
        type(DrivenCurrentResult) :: rc_result
        type(DrivenCurrent) :: positive_dc
        type(DrivenCurrent) :: negative_dc
    
        interface
            subroutine lhcurrent(driven_current, ispectr)
                use driven_current_module
                implicit none
                type(DrivenCurrent), intent(inout) :: driven_current
                integer,             intent(in)    :: ispectr
            end subroutine lhcurrent
        end interface
        inpt=NA1

        ! ---- starting LH current calculation 
        positive_dc = DrivenCurrent(NA1)
        negative_dc = DrivenCurrent(NA1)

        ! ---- positive spectrum:
        call lhcurrent(positive_dc, ispectr= 1)
        call positive_dc%evaluate(ROC)

        ! ---- negative spectrum:
        call lhcurrent(negative_dc, ispectr= -1)
        call negative_dc%evaluate(ROC)


        do i=1,inpt
            outj(i) = positive_dc%outj(i) + negative_dc%outj(i)
            sigmaj(i)=zero
            if (abs(afld(i)).gt.eps) then
                sigmaj(i) = (positive_dc%ohj(i) + negative_dc%ohj(i))/afld(i)
            end if
            !!!!       write(*,*) i,outj(i)
        end do

        call save_driven_current(positive_dc, negative_dc, time)

        rc_result = DrivenCurrentResult(positive_dc, negative_dc)
        call rc_result%print(time)
        call rc_result%save(time)

        call fokkerplanck_compute(time, TAU)

    end

    subroutine save_driven_current(pos_dc, neg_dc, time)
        use driven_current_module
        use rt_parameters, only : nr
        implicit none
        type(DrivenCurrent), intent(in) :: pos_dc
        type(DrivenCurrent), intent(in) :: neg_dc
        real(wp),            intent(in) :: time

        integer i
        character(132) fname
        integer :: iu

        write(fname,'("lhcd/driven_current/", f9.7,".dat")') time
        open(newunit=iu, file=FNAME, status="replace", action="write")
	    write(iu,'(10A16)'), 'index', 'pos_dc', 'neg_dc'


        do i=1, nr
            write (iu, *) i, pos_dc%current(i), neg_dc%current(i)
        end do

        !do i=1, pos_dc%grid_size
        !    write (iu, *) i, pos_dc%outj(i), neg_dc%outj(i)
        !end do

        close(iu)

    end subroutine

    subroutine lhcurrent(driven_current, ispectr)
        !subroutine lhcurrent(outj,ohj,cuj,cujoh,inpt,ispectr)
        !!      implicit real*8 (a-h,o-z)
        use plasma, only : rh, rh1, fn1,fn2, fvt, sk
        use lock_module
        use maxwell
        use nr_grid, only: vij, fij0, fij
        use rt_parameters, only : nr, inew
        use driven_current_module
        use math_module        
        implicit none
        type(DrivenCurrent), intent(inout) :: driven_current
        integer,             intent(in)    :: ispectr
        !real*8 outj(*),ohj(*),
        real*8 cuj,cujoh,curs,curs0,curdir
        real*8 currn,vt0,ccur,cfull,cfull0
        real*8 r,pn,fnr,fnrr,vt,vto!,rh1
        integer klo,khi,ierr,nrr,i,j,inpt,ismthout
        !common /a0ab/ nr
        !real*8 y2dn,y2tm,y2tmi
        !common /a0l3/ y2dn(501),y2tm(501),y2tmi(501)
        !integer inew
        !common /cnew/ inew !est !sav2008
        real(wp) :: fout

        !real*8 zv1,zv2
        !common/plosh/ zv1(100,2),zv2(100,2)!,sk(100)
        integer k
        !parameter(i0=1002)
        !real*8 vij,fij0,fij,dfij,dij,enorm,fst
        !common/lh/ vij(i0,100), fij0(i0,100,2), fij(i0,100,2), dfij(i0,100,2), dij(i0,100,2), enorm(100), fst(100)
        real*8,dimension(:),allocatable:: vj, fj, fj0, cur, cur0, currnt, rxx, wrk
        parameter(ismthout=1)

        inpt = driven_current%grid_size
        allocate(vj(i0),fj(i0),fj0(i0),cur(nr),cur0(nr),currnt(nr+2),rxx(nr+2),wrk(nr+2))
        !---------------------------------------------------
        ! initial constants
        !---------------------------------------------------
        !pqe=4.803e-10
        vt0 = fvt(zero)
        ccur = pqe*vt0*0.333d-9
        curdir = -dble(ispectr)

        cfull = zero
        cfull0 = zero
        k = (3-ispectr)/2
        do j=1, nr
            do i=1, i0
                vj(i)  = vij(i,j) !Vpar/Vt
                fj0(i) = fij0(i,j,k)
                fj(i)  = fij(i,j,k)-fij0(i,j,k)
            end do
            r=dble(j)/dble(nr+1)
            if(inew.eq.0) then !vardens
                pn=fn1(r,fnr)
            else
                pn=fn2(r,fnr,fnrr)
            end if
            vt  = fvt(r)
            vto = vt/vt0
            curs   = currlhcd(vj,fj)
            cur(j) = curs*pn*ccur*curdir*vto  !Ampere/cm2
            cfull  = cfull+cur(j)*sk(j)

            curs0   = currlhcd(vj,fj0)          
            cur0(j) = curs0*pn*ccur*curdir*vto  !Ampere/cm2
            cfull0  = cfull0+cur0(j)*sk(j)
        end do

        !cuj=cfull*1d-6   !driven current, MA
        driven_current%cu = cfull*1d-6

        !cujoh=cfull0*1d-6   !driven current, MA
        driven_current%cu0 = cfull0*1d-6 

        !!      write(*,*)
        !!      write(*,*)'ccur',ccur,' curdir=',curdir,' nr=',nr
        !!      write(*,*)'cu_out, MA=',cu_out,' cfull, A=',cfull

        do j=1,nr
            driven_current%current(j)=cur(j)*1.d-2    ! Jstoped, MA/m^2
        end do

        currn=cur(1)                   ! Jstoped, A/cm^2
        currnt(1)=currn*1.d-2          ! Jstoped, MA/m^2
        rxx(1)=zero
        do j=1,nr
            rxx(j+1)=dble(j)/dble(nr+1)
            currn=cur(j)                   ! Jstopped, A/cm^2
            currnt(j+1)=currn*1.d-2        ! Jstoped, MA/m^2
        end do
        nrr=nr+2
        rxx(nrr)=1.d0
        currnt(nr+2)=zero

        if(ismthout.ne.0) then
            do i=1,nrr
                wrk(i)=currnt(i)
            end do
            call fsmoth4(rxx,wrk,nrr,currnt)
        end if

        rh(1)=rh1
        if(rh(inpt).gt.1d0) rh(inpt)=1.d0
        do j=1,inpt
            call lock2(rxx,nrr,rh(j),klo,khi,ierr)
            if(ierr.ne.0) then
                write(*,*)'lock2 error in current profile for ASTRA'
                write(*,*)'ierr=',ierr,' j=',j,' rh(j)=',rh(j)
                write(*,*)'rxx(1)=',rxx(1),' rxx(nrr)=',rxx(nrr)
                pause
            end if
            call linf(rxx,currnt,rh(j),fout,klo,khi)
            driven_current%outj(j)=fout
        end do
        rh(1)=zero

        !------------------------------------------------------
        currn=cur0(1)                   ! Jstoped, A/cm^2
        currnt(1)=currn*1.d-2           ! Jstoped, MA/m^2
        rxx(1)=zero
        do j=1,nr
            rxx(j+1)=dble(j)/dble(nr+1)
            currn=cur0(j)                  ! Jstopped, A/cm^2
            currnt(j+1)=currn*1.d-2        ! Jstoped, MA/m^2
        end do
        nrr=nr+2
        rxx(nrr)=1.d0
        currnt(nr+2)=zero

        if(ismthout.ne.0) then
            do i=1,nrr
                wrk(i)=currnt(i)
            end do
            call fsmoth4(rxx,wrk,nrr,currnt)
        end if

        rh(1)=rh1
        if(rh(inpt).gt.1d0) rh(inpt)=1.d0
        do j=1,inpt
            call lock2(rxx,nrr,rh(j),klo,khi,ierr)
            if(ierr.ne.0) then
                write(*,*)'#2 lock2 error in current profile for ASTRA'
                write(*,*)'ierr=',ierr,' j=',j,' rh(j)=',rh(j)
                write(*,*)'rxx(1)=',rxx(1),' rxx(nrr)=',rxx(nrr)
                pause
            end if
            call linf(rxx,currnt,rh(j),fout,klo,khi)
            driven_current%ohj(j)=fout
        end do
        rh(1)=zero

        deallocate(vj,fj,fj0,cur,cur0,currnt,rxx,wrk)
      end

