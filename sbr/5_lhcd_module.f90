module lhcd_module
    !! LHCD модуль
    use kind_module
    use spectrum_mod
    implicit none
    
    type(Spectrum) :: full_spectrum
    type(Spectrum) :: pos_spectr, neg_spectr

contains
    subroutine ourlhcd2017(spectr, outpe, pe_out)      
        use constants, only: zero, xsgs
        use plasma, only: nspl, cltn, rh, r0, con, tcur
        use plasma, only: find_volums_and_surfaces
        use rt_parameters, only: pabs0, ipri, niterat
        use rt_parameters, only: nr, kv, ntet, pgiter, itend0
        use rt_parameters, only: iw, iw_start
        use trajectory_module, only: write_trajectories,  init_trajectories
        use spectrum_mod
        use manager_mod
        use nr_grid, only: ppv1, ppv2
        use nr_grid, only: pdprev1, pdprev2        
        use nr_grid, only: init_nr_grid_arrays, find_nevyazka
        use nr_grid, only: calculate_total_current_and_power
        use nr_grid, only: find_achieved_radial_points, renormalisation_power

        use iteration_result_mod
        use power, only: pnab, plost, psum4
        use small_vgrid, only: nvpt
        use small_vgrid, only: calculate_dfundv
        use small_vgrid, only: find_velocity_limits_and_initial_dfdv, recalculate_f_for_a_new_mesh
        use math_module,  only: integral
        use decrements,   only: kzero
        use source_new_mod
        
        implicit none
        real*8 outpe, pe_out 
        dimension outpe(*)

        type(Spectrum)        :: spectr
        type(IterationResult) :: iteration_result

        integer  :: iterat
        real(wp) :: cn1, avedens
        real(wp) :: anb, fuspow, o_da

        real(wp) :: q_rest, q_abs, q_cond
        real(wp) :: pchg

        real(wp) :: oi
        real(wp) :: ol, oc, oa, of
        real(wp) :: zff, cnyfoc, dconst, fout
        
        real(wp) :: plaun

        integer  :: ispectr
        integer  :: iww, iw0, izz

        plaun = spectr%input_power
        ispectr = spectr%direction
        iw = iw_start
        iw0=iw
    
        call find_volums_and_surfaces

        call init_nr_grid_arrays(cltn)
        call init_iteration_vars

        anb=zero
        fuspow=zero
        o_da=zero
        
        call find_velocity_limits_and_initial_dfdv(anb)
        call calculate_dfundv(ispectr)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(itend0.gt.0) then  ! begin alpha-source renormalisation
            call alpha_source_renormalisation(anb, fuspow)
        end if

        ! ------------------------------------
        !  set initial values of arrays
        ! ------------------------------------
        iww=0
        izz=zero
        kzero=kv

        call init_trajectories(iw0, spectr)
        call init_alphas
        
        ! ----------------------------------------------------------------------
        ! sign of driven current in right coordinate system {dro,dteta,dfi}:
        ! curdir=+1.0 for current drive in positive direction "dfi"
        ! curdir=-1.0 for current drive in negative direction "dfi"
        ! spectrum Nz>0 is along dfi>0 and Nz<0 is along dfi<0
        ! it is also OK if Npar is used instead of Nz, but for Btor>0, 
        ! that is along dfi>0
        ! curdir=-dble(ispectr)
        ! -----------------------------------------------------------------------

        !!!!!!!!!!!!!!! begin iterations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        q_rest=plaun
        iterat=0
        !! массивы для невзязки
        pdprev1=zero
        pdprev2=zero

    80  continue
        call manager(iterat, iw0, ntet, spectr)

        call find_achieved_radial_points(nvpt)

        call renormalisation_power
        
        pchg = find_nevyazka()

        call calculate_total_current_and_power(ol, oc, oa, of)

        !!!!!!!!! prepare to the next iteration !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        iterat=iterat+1
        q_abs=ol+oc+oa
        q_rest=pnab/xsgs
        q_cond=zero
        if(q_abs.ne.zero) q_cond=0.5d0*q_rest/q_abs
        !!!      if(q_cond.le.pabs0.and.pchg.lt.pgiter)

        call integral(1,nspl,rh,con,avedens) 

        iteration_result = IterationResult(number = iterat, &
            spectr_direction = ispectr, P_launched = plaun, &
            P_landau = ol, P_coll = oc, P_alph = oa, &
            alphas_power = fuspow, P_fast = of, &
            P_lost = plost/xsgs, P_not_accounted = pnab/xsgs, &
            P_landau_strong_absorption = ppv1/xsgs, &
            P_landau_weak_absorption = ppv2/xsgs, &
            P_turns = psum4/xsgs, efficiency = oi/plaun, &
            avedens = avedens*1.d19, r0 = r0*1.d-2, &
            eta_eff = 1.d17*avedens*r0*oi/plaun, &
            residual = pchg)

        if(iterat.gt.5.and.q_cond.le.pabs0.and.pchg.lt.pgiter) goto 110

        if(ipri.gt.1) then
            call iteration_result%print
            call iteration_result%save(tcur)
            !pause
        end if

        if(iterat.le.niterat) then
            call recalculate_f_for_a_new_mesh(ispectr, iterat)
            call init_iteration_vars
            call init_nr_grid_arrays(cltn)
            goto 80
        end if
        ! ------------------------------------------
        !  save results
        ! ------------------------------------------
110     continue

        if(ipri.gt.0) then
            write (*,*)
            write (*,*) 'RAY-TRACING RESULTS:'
            call iteration_result%print
            call iteration_result%save(tcur)
            write (*,*) '-------------------------------------------'
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call calculate_diffusion(ispectr)

        call write_trajectories(tcur,ispectr,spectr%size,ntet)  !writing trajectories into a file

        call write_lhcd_power(tcur, ispectr) 

        call calculate_out_power(outpe)

        pe_out=ol+oc
        
    end    
    
    subroutine alpha_source_renormalisation(anb, fuspow)
        use constants, only: zero, talfa, one_third
        use rt_parameters, only: nr, dra, factor
        use source_new_mod, only: rsou, sou, npta
        use plasma, only: fti, dn1, dn2, vk
        use nr_grid, only: dens, source
        use nr_grid, only: MAX_NR
        implicit none
        real(wp), intent(inout) :: anb
        real(wp), intent(inout) :: fuspow
        real(wp) :: r, hr, tt
        real(wp) :: anb0, aratio, sssour
        real(wp) :: ddens, tdens
        real(wp) :: sour(MAX_NR)
        integer j
        hr = 1.d0/dble(nr+1)
        fuspow=anb*talfa*1.6022d-19
        anb0=anb
        anb=zero
        do j=1,nr
            r=hr*dble(j)
            if(r.le.dra) then
                tt=fti(zero)**one_third
            else
                tt=fti(r-dra)**one_third    ! (shifted ti, kev)^1/3
            end if
            ddens=dn1*dens(j)
            tdens=dn2*dens(j)
            sour(j)=4d-12*factor*ddens*tdens*dexp(-20d0/tt)/tt**2
            anb=anb+sour(j)*vk(j)
        end do
        aratio=anb0/anb
        rsou(1)=zero
        sou(1)=aratio*sour(1)
        do j=1,nr
            r=hr*dble(j)
            rsou(j+1)=r
            sou(j+1)=aratio*sour(j)
            if(j.eq.nr) sssour=source(j)
            source(j)=sou(j+1)
        end do
        npta=nr+2
        rsou(npta)=1.d0
        sou(npta)=aratio*sour(nr)        
    end

    subroutine calculate_diffusion(ispectr)
        !! calculate diffusion
        use constants, only: c0, pme, zero
        use rt_parameters, only: nr, inew, ni1, ni2
        use plasma, only: zefff, fn1, fn2
        use plasma, only: vt0, fvt, cltn, cnye
        !use driven_current_module, only : zv1, zv2
        use nr_grid, only : dql
        use maxwell, only: i0
        use nr_grid, only: vij, dfij, dij
        use small_vgrid, only: ipt, ipt1, MAX_PT
        use small_vgrid, only: vrj, dj, vgrid
        use lock_module, only: lock
        implicit none
        integer, intent(in) :: ispectr
        integer  :: i,j,k
        integer  :: klo,khi,ierr
        real(wp) :: r, hr
        real(wp) :: vt, vto, vmax, zff, cnyfoc
        real(wp) :: pn, fnr, fnrr
        real(wp) :: dconst, ddout
        !real(wp) :: dijk(MAX_PT,100,2), vrjnew(MAX_PT,100,2)
        !встречает только один раз common/t01/dijk(101,100,2), vrjnew(101,100,2)
        !
        hr = 1.d0/dble(nr+1)
        k=(3-ispectr)/2
        do j=1,nr
            r=hr*dble(j)
            vt=fvt(r)
            vto=vt/vt0
            vmax=cltn/vto
            zff=(5d0+zefff(r))/5d0
            cnyfoc=zff*c0**4*cnye
            if(inew.eq.0) then !vardens
                pn=fn1(r,fnr)
            else
                pn=fn2(r,fnr,fnrr)
            end if
            dconst=vt0/(1.d-10*cnyfoc*pme*pn**2) !divided by 10^-10 here 
            !!!!!!!!                         and multiplied by 10^-10 in dfind()
            !!!old       dconst=vt0/(cnyfoc*pme*pn**2)
            !!!        dj(i)=dql(i,j)*dconst*vto !D_normir
            do i=1,ipt
                vrj(i)=vgrid(i,j)/vto      !Vpar/Vt
                dj(i)=dql(i,j)*dconst*vto  !D_normir
                !vrjnew(i,j,k)=vrj(i)
                !dijk(i,j,k)=dj(i)
            end do
            do i=1,i0
                if(vij(i,j).ge.vmax) then
                    ddout=zero
                else
                    call lock(vrj,ipt,vij(i,j),klo,khi,ierr)
                    if(ierr.eq.1) then
                        write(*,*)'lock error in output dql'
                        write(*,*)'j=',j,'ipt=',ipt
                        write(*,*)'vrj(1)=',vrj(1),' vrj(ipt)=',vrj(ipt)
                        write(*,*)'i=',i,' v=',vij(i,j),' vmax=',vmax
                        write(*,*)
                        pause'next key = stop'
                        stop
                    end if
                    !!!         call linf(vrj,dj,vij(i,j),ddout,klo,khi)
                    !!         if(ddout.le.1.d0) ddout=zero
                    ddout=dj(klo)
                end if
                dij(i,j,k)=ddout
            end do
            !zv1(j,k)=vrj(ipt1)
            !zv2(j,k)=vrj(ni1+ni2+ipt1)
        end do
    end 

    subroutine write_lhcd_power(time_stamp, ispectr) 
        use rt_parameters, only: nr
        use plasma, only: sk, vk
        use nr_grid, only: pdl, pdc        
        implicit none
        real(wp), intent(in) :: time_stamp
        integer,  intent(in) :: ispectr

        real(wp)      :: pwe
        integer       :: iu, i
        character(20) :: folder
        character(64) :: fname
    
        print *, 'write_lhcd_power time=', time_stamp
        !print *, name(m)
        if (ispectr>0) then
            folder = "lhcd/lhcd_power/pos/"
        else
            folder = "lhcd/lhcd_power/neg/"
        endif
    
        write(fname,'(A, f9.7,".dat")') folder, time_stamp
        print *, fname        

        open(newunit=iu, file=fname, status="replace", action="write")
	    write(iu,'(A5, 10A16)'), 'index', 'pwe', 'pdl', 'pdc', 'sk', 'vk'

        do i=1, nr
            pwe= (pdl(i)+pdc(i))/vk(i)
            write (iu, '(i6,5(ES22.14))') i, pwe, pdl(i), pdc(i), sk(i), vk(i)
        end do

        close(iu)

    end subroutine

    subroutine calculate_out_power(out_pe)
        use constants, only: zero, one
        use rt_parameters, only: nr, ismthout
        use plasma, only: rh, rh1, nspl, vk
        use nr_grid, only: pdl, pdc
        use math_module, only: fsmoth4
        use lock_module, only: lock2, linf
        implicit none
        real(wp), intent(inout) :: out_pe(*)
        real(wp) :: rxx(102)
        integer  :: i, j, nrr
        integer  :: klo,khi,ierr
        real(wp) :: hr, fout
        real(wp) :: pwe(102),wrk(102)

        hr = 1.d0/dble(nr+1)
        nrr=nr+2
        rxx(1)=zero
        rxx(nrr)=one
        do j=1,nr
            rxx(j+1)=hr*dble(j)
        end do

        do j=1,nr
            pwe(j+1)=(pdl(j)+pdc(j))/vk(j)
        end do
        pwe(1)=pwe(2)
        pwe(nr+2)=zero

        if(ismthout.ne.0) then
            do i=1,nrr
                wrk(i)=pwe(i)
            end do
            call fsmoth4(rxx,wrk,nrr,pwe)
        end if
        !
        rh(1)=rh1
        if(rh(nspl).gt.1.d0) rh(nspl)=1.d0
        do j=1,nspl
            call lock2(rxx,nrr,rh(j),klo,khi,ierr)
            if(ierr.ne.0) then
                write(*,*)'lock2 error in profiles for ASTRA'
                write(*,*)'ierr=',ierr,' j=',j,' rh(j)=',rh(j)
                write(*,*)'rxx(1)=',rxx(1),' rxx(nrr)=',rxx(nrr)
                pause
            end if
            call linf(rxx,pwe,rh(j),fout,klo,khi)
            out_pe(j)=fout
        end do
        rh(1)=zero
        !
    end 

    subroutine init_iteration_vars
        use constants, only : zero
        use power, only:  psum4, pnab, plost
        implicit none
        psum4=zero
        pnab=zero
        plost=zero
    end 

    subroutine init_alphas
        !! непонятная процедура. возможно результаты в dqi0 ???
        use constants, only: zero
        use rt_parameters, only: nr, itend0, kv
        use plasma, only:  vperp
        use nr_grid, only: dqi0, MAX_NR
        implicit none
        integer :: i, j
        real(wp) :: galfa(50,MAX_NR)  ! возможно массив должеб быть доступен еще где-то
        if(itend0.gt.0) then
            do j=1,nr           ! begin 'rho' cycle
                do i=1,50
                    dqi0(i,j)=zero
                end do
                call alphas(dqi0,vperp,j,kv,galfa)
            end do              ! end 'rho' cycle
        end if
        
    end subroutine

    subroutine alphas(d,u,j,kmax,g)
        use decrements, only: dgdu, kzero
        use constants, only : zero, one
        use nr_grid, only: MAX_NR
        implicit real*8 (a-h,o-z)
        integer, intent(in) :: j, kmax
        dimension d(50,MAX_NR),u(50,MAX_NR),g(50,MAX_NR)
        !common /arr/ dgdu(50,100),kzero(100)
        real(wp), parameter :: tiny=1.d-30
        integer :: k, km
        km=kzero(j)
        um=u(km,j)
        if(um.ge.one) then
            do k=1,kmax
                if(u(k,j).lt.one) then
                    uk=u(k,j)
                    uk2=uk**2
                    w=dsqrt(one-uk2)
                    g(k,j)=w/uk2
                    dgdu(k,j)=-one/(w*uk)-2.d0*w/(uk*uk2)
                else
                    g(k,j)=zero
                    dgdu(k,j)=zero
                end if
            end do
            return
        end if
  
        do k=1,km
            uk=u(k,j)
            uk2=uk**2
            w=dsqrt(one-uk2)
            g(k,j)=w/uk2
            dgdu(k,j)=-one/(w*uk)-2.d0*w/(uk*uk2)
        end do
  
        do k=km+1,kmax
            du=u(k,j)-u(k-1,j)
            if(u(k,j).lt.one) then
                beta=u(k,j)*dsqrt(one-u(k,j)**2)
            else
                beta=zero
            end if
            alfa=u(k,j)**3
            g(k,j)=(d(k,j)*g(k-1,j)+beta*du)/(d(k,j)+alfa*du)
            if(d(k,j).ne.zero) then
                dgdu(k,j)=(beta-alfa*g(k,j))/d(k,j)
            else
                dgdu(k,j)=(g(k,j)-g(k-1,j))/du
            end if
        end do
        do k=1,kmax
            if(g(k,j).lt.tiny) g(k,j)=zero
            if(dabs(dgdu(k,j)).lt.tiny) dgdu(k,j)=zero
        end do
        return
    end    
end module lhcd_module