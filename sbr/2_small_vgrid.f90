module small_vgrid
    use kind_module   
    use nr_grid, only: MAX_NR
    implicit none
    integer, parameter :: MAX_PT= 1000 ! размер температурной сетки
    real(wp) :: vmid(MAX_PT),vz1(MAX_PT),vz2(MAX_PT)
    !integer  :: ibeg(100),iend(100)

    real(wp) :: vrj(MAX_PT),dj(MAX_PT),djnew(MAX_PT)
    real(wp) :: dj2(MAX_PT),d2j(MAX_PT)

    real(wp) :: vgrid(MAX_PT,MAX_NR), dfundv(MAX_PT,MAX_NR)
    !!common/gridv/vgrid(101,100),dfundv(101,100)
    integer  :: nvpt
    !!common/gridv/nvpt
    integer :: ipt1, ipt2, ipt

    real(wp) :: vlf,vrt,dflf,dfrt
    !common /a0ghp/ vlf,vrt,dflf,dfrt
        
    integer, parameter :: kpt1=20, kpt3=20

contains

    subroutine distr(vz,j,ifound,fder)
        use lock_module      
        use rt_parameters, only: nr
        implicit none
        integer, intent(in) :: j
        integer, intent(inout) :: ifound
        real*8 vz,fder
        integer i,klo,khi,ierr,nvp
        real*8,dimension(:),allocatable:: vzj,dfdvj
        real(wp) :: dfout
        !real*8 vlf,vrt,dflf,dfrt
        !common /a0ghp/ vlf,vrt,dflf,dfrt
        !common/gridv/vgrid(101,100),dfundv(101,100),nvpt
        nvp=nvpt
        allocate(vzj(nvp),dfdvj(nvp))
        do i=1, nvp
            vzj(i)=vgrid(i,j)
            dfdvj(i)=dfundv(i,j)
        end do
        call lock2(vzj,nvp,vz,klo,khi,ierr)
        if(ierr.eq.0) then !vgrid(1,j) <= vz <= vgrid(nvpt,j)
            call linf(vzj,dfdvj,vz,dfout,klo,khi)
            ifound=klo
            vlf=vzj(klo)
            vrt=vzj(khi)
            fder=dfout
            dflf=dfdvj(klo)
            dfrt=dfdvj(khi)
        else if(ierr.eq.1) then !vz < vgrid(1,j)
            write(*,*)'exception: ierr=1 in distr()'
            pause'next key = stop'
            stop
        else if(ierr.eq.2) then !vz > vgrid(nvpt,j)
            !write(*,*)'exception: ierr=2 in distr()'
            print *, 'vz > vgrid(nvpt,j)'
            ifound=nvp
            vlf=vzj(nvp)
            vrt=vz
            fder=0
            dflf=0
            dfrt=0
            !pause'next key = stop'
            !stop
        else if(ierr.eq.3) then
            write(*,*)'exception in distr, klo=khi=',klo,' j=',j,' nvp=',nvp
            write(*,*)'vz=',vz,' v1=',vzj(1),' v2=',vzj(nvp)
            pause'next key = stop'
            stop
        end if
        deallocate(vzj,dfdvj)
    end   


    subroutine calculate_dfundv(ispectr)
        !! calculate dfundv что такое dfundv?
        use constants, only: zero
        use rt_parameters, only: nr
        use plasma, only: fvt, vt0, cltn
        use maxwell, only: i0
        use nr_grid, only: vij, dfij, dij
        use lock_module, only: lock, linf
        implicit none
        integer, intent(in) :: ispectr
        real(wp), dimension(:), allocatable:: vvj, vdfj
        integer  :: i, j, k 
        integer  :: klo,khi,ierr
        real(wp) :: dfout
        real(wp) :: r, hr
        real(wp) :: vt, vto
        hr = 1.d0/dble(nr+1)
        allocate(vvj(i0),vdfj(i0))
        k=(3-ispectr)/2
        do j=1,nr
            r=hr*dble(j)
            vt=fvt(r)
            vto=vt/vt0
            do i=1,i0
                vvj(i)=vij(i,j)
                vdfj(i)=dfij(i,j,k) !=dfundv(i,j)*vto**2
            end do
            do i=1,ipt
                vrj(i)=vgrid(i,j)/vto   !Vpar/Vt
                call lock(vvj,i0,vrj(i),klo,khi,ierr)
                if(ierr.eq.1) then
                    write(*,*)'lock error in read distribution function'
                    write(*,*)'j=',j,'i0=',i0
                    write(*,*)'vvj(1)=',vvj(1),' vvj(i0)=',vvj(i0)
                    write(*,*)'i=',i,' vrj(i)=',vrj(i),' vmax=',cltn/vto
                    write(*,*)
                    pause'next key = stop'
                    stop
                end if
                call linf(vvj,vdfj,vrj(i),dfout,klo,khi)
                dfundv(i,j)=dfout/vto**2
                if(dfundv(i,j).gt.zero) dfundv(i,j)=zero
            end do
        end do
        deallocate(vvj,vdfj)
    end

    subroutine gridvel(v1,v2,vmax,cdel,ni1,ni2,ipt1,kpt3,vrj)
        implicit none
        real(wp), intent(in) :: v1, v2, vmax, cdel
        integer,  intent(in) :: ni1, ni2, ipt1, kpt3
        real(wp), intent(inout) :: vrj(:)
        integer kpt1, kpt2, k
        real(wp) :: v12
        kpt1=ipt1-1
        kpt2=ni1+ni2+1
        do k=1,kpt1  !0<=v<v1
            vrj(k)=dble(k-1)*v1/dble(kpt1)
        end do
        v12=v1+(v2-v1)*cdel
        do k=1,ni1+1 !v1<=v<=v12
            vrj(k+kpt1)=v1+dble(k-1)*(v12-v1)/dble(ni1)
        end do
        do k=2,ni2+1 !!v12<v<=v2
            vrj(k+kpt1+ni1)=v12+dble(k-1)*(v2-v12)/dble(ni2)
        end do     
        do k=1,kpt3  !v2<v<=vmax
            vrj(k+kpt1+kpt2)=v2+dble(k)*(vmax-v2)/dble(kpt3)
        end do
    end    

    
    subroutine recalculate_f_for_a_new_mesh(ispectr, iterat)
        !!   recalculate f' for a new mesh
        use constants, only : zero
        use rt_parameters, only : nr, ni1, ni2, cdel
        use plasma, only: vt0, fvt, cltn
        use nr_grid, only: vzmin, vzmax
        use maxwell, only: i0
        use nr_grid, only: vij, dfij
        use lock_module        
        implicit none
        integer, intent(in) :: ispectr, iterat
        
        integer i, j, k
        real(wp) ::  dfout
        real(wp), dimension(:), allocatable:: vvj, vdfj
        integer :: klo,khi,ierr
        real(wp) :: r, hr, vt, vto, vmax
        real(wp) :: v1, v2, vp1, vp2

        allocate(vvj(i0),vdfj(i0))
        hr = 1.d0/dble(nr+1)
        k=(3-ispectr)/2
        do j=1,nr
            r=hr*dble(j)
            vt=fvt(r)
            vto=vt/vt0
            if(iterat.gt.0) then
                v1=dmin1(vzmin(j),vz1(j))
                v2=dmax1(vzmax(j),vz2(j))
            else
                v1=vzmin(j)
                v2=vzmax(j)
            end if
            vmax=cltn/vto
            vp1=v1/vto
            vp2=v2/vto
            call gridvel(vp1,vp2,vmax,cdel,ni1,ni2,ipt1,kpt3,vrj)
            do i=1,i0
                vvj(i)=vij(i,j)
                vdfj(i)=dfij(i,j,k) !=dfundv(i,j)*vto**2
            end do
            do i=1,ipt
                call lock(vvj,i0,vrj(i),klo,khi,ierr)
                if(ierr.eq.1) then
            !!!         if(vrj(i).gt.vvj(i0)) exit
                    write(*,*)'lock error in new v-mesh'
                    write(*,*)'j=',j,' i0=',i0
                    write(*,*)'vvj(1)=',vvj(1),' vvj(i0)=',vvj(i0)
                    write(*,*)'i=',i,' vrj(i)=',vrj(i)
                    write(*,*)
                    pause'next key = stop'
                    stop
                end if
                call linf(vvj,vdfj,vrj(i),dfout,klo,khi)
                vgrid(i,j)=vrj(i)*vto
                dfundv(i,j)=dfout/vto**2
                if(dfundv(i,j).gt.zero) dfundv(i,j)=zero
            end do
            vz1(j)=v1
            vz2(j)=v2
        end do
        deallocate(vvj,vdfj)
    end subroutine    

    subroutine find_velocity_limits_and_initial_dfdv(anb)
        use constants, only: c0, c1, zero, zalfa, xmalfa, xlog, one_third
        use rt_parameters, only: nr, inew, ni1, ni2, itend0, kv, factor
        use plasma !, only: fn1, fn2, fvt, vt0
        use nr_grid, only: dens, eta, fcoll
        use nr_grid, only: source !! нужна ли она тут вообще???
        use nr_grid, only: MAX_NR
        implicit none
        real(wp), intent(inout) :: anb
        integer  :: i, j, k
        real(wp) :: v, vt, vto, wpq, whe
        real(wp) :: u, u1, e1, e2, e3, tmp
        real(wp) :: cn1, cn2
        real(wp) :: tt, vmax, v1, v2
        real(wp) :: pn, fnr, fnrr
        real(wp) :: r, hr
        real(wp) :: dvperp, ddens, tdens
        real(wp) :: vpmin(MAX_NR), vcva(MAX_NR)
        hr = 1.d0/dble(nr+1)
        !c-------------------------------------------
        !c find velocity limits and initial dfdv
        !c--------------------------------------------
        ipt1=kpt1+1
        ipt2=ni1+ni2
        ipt=ipt1+ni1+ni2+kpt3
        if(ipt.gt.MAX_PT) then
            write(*,*)'ipt >MAX_PT'
            pause'stop program'
            stop
        end if
        nvpt=ipt

        do j=1,nr                  ! begin 'rho' cycle
            r=hr*dble(j)
            !!!!sav2008       pn=fn(r)
            !!       pn=fn1(r,fnr)
            !!       pn=fn2(r,fnr,fnrr) !sav2008
            if(inew.eq.0) then !vardens
                pn=fn1(r,fnr)
            else
                pn=fn2(r,fnr,fnrr)
            end if
            dens(j)=pn
            vt=fvt(r)
            vto=vt/vt0
            wpq=c0**2*pn
            whe=dabs(b_tor)*c1
            v=wpq/ww**2
            u1=whe/ww
            u=u1**2
            e1=1d0-v*(1d0/xmi-1d0/u)
            e2=v/u1
            e3=v
            tmp=ft(r)/0.16d-8 !Te, keV
            cn1=dsqrt(50d0/tmp)  !sav2008
            if(itend0.gt.0) then
                eta(j)=1d0-v
                vcva(j)=cnstvc*vt*dsqrt(2d0)/valfa
                vpmin(j)=2.0d0*dsqrt(tmp/(-eta(j)))
222             continue
                dvperp=(vpmax-vpmin(j))/dble(kv-1)
                if(dvperp.le.zero) then
                    vpmax=1.3d0*vpmax
                    go to 222
                end if
                do k=1,kv
                    vperp(k,j)=vpmin(j)+dble(k-1)*dvperp
                end do
                fcoll(j)=.5d-13*dens(j)*zalfa**2*xlog/xmalfa/tmp**1.5d0
                ddens=dn1*dens(j)
                tdens=dn2*dens(j)
                tt=fti(r)**one_third    ! (ti, keV)^1/3
                source(j)=4d-12*factor*ddens*tdens*dexp(-20d0/tt)/tt**2
                anb=anb+source(j)*vk(j)
            end if
            cn2=dsqrt(dabs(e1))+e2/dsqrt(e3) !sav2008
            !vz1(j)=cleft*cltn/cn1  !Vpar/Vt0
            !vz2(j)=cright*cltn/cn2  !Vpar/Vt0
            !if(vz2(j).gt.0.9d0*cltn) vz2(j)=0.9d0*cltn
            !v1=vz1(j)/vto !Vpar/Vt(rho)
            !v2=vz2(j)/vto !Vpar/Vt(rho)
            vmax=cltn/vto
            v1=4.d0  !Vpar/Vt(rho)
            v2=10.d0 !cright*cltn/cn2 !10.d0 !Vpar/Vt(rho)
            if(v2.ge.vmax) v2=0.5d0*vmax
            if(v1.ge.v2) v1=v2-2.d0
            call gridvel(v1,v2,vmax,0.5d0,ni1,ni2,ipt1,kpt3,vrj)
            vz1(j)=v1*vto !Vpar/Vt0
            vz2(j)=v2*vto !Vpar/Vt0
            if(vz2(j).gt.0.9d0*cltn) vz2(j)=0.9d0*cltn
            do i=1,ipt
                vgrid(i,j)=vrj(i)*vto
            end do
        end do                     ! end 'rho' cycle 


    end subroutine
end module small_vgrid