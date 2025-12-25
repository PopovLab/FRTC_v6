module plasma
    !! модуль параметров плазмы
    use kind_module  
    use nr_grid, only: MAX_NR
    implicit none
    integer ngrid, nspl
    !! ASTRA radial grid number
    real(wp) tcur
    !! время (придумать название для переменной получше)
    real(wp) rm
    !! minor radius in mid-plane, cm
    real(wp) b_tor0
    !! тороидальное магнитное поле
    !! временно нужно две переменных, тоже нужно исправить
    real(wp) b_tor
    !! тороидальное магнитное поле

    real(wp) r0
    real(wp) z0
    real(wp) rh1
    real(wp), dimension(:),allocatable:: con,tem,temi,zeff,afld
    real(wp), dimension(:),allocatable:: rh,rha,drhodr,delta,ell,gamm,amy

    real(wp) tet1, tet2, gap_tet_plus, gap_tet_minus
    !! GRILL parameters 

    real(wp) xmi,cnye,cnyi,xsz,vt0 
    !!/a0ef3/ xmi,cnye,cnyi,xsz,vt0 
    real(wp) cnstvc

    real(wp) ww
    !! частота падающей волны 

    real(wp) cltn
    !!common /a0ef1/ cltn    

    real(wp) vperp(50,MAX_NR),cnstal,zza,zze,valfa!,kv
    !common /a0i5/ vperp(50,100),cnstal,zza,zze,valfa!,kv

    real(wp) vpmax

    real(wp) vk(MAX_NR), sk(MAX_NR)
    !common /a0i2/ vk(100)

    integer, parameter :: ipsy = 5, ncoef = 5
    !!   ipsy = number of polinomial decomposition coefficients
    !!   used for interpolation of Zakharov's moments.
    real(wp), dimension(ipsy) :: cdl,cly,cgm,cmy,coeffs


    real(wp) y2dn(501),y2tm(501),y2tmi(501)
    !!common /a0l3/
    real(wp) y2zeff(501)
    !!common /a0l5/ 

    integer ncheb
    real(wp) chebne(50),chebdne(50),chebddne(50)    
    !!common/ne_cheb

    real(wp) enorm(MAX_NR), fst(MAX_NR)
    !! em поле и еще что-то
    real(wp) dn1, dn2, dn3
contains
    subroutine init_plasma(NA1, ABC, BTOR, RTOR, UPDWN, GP2, AMETR, RHO, SHIF, ELON, TRIA,MU, NE, TE, TI, ZEF, UPL)
        use constants
        use approximation
        use rt_parameters
        use spline_module
        use chebyshev
        use math_module
        implicit none
        integer, intent(in)  :: NA1
        real(wp), intent(in) :: ABC, BTOR, RTOR, UPDWN, GP2
        real(wp), dimension(*) :: AMETR, RHO, SHIF, ELON, TRIA,MU,  NE, TE, TI, ZEF, UPL
        integer i, k
        integer, parameter :: N  = 501
        real(wp) :: znak_tor, znak_pol, fpol, dfmy

        ngrid = NA1
        nspl = ngrid
        if (.not. allocated(rh)) then
            allocate(rh(N),rha(N),drhodr(N),con(N),tem(N), source=0.0_wp)
            allocate(temi(N),zeff(N), afld(N), source=0.0_wp)
            allocate(delta(N),ell(N),gamm(N),amy(N), source=0.0_wp)
        end if
        do i=1, ngrid
            rh(i)=AMETR(i)/ABC
            rha(i)=RHO(i)/ABC  !/ABC instead of /ROC is not a mistake!
            delta(i)=(SHIF(1)-SHIF(i))/ABC  !FRTC Shafr. shift. defin.
            ell(i)=ELON(i)
            gamm(i)=rh(i)*TRIA(i)
            con(i)=NE(i)
            tem(i)=TE(i)
            temi(i)=TI(i)
            zeff(i)=ZEF(i)
            if (upl_fix) then
                afld(i)=upl_value/RTOR/GP2 
            else    
                afld(i)=UPL(i)/RTOR/GP2 
            endif
        end do
        rh(ngrid)=1.d0
        rh1=rh(1)          !saving the first ASTRA radial grid element
        rh(1) = 0.0d0         !shifting the first element to zero
        rha(1) = 0.0d0        !shifting the first element to zero
        delta(1) = 0.0d0      !putting delta(rh=0.)=0.
        gamm(1) = 0.0d0       !putting gamm(rh=0.)=0.
   
        b_tor0=1.d4*BTOR*RTOR/(RTOR+SHIF(1)) !B_tor_(magnetic axis), Gauss
        rm=1.d2*ABC                       !minor radius in mid-plane, cm
        r0=1.d2*(RTOR+SHIF(1))     !x-coordinate of the magnetic axis, cm
        z0=1.d2*UPDWN              !z-coordinate of the magnetic axis, cm


    !   spline approximation of plasma profiles         
    !
    !   shift as a function of "minor radius":
        call approx(rh,delta,ngrid,polin1,ipsy-1,coeffs)
        cdl(1)=0.0d0
        do k=2,ipsy
            cdl(k)=coeffs(k-1)
        end do
 
    !   triangularity as a function of "minor radius":
        call approx(rh,gamm,ngrid,polin1,ipsy-1,coeffs)
        cgm(1)=0.0d0
        do k=2,ipsy
            cgm(k)=coeffs(k-1)
        end do
 
    !   ellipticity as a function of "minor radius":
        call approx(rh,ell,ngrid,polin,ipsy,cly)            

    !  "poloidal magnetic field":
        call diff(rh,rha,ngrid,drhodr)
 
        do i=2,ngrid
            amy(i)=1.d4*BTOR*MU(i)*rha(i)*drhodr(i)
            !print *, amy(i), BTOR, MU(i)
        end do
        !print *, '----------------'
        amy(1)=0.d0  
        
    !! amy=(btor/q)*rho*(drho/dr) is a function of "minor radius" r=rh(i).
    !! Poloidal magnetic field: B_pol=amy(r)*sqrt(g22/g), where g is
    !! determinant of 3D metric tensor and g22 is the (22) element of
    !! the tensor, normalized on ABC^4 and ABC^2, correspondingly.
    !!
    !!  Polinomial approximation of the amy(r):
    !    inpt2=ngrid-3
        call approx(rh,amy,ngrid-3,polin1,ipsy-1,coeffs)
        cmy(1)=0.d0
        do k=2,ipsy
         cmy(k)=coeffs(k-1)
        end do
  
        ! зачем-то меняет знак коэффициентов????
        znak_tor=dsign(1.d0,dble(itor))
        b_tor=znak_tor*dabs(b_tor0)
        fpol=fdf(1.d0,cmy,ncoef,dfmy)
        znak_pol=dsign(1.d0,dble(i_pol))*dsign(1.d0,fpol)
        do i=1,ncoef
         cmy(i)=znak_pol*cmy(i)
        end do

        
    !!!!!!!!!!!!!!! spline approximation of plasma profiles !!!!!!!!!!!!!!!!
        call splne(rh,con,nspl,y2dn)
        call splne(rh,tem,nspl,y2tm)
        call splne(rh,zeff,nspl,y2zeff)
        call splne(rh,temi,nspl,y2tmi)

        if(inew.ne.0) then
            ncheb=20
            call chebft1(zero,1.d0,chebne,ncheb,fn)
            call chder(zero,1.d0,chebne,chebdne,ncheb)
            call chder(zero,1.d0,chebdne,chebddne,ncheb)
        end if    
        
        call init_parameters
        call find_volums_and_surfaces

    end subroutine

    function calc_theta(z, xly) result(tet)
        use constants
        implicit none
        real(wp) :: z, xly
        real(wp) :: arg, tet
        arg=(z-z0)/(xly*rm)
        if(dabs(arg).lt.1.d0) then
            tet=dasin(arg)      ! upper grill corner poloidal coordinate
        else
            tet=0.5d0*pi         ! upper grill corner poloidal coordinate
        end if
    end function

    subroutine init_parameters
        use constants
        use approximation
        use rt_parameters
        implicit none
        real(wp) :: xly, xlyp, arg1, arg2  
        real(wp) :: hr, sss
    !!!   
        xly = fdf(one,cly,ncoef,xlyp)
        tet1 = calc_theta(zplus, xly)
        tet2 = calc_theta(zminus, xly)
        gap_tet_plus = calc_theta(ZGapPlus, xly)
        gap_tet_minus = calc_theta(ZGapMinus, xly)
        
        !------------------------------------------------------------
        ! calculate constants
        !---------------------------------------
        hr = 1.d0/dble(nr+1)        
        dn1=1d0/(zi1+dni2*zi2+dni3*zi3)
        dn2=dni2*dn1
        dn3=dni3*dn1
        sss=zi1**2*dn1/xmi1+zi2**2*dn2/xmi2+zi3**2*dn3/xmi3
        xmi=1836.d0/sss
        cnstvc=(.75d0*piq*sss/1836.d0)**one_third
        ww=freq*pi2*1.0d+09 
        cnye=xlog/pi4
        cnyi=dsqrt(2d0)/(3d0*piq) !%for Vt=sqrt(Te/m)
        vt0=fvt(zero)
        !!!!!!!!      ptkev=ft(zero)/0.16d-8  !Te in keV
        cltn=clt/vt0
        xsz=clt/ww/rm
        !ccur=pqe*vt0*0.333d-9
        !!      ccurnr=pqe*pqe*0.333d-9/pme
        rrange=rrange*hr !ToDo если вызывается несколько раз то будут проблемы
        
        valfa=1.d9*dsqrt(1.91582d0*talfa/xmalfa)
        !  valfa (cgs units) = birth velocity
        zza=cnst1*(zalfa/xmalfa/valfa)**2*(clt/valfa)**3/pi
        zze=cnst2*2.d9*freq
        cnstal=(dsqrt(cnst1)/xmalfa/pi)*(zalfa*vt0/valfa)**2*clt/valfa
        vpmax=dsqrt(energy/talfa)
        !  "vpmax" in valfa velocity units !        
    end subroutine

    subroutine write_plasma(time_stamp)
        real(wp), intent(in) :: time_stamp
        character(120) fname
        integer, parameter :: iu = 20
        integer i
        write(fname,'("lhcd/plasma/", f9.7,".dat")')  time_stamp
        print *, fname
        open(iu, file=fname, status = 'new')
        write (iu, '(A)'), '#vars'
        write (iu, '(14A23)') 'b_tor0', 'rm', 'r0', 'z0'
        write (iu, ' (14(ES23.14))') b_tor0, rm, r0, z0
        write (iu, *) 

        write (iu, '(A)'), '#approx'
        write (iu, '(14A23)') 'cdl', 'cly', 'cgm', 'cmy'
        do i=1, ncoef
            write (iu, ' (14(ES23.14))') cdl(i),cly(i),cgm(i),cmy(i)
        end do  
        write (iu, *) 

        write (iu, '(A)'), '#radial_data'
        write (iu, '(14A23)') 'rh', 'rha', 'delta', 'ell', 'gamm', 'con', 'tem', 'temi', 'zeff', 'afld' 
        do i=1, ngrid
            write (iu, ' (14(ES23.14))') rh(i), rha(i), delta(i), ell(i), gamm(i), con(i), tem(i), temi(i), zeff(i), afld(i)
        end do                    

        close(iu)

    end subroutine

     subroutine write_equilibrium(time_stamp)
        real(wp), intent(in) :: time_stamp
        character(120) file_name
        integer :: unit_num, iostat
        integer i
        ! NAMELIST definition
        namelist /geometry/  rm, r0, z0
        namelist /fields/  b_tor0
        namelist /profile_approx/  cdl,cly,cgm,cmy
        namelist /radial_data/  rh, rha, delta, ell, gamm, con, tem, temi, zeff, afld
        
        write(file_name,'("lhcd/plasma/", f9.7,".nml")')  time_stamp
        print *, file_name
        open(newunit= unit_num, file= file_name, status= 'new')
        write(unit_num, nml=geometry, iostat=iostat)
        write(unit_num, nml=fields, iostat=iostat)
        write(unit_num, nml=profile_approx, iostat=iostat)
        write(unit_num, nml=radial_data, iostat=iostat)
        close(unit_num)
    end subroutine

    subroutine write_lcms
        !! write lcms
        use constants
        use approximation
        integer, parameter :: iu = 20
        integer  :: i
        real(wp) :: xr, th
        real(wp) :: xdl, xdlp, xly, xlyp, xgm, xgmp
        real(wp) :: x, xx, z, zz, pl, pc, pa
        real(wp) :: cotet, sitet
        open(iu,file='lhcd/out/lcms.dat')
        write(iu,*)'     R(m)            Z(m)'
        write(iu,*)
        xr=1.d0
        xdl=fdf(xr,cdl,ncoef,xdlp)
        xly=fdf(xr,cly,ncoef,xlyp)
        xgm=fdf(xr,cgm,ncoef,xgmp)
        do i=1,101
            th=dble(i-1)*pi2/dble(100)
            cotet=dcos(th)
            sitet=dsin(th)
            xx=-xdl+xr*cotet-xgm*sitet**2
            zz=xr*xly*sitet
            x=(r0+rm*xx)/1d2
            z=(z0+rm*zz)/1d2
            write(iu,'(6(e13.6,3x))') x,z
        end do
        close(iu)
    end subroutine

    subroutine find_volums_and_surfaces
        use constants
        use rt_parameters, only: nr
        use math_module, only: gaussint
        implicit none
        integer j
        real(wp) hr, rxx, vk0, sk0
        real(wp), parameter :: eps=1.d-6
        !--------------------------------------------------------
        ! find volums and surfaces
        !--------------------------------------------------------
        hr = 1.d0/dble(nr+1)
        vk0=pi2*hr*rm**3
        sk0=hr*rm**2
        do j=1,nr
            rxx=hr*dble(j)
            vk(j)=vk0*gaussint(obeom, zero, pi2, rxx, eps)
            sk(j)=sk0*gaussint(ploshad, zero, pi2, rxx, eps)
        end do        
    end subroutine

    real(wp) function fn(x)
    !! plasma  density,  cm^-3
        use constants, only: zero
        use spline_module      
        real(wp), intent(in) :: x
        real(wp) :: pa, r, y, dy
        real(wp), parameter :: alfa=4.d0, dr=.02d0
        pa=dabs(x)
        if(pa.le.rh(nspl)) then
            call splnt(rh,con,y2dn,nspl,pa,y,dy)
        else
            r=pa-rh(nspl)
            y=con(nspl)*dexp(-alfa*(r/dr)**2)
        end if
        fn=y*1.d+13    !cm^-3
    end    

    real(wp) function fvt(r)
    !! нет описания
        real(wp), intent(in) :: r
        real(wp) :: pt
        pt=ft(r)
        fvt=sqrt(pt/9.11d-28)
    end

    real(wp) function fn1(x,fnp)
    !! plasma density and its derivative
        use constants, only: zero
        use spline_module      
        real(wp), intent(in) :: x
        real(wp), intent(out) :: fnp
        real(wp) :: r, pa, y1, y, s, dy, dy1 
        real(wp), parameter :: alfa=4.d0, dr=.02d0
        pa=abs(x)
        if(pa.le.rh(nspl)) then
            call splnt(rh,con,y2dn,nspl,pa,y,dy)
        else
            call splnt(rh,con,y2dn,nspl,rh(nspl),y1,dy1)
            r=pa-rh(nspl)
            y=rh(nspl)*exp(-alfa*(r/dr)**2)
            dy=-2.d0*alfa*y*r/dr**2 !corrected
        end if
        fn1=y*1.d+13    !cm^-3
        fnp=dy*1.d+13
    end

    real(wp) function fn2(r, fnp, fnpp)
    !! plasma density and its first and second derivatives
        use constants, only: zero
        use chebyshev
        real(wp), intent(in) :: r
        real(wp), intent(out) :: fnp, fnpp
        real(wp) :: x, y1, y, s, dy, ddy 
        real(wp), parameter :: alfa=4.d0, dr=.02d0
        x=abs(r)
        if(x.le.1.d0) then
            y=chebev(zero,1.d0,chebne,ncheb,x)
            dy=chebev(zero,1.d0,chebdne,ncheb,x)
            ddy=chebev(zero,1.d0,chebddne,ncheb,x)
        else
            y1=chebev(zero,1.d0,chebne,ncheb,1.d0)
            s=x-1.d0
            y=y1*exp(-alfa*(s/dr)**2)
            dy=-2.d0*alfa*y*s/dr**2
            ddy=-2.d0*alfa*y*(1.d0-2.d0*alfa*(s/dr)**2)/dr**2
        end if
        fn2=y    !cm^-3
        fnp=dy
        fnpp=ddy
    end

    real(wp) function ft(x)
    !! electron temperature, erg
        use constants, only: zero
        use spline_module
        real(wp), intent(in) :: x
        real(wp) :: pa, r, y, dy
        real(wp), parameter :: alfa=4.d0, dr=.02d0
        pa=abs(x) !#@sav
        if(pa.le.rh(nspl)) then
            call splnt(rh,tem,y2tm,nspl,pa,y,dy)
        else
            r=pa-rh(nspl)
            y=tem(nspl)*exp(-alfa*(r/dr)**2)
        end if
        !!      ft=y            ! kev
        ft=y*0.16d-8      ! erg
    end    

    real(wp) function fti(x)
    !! ion temperature, kev
        use constants, only: zero
        use spline_module      
        real(wp), intent(in) :: x
        real(wp) :: pa, r, y, dy
        real(wp), parameter :: alfa=4.d0, dr=.02d0
        pa=abs(x) !#@sav
        if(pa.le.rh(nspl)) then
            call splnt(rh,temi,y2tmi,nspl,pa,y,dy)
        else
            r=pa-rh(nspl)
            y=temi(nspl)*exp(-alfa*(r/dr)**2)
        end if
        fti=y              ! kev
    end
    
    real(wp) function zefff(x)
    !! z_effective profile
        use constants, only: zero    
        use spline_module      
        real(wp), intent(in) :: x
        real(wp) :: pa, r, y, dy
        real(wp), parameter :: alfa=4.d0, dr=.02d0
        pa=abs(x) !#@sav
        if(pa.le.rh(nspl)) then
            call splnt(rh,zeff,y2zeff,nspl,pa,y,dy)
        else
            r=pa-rh(nspl)
            y=zeff(nspl)*exp(-alfa*(r/dr)**2)
        end if
        zefff=y
    end


    subroutine calc_enorm
        use constants
        use rt_parameters, only: nr, inew
        use spline_module
        use maxwell
        use lock_module
        implicit none
        integer j, klo,khi,ierr
        real(wp) :: efld
        real(wp) :: r, pn, vt, tmp, xlogj,vmax
        real(wp) :: fnr,fnrr, dens
        !real*8 fn1,fn2
        do j=1,nr
            r=dble(j)/dble(nr+1)
            call lock(rh,nspl,r,klo,khi,ierr)
            if(ierr.eq.1) then
                write(*,*)'lock error in saveprofiles, Efield'
                write(*,*)'j=',j,' rh(j)=',rh(j),' r=',r
                pause
                stop
            end if
            call linf(rh,afld,r,efld,klo,khi)
            if(inew.eq.0) then !vardens
                pn=fn1(r,fnr)
            else
                pn=fn2(r,fnr,fnrr)
            end if
            vt=fvt(r)
            tmp=ft(r)/0.16d-8  !Te,  KeV
            dens=pn/1.d+13     !10^13 cm^-3
            xlogj=dlog(5.1527d7*tmp*16.d0*dsqrt(tmp)/dsqrt(dens))
            enorm(j)=(3.835d0/xlogj)*efld*tmp/dens
            enorm(j)=enorm(j)*5.d0/(5.d0+zefff(r))
            !!fst(j)=pn*xlogj*c0**4/pi4/vt**3
            fst(j)=((5.d0+zefff(r))/5.d0)*pn*xlogj*c0**4/pi4/vt**3
        end do        
    end subroutine

    subroutine init_maxwell
        use constants
        use rt_parameters, only: nr, inew
        use spline_module
        use maxwell        
        use nr_grid, only: vij, fij, dfij, fij0, dij
        implicit none
        integer j
        real(wp) r, vclt
        do j=1,nr
            r=dble(j)/dble(nr+1)
            vclt=3.d10/fvt(r)
            !print *, vclt
            !call init_vi(vclt, vij(:,j))
            vij(1:i0,j) = create_vt_grid(vclt)
            call init_fmaxw_classic(vclt,enorm(j),fij(:,j,1),dfij(:,j,1)) ! positive
            call init_fmaxw_ext(vclt,enorm(j),fij(:,j,2),dfij(:,j,2))     ! negative
          end do
          fij0(:,:,:)=fij(:,:,:)
          dij(:,:,:)=zero

    end subroutine

    function obeom(ptet,pa) result(res)
        use constants
        use approximation
        implicit real*8 (a-h,o-z)
        real(wp), intent(in)    :: ptet
        real(wp), intent(in)    :: pa
        real(wp)                :: res
        !common /a0befr/ pi,pi2
        !common /a0ef1/ cltn
        !common /a0k/ cdl(10),cly(10),cgm(10),cmy(10),ncoef
        parameter(pa0=0.d0)
        xdl=fdf(pa,cdl,ncoef,xdlp)
        xly=fdf(pa,cly,ncoef,xlyp)
        xgm=fdf(pa,cgm,ncoef,xgmp)
        xlyv=xlyp*pa+xly
        cotet=dcos(ptet)
        sitet=dsin(ptet)
        dxdr=-xdlp+cotet-xgmp*sitet**2
        dxdt=-(pa+two*xgm*cotet)*sitet
        dzdr=xlyv*sitet
        dzdt=xly*pa*cotet
        x0=r0/rm-xdl+pa*cotet-xgm*sitet**2
        dxdrdt=-sitet-two*xgmp*sitet*cotet
        dzdrdt=xlyv*cotet
        dxdtdt=-pa*cotet-two*xgm*(cotet**2-sitet**2)
        dzdtdt=-xly*pa*sitet
        x0t=dxdt
!--------------------------------------
! components of metric tensor
!--------------------------------------
        g11=dxdr**2+dzdr**2
        g22=dxdt**2+dzdt**2
        g12=dxdr*dxdt+dzdr*dzdt
        g33=x0**2
        xj=(dzdr*dxdt-dxdr*dzdt)**2  !gg=g11*g22-g12*g12
        g=xj*g33
        res=dsqrt(g)
    end

    function ploshad(ptet,pa) result(res)
        use constants
        use approximation
        implicit real*8 (a-h,o-z)
        real(wp), intent(in)    :: ptet
        real(wp), intent(in)    :: pa
        real(wp)                :: res
        !common /a0befr/ pi,pi2
        !common /a0ef1/ cltn
        !common /a0k/ cdl(10),cly(10),cgm(10),cmy(10),ncoef
        parameter(pa0=0.d0)
        xdl=fdf(pa,cdl,ncoef,xdlp)
        xly=fdf(pa,cly,ncoef,xlyp)
        xgm=fdf(pa,cgm,ncoef,xgmp)
        xlyv=xlyp*pa+xly
        cotet=dcos(ptet)
        sitet=dsin(ptet)
        dxdr=-xdlp+cotet-xgmp*sitet**2
        dxdt=-(pa+two*xgm*cotet)*sitet
        dzdr=xlyv*sitet
        dzdt=xly*pa*cotet
        x0=r0/rm-xdl+pa*cotet-xgm*sitet**2
        dxdrdt=-sitet-two*xgmp*sitet*cotet
        dzdrdt=xlyv*cotet
        dxdtdt=-pa*cotet-two*xgm*(cotet**2-sitet**2)
        dzdtdt=-xly*pa*sitet
        x0t=dxdt
        !--------------------------------------
        ! components of metric tensor
        !--------------------------------------
        g11=dxdr**2+dzdr**2
        g22=dxdt**2+dzdt**2
        g12=dxdr*dxdt+dzdr*dzdt
        xj=(dzdr*dxdt-dxdr*dzdt)**2  !gg=g11*g22-g12*g12
        res=dsqrt(xj)
    end


end module plasma
