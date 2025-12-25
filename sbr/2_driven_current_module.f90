module driven_current_module
    !! Driven Current Module
    use kind_module
    implicit none

    !real(wp)  :: zv1(100,2), zv2(100,2)
    !common/plosh/ zv1(100,2),zv2(100,2)!,sk(100)
    
    type DrivenCurrent
    real(wp) :: cu 
    !! ??  может лучше cuj

    real(wp) :: cu0    
    !! ??              cujoh

    real(wp) :: c
    !! ??

    real(wp) :: c0
    !! ??

    real(wp), dimension(:), allocatable :: current
    !! current(nr)

    real(wp), dimension(:), allocatable :: outj
    !! outj(i)  = LH driven current density, MA/m^2

    real(wp), dimension(:), allocatable :: ohj
    !! 
    integer   :: grid_size

    contains
        procedure :: evaluate => DrivenCurrent_evaluate
    end type

    interface DrivenCurrent
        module procedure :: DrivenCurrent_constructor
    end interface DrivenCurrent

    type DrivenCurrentResult
        real(wp) :: cup, cp
        !!
        real(wp) :: cum, cm
        !!
        real(wp) :: cup0, cp0
        !!
        real(wp) :: cum0, cm0
        !!

    contains
        procedure :: print => driven_current_result_print
        procedure :: save  => driven_current_result_save

    end type DrivenCurrentResult

    interface DrivenCurrentResult
        module procedure :: DrivenCurrentResult_constructor
    end interface DrivenCurrentResult    
    
contains

    function DrivenCurrent_constructor(size) result(this)
        !- конструктор для DrivenCurrent
        use rt_parameters, only : nr
        implicit none
        type(DrivenCurrent) :: this
        integer, value :: size

        this%cu  = 0
        this%cu0 = 0   
        this%c   = 0
        this%c0  = 0

        allocate(this%current(nr))

        this%grid_size = size
        allocate(this%outj(size), this%ohj(size))

    end function DrivenCurrent_constructor

    subroutine DrivenCurrent_evaluate(this, ROC)
        !- заключительное вычисление DrivenCurrent
        use constants, only: zero
        implicit none
        class(DrivenCurrent), intent(inout) :: this
        real(wp),             intent(in)    :: ROC
        external aiint
        real(wp) aiint
        integer   :: i
        if (this%cu0.ne.zero) then
            this%c0 = aiint(this%ohj, ROC)
            if (this%c0.ne.zero) then
                do i=1, this%grid_size
                    this%ohj(i)=this%cu0*this%ohj(i)/this%c0
                end do
            end if
        end if
        if (this%cu.ne.zero) then
            this%c = aiint(this%outj,ROC)
            if (this%c.ne.zero) then
                do i=1, this%grid_size
                    this%outj(i) = this%cu*this%outj(i)/this%c
                end do
            end if
        end if

    end subroutine DrivenCurrent_evaluate


    function DrivenCurrentResult_constructor(positive_dc, negative_dc) result(this)
        !- конструктор для DrivenCurrentResult
        implicit none
        type(DrivenCurrent), intent(in) :: positive_dc
        type(DrivenCurrent), intent(in) :: negative_dc
        type(DrivenCurrentResult) :: this

        this%cup  = positive_dc%cu
        this%cum  = negative_dc%cu
        this%cp   = positive_dc%c
        this%cm   = negative_dc%c
        this%cup0 = positive_dc%cu0
        this%cum0 = negative_dc%cu0
        this%cp0  = positive_dc%c0
        this%cm0  = negative_dc%c0

    end function DrivenCurrentResult_constructor

    subroutine driven_current_result_print(this, time)
        class(DrivenCurrentResult), intent(in) :: this
        real(wp), intent(in)  :: time
        print *, '------- driven current ---------'
        print *, 'time=',time
        print *, 'cup=',this%cup,' cp=',this%cp
        print *, 'cum=',this%cum,' cm=',this%cm
        print *, 'cup0=',this%cup0,' cp0=',this%cp0
        print *, 'cum0=',this%cum0,' cm0=',this%cm0
        print *, 'sigma driven current, MA=', this%cp0 + this%cm0
        print *, 'driven current, MA=', this%cup + this%cum
        print *, '--------------------------------'
    end subroutine driven_current_result_print   

    subroutine driven_current_result_save(this, time)
        class(DrivenCurrentResult), intent(in) :: this
        real(wp), intent(in)  :: time        
        logical, save :: first_time=.TRUE.
        character(132) FNAME
        integer :: io
        FNAME = "lhcd/dc_result.dat"

		if (first_time) then
            open(newunit=io, file=FNAME, status="replace", action="write")
			write(io,'(100A22)'), "Time", 'cup', 'cp', 'cum', 'cm', 'cup0', 'cp0', 'cum0', 'cm0'
            close(io)
			first_time = .FALSE.
        else
            open(newunit=io, file=FNAME, position="append")
            write(io,'(100f22.14)') time,  this%cup, this%cp, this%cum, this%cm, this%cup0, this%cp0, this%cum0, this%cm0
            close(io)
		end if

    end subroutine driven_current_result_save       
end module driven_current_module