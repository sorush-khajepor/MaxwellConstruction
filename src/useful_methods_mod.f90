module useful_methods_mod
    use method_mod
    use eos_mod
    implicit none

!   Methods are the linkers between solvers and other classes
!   for more details refer to method_mod.f90

!   m4p_p0  = method  4 p-p0 or EOS-p0
    type,extends(method_t) :: m4p_p0_t
        class(eos_t),pointer :: eos
        real(kind=wp) :: p0
    contains
        procedure :: func => m4p_p0_func
        procedure :: init => init_m4p_p0
    end type

!   m4mi = method 4 maxwell integral
    type,extends(method_t) :: m4mi_t
        class(eos_t),pointer :: eos
        real(kind=wp) :: p0
    contains
        procedure :: func => m4mi_func
        procedure :: init => init_m4mi
    end type

!   m4dp_drho = mathod for dp/drho
    type,extends(method_t) :: m4dp_drho_t
        class(eos_t),pointer :: eos
    contains
        procedure :: init => init_m4dp_drho
        procedure :: func => m4dp_drho_func
    end type

!   m4pres = method for pressure
    type,extends(m4dp_drho_t) :: m4pres_t
    contains
        procedure :: func => m4pres_func
    end type

contains
subroutine init_m4p_p0(this,eos)
    class(m4p_p0_t) :: this
    class(eos_t),target   :: eos
    this%eos=>eos
end subroutine
function m4p_p0_func(this,x) result (f)
    class(m4p_p0_t) :: this
    real(kind=wp),intent(in):: x
    real(kind=wp) :: f
    f=this%eos%getp(x)-this%p0
end function

subroutine init_m4mi(this,eos)
    class(m4mi_t) :: this
    class(eos_t),target   :: eos
    this%eos=>eos
end subroutine

function m4mi_func(this,x) result (f)
    class(m4mi_t) :: this
    real(kind=wp),intent(in):: x
    real(kind=wp) :: f
    f=(this%eos%getp(x)-this%p0)/x**2
end function

subroutine init_m4dp_drho(this,eos)
    class(m4dp_drho_t) :: this
    class(eos_t),target   :: eos
    this%eos=>eos
end subroutine

function m4dp_drho_func(this,x) result (f)
    class(m4dp_drho_t) :: this
    real(kind=wp),intent(in):: x
    real(kind=wp) :: f
    f=this%eos%getdp_drho(x)
end function

function m4pres_func(this,x) result (f)
    class(m4pres_t) :: this
    real(kind=wp),intent(in):: x
    real(kind=wp) :: f
    f=this%eos%getp(x)
end function

end module
