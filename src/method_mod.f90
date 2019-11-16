module method_mod
    use kinds_f
    implicit none
!   method is an standard framework with known interface to link
!   complex classes to one dimensional math solvers (of y=f(x)).
!   In other words, all different math calculation like
!   root finding, integration, differentiation can be acheived
!   through this framework.
!   therefore, while solvers and other class being connected, they
!   can be compiled seperately ;).
    type, abstract:: method_t
    contains
        procedure(func),deferred :: func
    end type
    abstract interface
        function func(this,x)  result (f)
            import method_t
            import wp
            class(method_t) :: this
            real(kind=wp),intent(in) :: x
            real(kind=wp) :: f
        end function
    end interface

end module
