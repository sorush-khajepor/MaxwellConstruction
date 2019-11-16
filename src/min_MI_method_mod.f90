! Module for minimization of Maxwell Integral
module min_MI_method_mod
    use useful_methods_mod
    use brent_mod
    use quanc8_mod
    implicit none

!   m4minmi = method for minimization of Maxwell Integral
!   Chek method_mod for more details on methods :).
!   initialization of the class pointer with null() makes
!   us sure they are not associated and it is necessary.
!   Otherwise compiler think they are associated!
    type,extends(method_t) :: m4minmi_t
        class(brent_t),     pointer  :: brent_mml=>null(),brent_sping=>null(),brent_spinl=>null()
        class(quanc8_t),    pointer  :: quanc=>null()
        class(eos_t),       pointer  :: eos=>null()
        class(m4dp_drho_t), pointer  :: m4dp_drho=>null()
        class(m4p_p0_t),    pointer  :: m4p_p0=>null()
        class(m4mi_t),      pointer  :: m4mi=>null()
        real(kind=wp)   :: rhomin,rhomax,pmin,pmax
        real(kind=wp)   :: rho_st,rho_nd,tol
    contains
        procedure :: func => m4minmi_func
        procedure :: init => init_m4minmi
        procedure :: getrhog => getrhog
        procedure :: getrhol => getrhol
        procedure :: getintegral_err => getintegral_err
    end type
contains
subroutine init_m4minmi(this,eos,rho_st,rho_nd,tol)
    class(m4minmi_t)     :: this
    class(eos_t),target :: eos
    real(kind=wp) :: rho_st,rho_nd
    real(kind=wp),intent(in),optional :: tol
    associate(brent_mml=>this%brent_mml,brent_sping=>this%brent_sping,brent_spinl=>this%brent_spinl,&
              quanc=>this%quanc,m4dp_drho=>this%m4dp_drho,m4p_p0=>this%m4p_p0,m4mi => this%m4mi)

    this%eos=>eos
    this%rho_st = rho_st
    this%rho_nd = rho_nd
    if (present(tol)) then
        this%tol=tol
    else
        this%tol = epsilon(1._wp)
    end if
    if (.not.associated(this%brent_mml)) &
    allocate(this%brent_mml,this%brent_sping,this%brent_spinl,this%quanc,this%m4dp_drho,this%m4p_p0,this%m4mi)

    call m4dp_drho%init(eos)
    call m4p_p0%init(eos)
    call m4mi%init(eos)

    call brent_mml%init(method=m4dp_drho,mbound_st=rho_st,mbound_nd=rho_nd,tol=this%tol)
    call brent_mml%find_all_roots()

    this%rhomax = brent_mml%roots(1)
    this%rhomin = brent_mml%roots(2)
    this%pmin = eos%getp(this%rhomin)
    this%pmax = eos%getp(this%rhomax)


    call brent_sping%init (method=m4p_p0, mbound_st=this%rho_st,mbound_nd=this%rhomax,tol=this%tol)
    call brent_spinl%init (method=m4p_p0, mbound_st=this%rhomin,mbound_nd=rho_nd,tol=this%tol)
    call quanc%init(method=m4mi,mbound_st=this%rho_st,mbound_nd=this%rho_nd)

    end associate
end subroutine

function m4minmi_func(this,x) result (f)
    class(m4minmi_t) :: this
    real(kind=wp),intent(in):: x
    real(kind=wp) :: f
    real(kind=wp) :: rhog,rhol,p0

    p0 = x

    this%m4mi%p0 = p0

    rhog = this%getrhog(p0)
    rhol = this%getrhol(p0)
    f = this%quanc%integrate(sbound_st=rhog,sbound_nd=rhol)
end function
! finds the intersection of y=p0 and vapor branch of eos
function getrhog(this,p0) result(rho)
    class(m4minmi_t) :: this
    real(kind=wp) :: p0,rho
    logical :: is_root
    this%m4p_p0%p0 = p0
    rho = this%brent_sping%getroot(is_root=is_root)
    if (.not.is_root) stop "Error: no density root found in p-p0 for vapor branch."
end function

! finds the intersection of y=p0 and liquid branch of eos
function getrhol(this,p0) result(rho)
    class(m4minmi_t) :: this
    real(kind=wp) :: p0,rho
    logical :: is_root
    this%m4p_p0%p0 = p0
    rho = this%brent_spinl%getroot(is_root=is_root)
    if (.not.is_root) stop "Error: no density root found in p-p0 for liquid branch."
end function

function getintegral_err(this,rhog,rhol,p0) result(err)
    class(m4minmi_t) :: this
    real(kind=wp),intent(in) :: rhog,rhol,p0
    real(kind=wp) :: err

    this%m4mi%p0 = p0

    err = this%quanc%integrate(sbound_st=rhog,sbound_nd=rhol)

end function

end module













