module eos_mod
    use kinds_f
    implicit none

    type, abstract :: eos_t
        real(kind=wp) :: a,b,c,Pc,T,Tr,rho_c
        real(kind=wp), allocatable :: Tc
        character(len=5) :: tname
    contains
        procedure :: init_eos => init_eos
        generic   :: init => init_eos
        procedure :: getp_f_rho => getp_f_rho_eos
        procedure :: getp_f_rho_Tr => getp_f_rho_Tr_eos
        generic   :: getp => getp_f_rho, getp_f_rho_Tr
        procedure :: setcritpoint => setcritpoint_eos
        procedure :: getdp_drho => getdp_drho_f_rho_eos
        procedure :: getd2p_drho2 => getd2p_drho2_f_rho_eos
        procedure :: setT => setTemp_eos
        procedure :: getmyname => getmyname_eos
    end type

    type,extends(eos_t) :: vw_t
    contains
        procedure :: init_vw => init_vw
        generic   :: init=> init_vw
        procedure :: setcritpoint => setcritpoint_vw
        procedure :: getp_f_rho_Tr => getp_f_rho_Tr_vw
        procedure :: getp_f_rho => getp_f_rho_vw
        procedure :: getdp_drho => getdp_drho_f_rho_vw
        procedure :: getd2p_drho2 => getd2p_drho2_f_rho_vw
    end type

    type,extends(eos_t) :: cs_t
    contains
        procedure :: init_cs => init_cs
        generic   :: init=> init_cs
        procedure :: setcritpoint => setcritpoint_cs
        procedure :: getp_f_rho => getp_f_rho_cs
        procedure :: getdp_drho => getdp_drho_f_rho_cs
        procedure :: getd2p_drho2 => getd2p_drho2_f_rho_cs
    end type

    type,extends(eos_t) :: pr_t
        real(kind=wp) :: omega,alpha
    contains
        procedure :: init_pr => init_pr
        generic   :: init=> init_pr
        procedure :: setcritpoint => setcritpoint_pr
        procedure :: getp_f_rho => getp_f_rho_pr
        procedure :: getdp_drho => getdp_drho_f_rho_pr
        procedure :: getd2p_drho2 => getd2p_drho2_f_rho_pr
        procedure :: setT => setTemp_pr
    end type

!   in R-K EOS we assumed R  universal gas constant (specific) is 1
!   As term RT and T^(1/2) apear seperately.
    type,extends(eos_t) :: rk_t
    contains
        procedure :: init_rk => init_rk
        generic   :: init=> init_rk
        procedure :: setcritpoint => setcritpoint_rk
        procedure :: getp_f_rho => getp_f_rho_rk
        procedure :: getdp_drho => getdp_drho_f_rho_rk
        procedure :: getd2p_drho2 => getd2p_drho2_f_rho_rk
    end type

    type,extends(eos_t) :: srk_t
        real(kind=wp) :: omega,alpha
    contains
        procedure :: init_srk => init_srk
        generic   :: init=> init_srk
        procedure :: setcritpoint => setcritpoint_srk
        procedure :: getp_f_rho => getp_f_rho_srk
        procedure :: getdp_drho => getdp_drho_f_rho_srk
        procedure :: getd2p_drho2 => getd2p_drho2_f_rho_srk
        procedure :: setT => setTemp_srk
    end type

contains
! --------------------- eos_t --------------------------


subroutine init_eos(this)
    class(eos_t) :: this
    stop "Error:init_eos function is an abstract function"
end subroutine

function getp_f_rho_eos(this,rho) result(p)
    class(eos_t) :: this
    real(kind=wp),intent(in) :: rho
    real(kind=wp) :: p
    stop "Error:getp_f_rho function is an abstract function"
end function
function getp_f_rho_Tr_eos(this,rho,Tr) result(p)
    class(eos_t) :: this
    real(kind=wp),intent(in) :: rho,Tr
    real(kind=wp) :: p
    stop "Error: getp_f_rho_Tr function is an abstract function"
end function
subroutine setcritpoint_eos(this)
    class(eos_t) :: this
    stop "Error: setcritpoint function is an abstract function"
end subroutine
function getdp_drho_f_rho_eos(this,rho) result(dp)
    class(eos_t) :: this
    real(kind=wp),intent(in) :: rho
    real(kind=wp) :: dp
    stop "Error: getdp_drho_f_rho_eos function is an abstract function"
end function
function getd2p_drho2_f_rho_eos(this,rho) result(d2p)
    class(eos_t) :: this
    real(kind=wp),intent(in) :: rho
    real(kind=wp) :: d2p
    stop "Error: getd2p_drho2_f_rho_eos function is an abstract function"
end function

function getmyname_eos(this) result(tname)
    class(eos_t) :: this
    character(len=5) :: tname
    tname = this%tname
end function

subroutine setTemp_eos(this,Tr)
    class(eos_t) :: this
    real(kind=wp),intent(in) :: Tr

    if(.not.allocated(this%Tc)) stop "Error in setTemp_eos: Tc is not allocated!"

    this%Tr=Tr
    This%T=Tr*this%Tc

end subroutine


!---------------------VW_t ---------------------------------------

subroutine init_vw(this,a,b,Tr)
    class(vw_t) :: this
    real(kind=wp),intent(in) :: a,b,Tr
    this%tname = "VW"
    this%a = a
    this%b= b
    call this%setcritpoint()
    call this%setT(Tr=Tr)
end subroutine

subroutine setcritpoint_vw(this)
    class(vw_t) :: this
    allocate(this%Tc)
    this%Pc = this%a/(27._wp*this%b**2)
    this%Tc = 8._wp*this%a/(27._wp*this%b)
    this%rho_c = 1._wp/(3._wp*this%b)
end subroutine
function getp_f_rho_Tr_vw(this,rho,Tr) result(p)
    class(vw_t) :: this
    real(kind=wp),intent(in) :: rho,Tr
    real(kind=wp) :: p,T
    if(allocated(this%Tc)) then
        T = Tr*this%Tc
    else
        stop "Tc in VW class is not set!"
    end if
    p = rho*T/(1._wp-this%b*rho) - this%a*rho*rho
end function
function getp_f_rho_vw(this,rho) result(p)
    class(vw_t) :: this
    real(kind=wp),intent(in) :: rho
    real(kind=wp) :: p
    p = rho*this%T/(1._wp-this%b*rho) - this%a*rho*rho
end function
function getdp_drho_f_rho_vw(this,rho) result(dp)
    class(vw_t) :: this
    real(kind=wp),intent(in) :: rho
    real(kind=wp) :: dp
    associate(b=>this%b,a=>this%a,T=>this%T)
    dp = T/(1._wp-b*rho)+rho*T*b/(1._wp-b*rho)**2-2._wp*a*rho
    end associate
end function

function getd2p_drho2_f_rho_vw(this,rho) result(d2p)
    class(vw_t) :: this
    real(kind=wp),intent(in) :: rho
    real(kind=wp) :: d2p
    associate(b=>this%b,a=>this%a,T=>this%T)
    d2p = 2._wp * T / (1._wp - b * rho) ** 2 * b + 2._wp * rho * T / (1._wp - b * rho) &
          ** 3 * b ** 2 - 2._wp * a
    end associate
end function

!-----------CS ------------------------------------------------

subroutine init_cs(this,a,b,Tr)
    class(cs_t) :: this
    real(kind=wp),intent(in) :: a,b,Tr
    this%tname = "CS"
    this%a = a
    this%b= b
    call this%setcritpoint()
    call this%setT(Tr=Tr)
end subroutine

subroutine setcritpoint_cs(this)
    class(cs_t) :: this
    real(kind=wp) :: ca,cb,abt,r
    allocate(this%Tc)
    r = 0.52177553676981578055344660_wp
    ca = 0.49638805772940987326931762_wp
    cb = 0.18729456694673304069903652_wp

    this%Tc = cb*this%a/(ca*this%b)
    this%pc = cb/this%b*this%Tc
    this%rho_c = r/this%b
end subroutine

function getp_f_rho_cs(this,rho) result(p)
    class(cs_t) :: this
    real(kind=wp),intent(in) :: rho
    real(kind=wp) :: p,BB,r
    associate(T=>this%T,a=>this%a)
    BB = this%b/4._wp
    r = BB*rho
    p = rho*T*(1._wp+r+r**2-r**3)/(1._wp-BB*rho)**3 - a*rho**2
    end associate
end function
function getdp_drho_f_rho_cs(this,rho) result(dp)
    class(cs_t) :: this
    real(kind=wp),intent(in) :: rho
    real(kind=wp) :: dp,r,BB
    associate(a=>this%a,T=>this%T)

    BB = this%b/4._wp
    r = BB*rho

    dp = T * (1._wp + r + r ** 2 - r ** 3) / (1._wp - r) ** 3 + r / BB * T * (BB &
      + 2._wp * BB * r - 3._wp * BB * r ** 2) / (1._wp - r) ** 3 + 3 * r * T * (1._wp + r &
      + r ** 2 - r ** 3) / (1._wp - r) ** 4 - 2._wp * a * r / BB
    end associate
end function

function getd2p_drho2_f_rho_cs(this,rho) result(d2p)
    class(cs_t) :: this
    real(kind=wp),intent(in) :: rho
    real(kind=wp) :: d2p,r,BB
    associate(a=>this%a,T=>this%T)

    BB = this%b/4._wp
    r = BB*rho

    d2p = 2._wp * T * (BB + 2._wp * BB * r - 3._wp * BB * r ** 2) / (1._wp - r) ** 3 + 6._wp    &
     * T * (1._wp + r + r ** 2 - r ** 3) / (1._wp - r) ** 4 * BB + r / BB * T * (               &
     2._wp * BB ** 2 - 6._wp * BB ** 2 * r) / (1._wp - r) ** 3 + 6._wp * r * T * (BB + 2 *      &
      BB * r - 3._wp * BB * r ** 2) / (1._wp - r) ** 4 + 12._wp * r * BB * T * (1._wp + r       &
     + r ** 2 - r ** 3) / (1._wp - r) ** 5 - 2._wp * a
    end associate
end function
!--------------------------------PPPPPPP----RRRRRRRRR---------------------------

subroutine init_pr(this,a,b,omega,Tr)
    class(pr_t) :: this
    real(kind=wp),intent(in) :: a,b,omega,Tr
    this%tname = "PR"
    this%a = a
    this%b= b
    this%omega=omega
    call this%setcritpoint()
    call this%setT(Tr=Tr)

end subroutine

subroutine setTemp_pr(this,Tr)
    class(pr_t) :: this
    real(kind=wp),intent(in) :: Tr
    associate(alpha=>this%alpha,omega=>this%omega)

    if(.not.allocated(this%Tc)) stop "Error in setTemp_eos: Tc is not allocated!"

    this%Tr=Tr
    This%T=Tr*this%Tc
    alpha = (1._wp+ (0.37464_wp + 1.54226_wp*omega-0.26992_wp*omega**2)*(1._wp-sqrt(Tr)))**2

    end associate
end subroutine

subroutine setcritpoint_pr(this)
    class(pr_t) :: this
    real(kind=wp) :: c1,c2,rc
    allocate(this%Tc)
    rc = 0.25307658654159946227082744_wp    ! rc=r_critical, r = rho*b
    c1 = 0.07779607390388845597184472_wp
    c2 = 0.45723552892138218938346024_wp

    this%Pc = this%a/this%b**2 * c1**2/c2
    this%Tc = c1*this%a/(c2*this%b)
    this%rho_c = rc/this%b
end subroutine

function getp_f_rho_pr(this,rho) result(p)
    class(pr_t) :: this
    real(kind=wp),intent(in) :: rho
    real(kind=wp) :: p
    associate(T=>this%T,a=>this%a,b=>this%b,alpha=>this%alpha)
    p = rho*T/(1._wp-b*rho) - alpha*a*rho**2/(1._wp+2._wp*b*rho-(b*rho)**2)
    end associate
end function

function getdp_drho_f_rho_pr(this,rho) result(dp)
    class(pr_t) :: this
    real(kind=wp),intent(in) :: rho
    real(kind=wp) :: dp
    associate(T=>this%T,a=>this%a,b=>this%b,alpha=>this%alpha)


    dp = T / (1._wp - b * rho) + rho * T / (1._wp - b * rho) ** 2 * b - 2._wp * alpha &
      * a * rho / (1._wp + 2._wp * b * rho - b ** 2 * rho ** 2) + alpha * a       &
     * rho ** 2 / (1._wp + 2._wp * b * rho - b ** 2 * rho ** 2) ** 2 * (2._wp * b -   &
      2._wp * b ** 2 * rho)


    end associate
end function

function getd2p_drho2_f_rho_pr(this,rho) result(d2p)
    class(pr_t) :: this
    real(kind=wp),intent(in) :: rho
    real(kind=wp) :: d2p
    associate(T=>this%T,a=>this%a,b=>this%b,alpha=>this%alpha)

    d2p = 2._wp * T / (1._wp - b * rho) ** 2 * b + 2._wp * rho * T / (1._wp - b * rho) &
      ** 3 * b ** 2 - 2._wp * alpha * a / (1._wp + 2._wp * b * rho - b ** 2 * rho **2) &
      + 4._wp * alpha * a * rho / (1._wp + 2._wp * b * rho - b ** 2 * rho ** 2) &
      ** 2 * (2._wp * b - 2._wp * b ** 2 * rho) - 2._wp * alpha * a * rho ** 2 / (1._wp &
      + 2._wp * b * rho - b ** 2 * rho ** 2) ** 3 * (2._wp * b - 2._wp * b ** 2 *      &
     rho) ** 2 - 2._wp * alpha * a * rho ** 2 / (1._wp + 2._wp * b * rho - b ** 2 * &
     rho ** 2) ** 2 * b ** 2

    end associate
end function

!--------------------------------RRRR----kkkkk---------------------------

subroutine init_rk(this,a,b,Tr)
    class(rk_t) :: this
    real(kind=wp),intent(in) :: a,b,Tr
    this%tname = "RK"
    this%a = a
    this%b= b
    call this%setcritpoint()
    call this%setT(Tr=Tr)

end subroutine


subroutine setcritpoint_rk(this)
    class(rk_t) :: this
    real(kind=wp) :: c1,c2,rc
    allocate(this%Tc)
    rc = 0.25992104989487316476721061_wp    ! rc=r_critical, r = rho*b
    c1 = 0.08664034996495772158907019_wp
    c2 = 0.42748023354034140439099065_wp

    this%Tc = (c1*this%a/(c2*this%b))**(2._wp/3._wp)
    this%Pc =  this%Tc*c1/this%b
    this%rho_c = rc/this%b

end subroutine

function getp_f_rho_rk(this,rho) result(p)
    class(rk_t) :: this
    real(kind=wp),intent(in) :: rho
    real(kind=wp) :: p
    associate(T=>this%T,a=>this%a,b=>this%b)
    p = rho*T/(1._wp-b*rho) - a*rho**2/(sqrt(T)*(1._wp+b*rho))
    end associate
end function

function getdp_drho_f_rho_rk(this,rho) result(dp)
    class(rk_t) :: this
    real(kind=wp),intent(in) :: rho
    real(kind=wp) :: dp
    associate(T=>this%T,a=>this%a,b=>this%b)


    dp = T / (1._wp - b * rho) + rho * T * b / (1._wp - b * rho) ** 2 - 2._wp * a &
     * rho * T ** (-0.5_wp) / (1._wp + b * rho) + a * rho ** 2 * b *  &
      T ** (-0.5_wp) / (1._wp + b * rho) ** 2


    end associate
end function

function getd2p_drho2_f_rho_rk(this,rho) result(d2p)
    class(rk_t) :: this
    real(kind=wp),intent(in) :: rho
    real(kind=wp) :: d2p
    associate(T=>this%T,a=>this%a,b=>this%b)

    d2p = 2._wp * T * b / (1._wp - b * rho) ** 2 + 2._wp * rho * T * b ** 2 / (1._wp -   &
      b * rho) ** 3 - 2._wp * a * T ** (-0.5_wp) / (1._wp + b * rho) + 4._wp  &
      * a * rho * b * T ** (-0.5_wp) / (1._wp + b * rho) ** 2 - 2._wp *   &
      a * rho ** 2 * b ** 2 * T ** (-0.5_wp) / (1._wp + b * rho) ** 3

    end associate
end function

!-------------------SSSSSSS RRRR KKKKK ------------------------

!--------------------------------RRRR----kkkkk---------------------------

subroutine init_srk(this,a,b,omega,Tr)
    class(srk_t) :: this
    real(kind=wp),intent(in) :: a,b,Tr,omega
    this%tname = "SRK"
    this%a = a
    this%b= b
    this%omega = omega
    call this%setcritpoint()
    call this%setT(Tr=Tr)

end subroutine

subroutine setTemp_srk(this,Tr)
    class(srk_t) :: this
    real(kind=wp),intent(in) :: Tr
    associate(alpha=>this%alpha,omega=>this%omega)

    if(.not.allocated(this%Tc)) stop "Error in setTemp_eos: Tc is not allocated!"

    this%Tr=Tr
    This%T=Tr*this%Tc
    alpha = (1._wp+ (0.480_wp + 1.574_wp*omega-0.176_wp*omega**2)*(1._wp-sqrt(Tr)))**2

    end associate
end subroutine


subroutine setcritpoint_srk(this)
    class(srk_t) :: this
    real(kind=wp) :: ca,cb,rc
    allocate(this%Tc)
    rc = 0.25992104989487316476721061_wp    ! rc=r_critical, r = rho*b
    ca = 0.42748023354034140439099065_wp
    cb = 0.08664034996495772158907019_wp


    this%Tc = cb*this%a/(ca*this%b)
    this%Pc =  this%Tc*cb/this%b
    this%rho_c = rc/this%b

end subroutine

function getp_f_rho_srk(this,rho) result(p)
    class(srk_t) :: this
    real(kind=wp),intent(in) :: rho
    real(kind=wp) :: p
    associate(T=>this%T,a=>this%a,b=>this%b,alpha=>this%alpha)
    p = rho*T/(1._wp-b*rho) - alpha*a*rho**2/(1._wp+b*rho)
    end associate
end function

function getdp_drho_f_rho_srk(this,rho) result(dp)
    class(srk_t) :: this
    real(kind=wp),intent(in) :: rho
    real(kind=wp) :: dp
    associate(T=>this%T,a=>this%a,b=>this%b,alpha=>this%alpha)


    dp =  T / (1 - b * rho) + rho * T / (1 - b * rho) ** 2 * b - 2 * alpha &
      * a * rho / (1 + b * rho) + alpha * a * rho ** 2 / (1 + b * rho  &
     ) ** 2 * b


    end associate
end function

function getd2p_drho2_f_rho_srk(this,rho) result(d2p)
    class(srk_t) :: this
    real(kind=wp),intent(in) :: rho
    real(kind=wp) :: d2p
    associate(T=>this%T,a=>this%a,b=>this%b,alpha=>this%alpha)

    d2p = 2 * T / (1 - b * rho) ** 2 * b + 2 * rho * T / (1 - b * rho) ** 3 &
         * b ** 2 - 2 * alpha * a / (1 + b * rho) + 4 * alpha * a * rho &
         / (1 + b * rho) ** 2 * b - 2 * alpha * a * rho ** 2 / (1 + b *  &
         rho) ** 3 * b ** 2

    end associate
end function

end module
