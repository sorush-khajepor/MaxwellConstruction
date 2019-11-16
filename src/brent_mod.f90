!------------------------------------------------------------
! Using Brent's method, find the root of a function func known to lie between x1 and x2.
! The root, returned as root, will be refined until its accuracy is tol.
! Parameters: Maximum allowed number of iterations, and machine ?oating-point precision.
!-------------------------------------------------------------

module brent_mod
    use method_mod
    use kinds_f
    implicit none



!   yroot = f(root)
    type brent_t
        class(method_t), pointer :: method
        real(kind=wp) :: mbound_st,mbound_nd,tol,eps
        real(kind=wp),allocatable :: bounds(:,:),roots(:)
        integer :: itmax,res,numoroots,maxroots
    contains
        procedure :: init => init_brent
        procedure :: setmethod => setmethod_brent
        procedure :: getroot => get_Brent_root
        procedure :: f => curve_function
        procedure :: analyse_curve => analyse_curve
        procedure :: find_all_roots => find_all_roots
    end type
    private :: curve_function
contains
!---------------------------------------------------------------------


subroutine setmethod_brent(this,method)
    class(brent_t) :: this
    class(method_t),target :: method
    this%method => method
end subroutine
subroutine init_brent(this,method,mbound_st,mbound_nd,tol,itmax,eps,res,maxroots)
    class(brent_t):: this
    class(method_t),target :: method
    real(kind=wp), intent(in):: mbound_st,mbound_nd
    real(kind=wp), intent(in),optional:: eps,tol
    integer, intent(in),optional:: itmax,res,maxroots
    call this%setmethod(method)
    this%mbound_st = mbound_st
    this%mbound_nd = mbound_nd
    if (present(itmax)) then
        this%itmax = itmax
    else
        this%itmax = 1000
    end if
    if (present(eps)) then
        this%eps = eps
    else
        this%eps =   epsilon(1.0_wp)
    end if
    if (present(tol)) then
        this%tol = tol
    else
        this%tol =   epsilon(1.0_wp)
    end if
    if (present(res)) then
        this%res = res
    else
        this%res = 1000
    end if
    if (present(maxroots)) then
        this%maxroots = maxroots
    else
        this%maxroots = 50
    end if

end subroutine

function curve_function (this,x) result (y)
    class(brent_t) :: this
    real(kind=wp) :: x,y
    y = this%method%func(x)
end function

subroutine analyse_curve(this)
    class(brent_t):: this
    real(kind=wp) :: dx,x,tmp_bounds(1:2,1:this%maxroots)
    integer :: i,signchange
    associate (mbound_st=>this%mbound_st,mbound_nd=>this%mbound_nd,res=>this%res,nroot=>this%numoroots)
    dx = (mbound_nd-mbound_st)/res
    x = mbound_st
    signchange = 0
    do i = 1, res
        x=x+dx
        if (this%f(x)*this%f(x-dx).lt.0._wp) then
            signchange = signchange + 1
            tmp_bounds(1,signchange) = x-dx
            tmp_bounds(2,signchange) = x
        end if
    end do

    nroot = signchange
    if (signchange.eq.0) stop "Error in analyse_curve : no roots found!"

    if(allocated(this%bounds)) deallocate(this%bounds)
    if(allocated(this%roots))  deallocate(this%roots)
    allocate (this%bounds(1:2,1:nroot))
    allocate(this%roots(1:nroot))

    this%bounds(1:2,1:nroot) = tmp_bounds(1:2,1:nroot)

    end associate
end subroutine
! call this function when there are many roots in betweein
! main interval [mbound_st,mbound_nd]
subroutine find_all_roots(this)
    class(brent_t) :: this
    integer :: i
    associate(nroot=>this%numoroots)
    call this%analyse_curve()
    do i=1,nroot
        this%roots(i) = this%getroot(this%bounds(1,i),this%bounds(2,i))
    enddo
    end associate
end subroutine

! This function originally writtern by by Richard Brent and John Burkardt
! I slightly modify it to be object oriented. See the explanation below.

function get_Brent_root (this, a, b, is_root ) result(zero)
!*****************************************************************************80
!
!! ZERO seeks the root of a function F(X) in an interval [A,B].
!
!  Discussion:
!
!    The interval [A,B] must be a change of sign interval for F.
!    That is, F(A) and F(B) must be of opposite signs.  Then
!    assuming that F is continuous implies the existence of at least
!    one value C between A and B for which F(C) = 0.
!
!    The location of the zero is determined to within an accuracy
!    of 6 * MACHEPS * abs ( C ) + 2 * T.
!
!    Thanks to Thomas Secretin for pointing out a transcription error in the
!    setting of the value of P, 11 February 2013.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2013
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the change of 
!    sign interval.
!
!    Input, real ( kind = 8 ) MACHEP, an estimate for the relative machine
!    precision.
!
!    Input, real ( kind = 8 ) T, a positive error tolerance.
!
!    Input, external real ( kind = 8 ) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose zero is being sought.
!
!    Output, real ( kind = 8 ) ZERO, the estimated value of a zero of

!    the function F.
!
  implicit none

  class(brent_t) :: this
  real ( kind = wp ) c
  real ( kind = wp ) d
  real ( kind = wp ) e
  real ( kind = wp ) f
  real ( kind = wp ) fa
  real ( kind = wp ) fb
  real ( kind = wp ) fc
  real ( kind = wp ) m
  real ( kind = wp ) p
  real ( kind = wp ) q
  real ( kind = wp ) r
  real ( kind = wp ) s
  real ( kind = wp ) sa
  real ( kind = wp ) sb
  real ( kind = wp ) tol
  real ( kind = wp ) zero
  logical, intent(out), optional :: is_root
  real(kind=wp), intent(in),optional :: a, b

!
!  Make local copies of A and B.
!

    if (.not.(present(a).and.present(b))) then
        sa = this%mbound_st
        sb = this%mbound_nd
    else
        sa = a 
        sb = b  
    end if

  fa = this%f(sa)
  fb = this%f(sb)


        if (present(is_root)) is_root=.false.
        if((fa.gt.0.0_wp.and.fb.gt.0.0_wp).or.(fa.lt.0.0_wp.and.fb.lt.0.0_wp)) then
                return
        endif
        if (present(is_root)) is_root=.true.


  c = sa
  fc = fa
  e = sb - sa
  d = e

  do

    if ( abs ( fc ) < abs ( fb ) ) then

      sa = sb
      sb = c
      c = sa
      fa = fb
      fb = fc
      fc = fa

    end if

    tol = 2.0D+00 * this%eps * abs ( sb ) + this%tol 
    m = 0.5D+00 * ( c - sb )

    if ( abs ( m ) <= tol .or. fb == 0.0D+00 ) then
      exit
    end if

    if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then

      e = m
      d = e

    else

      s = fb / fa

      if ( sa == c ) then

        p = 2.0D+00 * m * s
        q = 1.0D+00 - s

      else

        q = fa / fc
        r = fb / fc
        p = s * ( 2.0D+00 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0D+00 ) )
        q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

      end if

      if ( 0.0D+00 < p ) then
        q = - q
      else
        p = - p
      end if

      s = e
      e = d

      if ( 2.0D+00 * p < 3.0D+00 * m * q - abs ( tol * q ) .and. &
        p < abs ( 0.5D+00 * s * q ) ) then
        d = p / q
      else
        e = m
        d = e
      end if

    end if

    sa = sb
    fa = fb

    if ( tol < abs ( d ) ) then
      sb = sb + d
    else if ( 0.0D+00 < m ) then
      sb = sb + tol
    else
      sb = sb - tol
    end if

    fb = this%f(sb) 

    if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or. &
         ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
      c = sa
      fc = fa
      e = sb - sa
      d = e
    end if

  end do

  zero = sb

  return
end


end module

