! Module for numerical integration
! Method is the function will be integrated.
module quanc8_mod
    use kinds_f
    use brent_mod, only : method_t
    implicit none

    type quanc8_t
        real(kind=wp) :: mbound_st,mbound_nd
        real(kind=wp) :: abserr,relerr,errest,flag
        integer :: nofun
        class(method_t), pointer :: method
    contains
        procedure :: init => init_quanc8
        procedure :: integrate => quanc8
        procedure :: f => curve_function
    end type

contains

subroutine init_quanc8(this,method,mbound_st,mbound_nd,abserr,relerr)
    class(quanc8_t) :: this
    class(method_t),target :: method
    real(kind=wp),intent(in) :: mbound_st,mbound_nd
    real(kind=wp),intent(in),optional:: abserr,relerr
    this%method => method
    this%mbound_st=mbound_st
    this%mbound_nd=mbound_nd
    if (present(abserr)) then
        this%abserr = abserr
    else
        this%abserr = 0._wp
    end if
    if (present(relerr)) then
        this%relerr = relerr
    else
        this%relerr=epsilon(1.0_wp)
    end if

end subroutine

function curve_function (this,x) result (y)
    class(quanc8_t) :: this
    real(kind=wp) :: x,y

    y = this%method%func(x)

end function

! Based on quanc8.f file from http://www.netlib.no/
! mbound = main bound boundary is set in quanc8_t object.
! but integration can be taken on smaller intervals as:
! sbound = short (local) boundary for this function
! if sbound is not set, mbound is used.
function quanc8(this,sbound_st,sbound_nd) result(ans)
      class (quanc8_t) :: this
      real(kind=wp)   :: ans
      real(kind=wp),optional::   sbound_st,sbound_nd
      real(kind=wp)  :: a,b


      !abserr=0.0, relerr=1.0e-16)
!
!    estimate the integral of this%f(x) from a to b
!    to a user provided tolerance.
!    an automatic adaptive routine based on
!    the 8-panel newton-cotes rule.
!
!    input ..
!
!    fun     the name of the integrand function subprogram this%f(x).
!    a       the lower limit of integration.
!    b       the upper limit of integration.(b may be less than a.)
!    relerr  a relative error tolerance. (should be non-negative)
!    abserr  an absolute error tolerance. (should be non-negative)
!
!    output ..
!
!    ans  an approximation to the integral hopefully satisfying the
!            least stringent of the two error tolerances.
!    errest  an estimate of the magnitude of the actual error.
!    nofun   the number of function values used in calculation of result.
!    flag    a reliability indicator.  if flag is zero, then result
!            probably satisfies the error tolerance.  if flag is
!            xxx.yyy , then  xxx = the number of intervals which have
!            not converged and  0.yyy = the fraction of the interval
!            left to do when the limit on  nofun  was approached.
!
      double precision w0,w1,w2,w3,w4,area,x0,f0,stone,step,cor11,temp
      double precision qprev,qnow,qdiff,qleft,esterr,tolerr
      double precision qright(31),fff(16),x(16),fsave(8,30),xsave(8,30)
      integer levmin,levmax,levout,nomax,nofin,lev,nim,i,j
      associate(abserr=>this%abserr,relerr=>this%relerr,errest=>this%errest,flag=>this%flag,nofun=>this%nofun)
      if (.not.(present(sbound_st).and.present(sbound_nd))) then
        a = this%mbound_st
        b = this%mbound_nd
      else
        a = sbound_st
        b = sbound_nd
      end if

!
!    ***   stage 1 ***   general initialization
!    set constants.
      levmin = 1
      levmax = 30
      levout = 6
      nomax = 5000
      nofin = nomax - 8*(levmax-levout+2**(levout+1))
!
!    trouble when nofun reaches nofin
!
      w0 =   3956.0_wp / 14175.0_wp
      w1 =  23552.0_wp / 14175.0_wp
      w2 =  -3712.0_wp / 14175.0_wp
      w3 =  41984.0_wp / 14175.0_wp
      w4 = -18160.0_wp / 14175.0_wp
!
!    initialize running sums to zero.
!
      flag = 0.0_wp
      ans = 0.0_wp
      cor11  = 0.0_wp
      errest = 0.0_wp
      area   = 0.0_wp
      nofun = 0
      if (a .eq. b) return
!
!    ***   stage 2 ***   initialization for first interval
!
      lev = 0
      nim = 1
      x0 = a
      x(16) = b

      qprev  = 0.0_wp
      f0 = this%f(x0)
      stone = (b - a) / 16.0_wp
      x(8)  =  (x0  + x(16)) / 2.0_wp
      x(4)  =  (x0  + x(8))  / 2.0_wp
      x(12) =  (x(8)  + x(16)) / 2.0_wp
      x(2)  =  (x0  + x(4))  / 2.0_wp
      x(6)  =  (x(4)  + x(8))  / 2.0_wp
      x(10) =  (x(8)  + x(12)) / 2.0_wp
      x(14) =  (x(12) + x(16)) / 2.0_wp
      do 25 j = 2, 16, 2
         fff(j) = this%f(x(j))
   25 continue
      nofun = 9
!
!    ***   stage 3 ***   central calculation
!    requires qprev,x0,x2,x4,...,x16,f0,f2,f4,...,f16.
!    calculates x1,x3,...x15, f1,f3,...f15,qleft,qright,qnow,qdiff,area.
!
   30 x(1) = (x0 + x(2)) / 2.0_wp
      fff(1) = this%f(x(1))
      do 35 j = 3, 15, 2
         x(j) = (x(j-1) + x(j+1)) / 2.0_wp
         fff(j) = this%f(x(j))
   35 continue
      nofun = nofun + 8
      step = (x(16) - x0) / 16.0_wp
      qleft  =  (w0*(f0 + fff(8))  + w1*(fff(1)+fff(7))  + w2*(fff(2)+fff(6)) &
       + w3*(fff(3)+fff(5))  +  w4*fff(4)) * step
      qright(lev+1)=(w0*(fff(8)+fff(16))+w1*(fff(9)+fff(15))+w2*(fff(10)+fff(14)) &
       + w3*(fff(11)+fff(13)) + w4*fff(12)) * step
      qnow = qleft + qright(lev+1)
      qdiff = qnow - qprev
      area = area + qdiff
!
!    ***   stage 4 *** interval convergence test
!
      esterr = abs(qdiff) / 1023.0_wp
      tolerr = max(abserr,relerr*abs(area)) * (step/stone)
      if (lev .lt. levmin) go to 50
      if (lev .ge. levmax) go to 62
      if (nofun .gt. nofin) go to 60
      if (esterr .le. tolerr) go to 70
!
!    ***   stage 5   ***   no convergence
!    locate next interval.
!
   50 nim = 2*nim
      lev = lev+1
!
!    store right hand elements for future use.
!
      do 52 i = 1, 8
         fsave(i,lev) = fff(i+8)
         xsave(i,lev) = x(i+8)
   52 continue
!
!    assemble left hand elements for immediate use.
!
      qprev = qleft
      do 55 i = 1, 8
         j = -i
         fff(2*j+18) = fff(j+9)
         x(2*j+18) = x(j+9)
   55 continue
      go to 30
!
!    ***   stage 6   ***   trouble section
!    number of function values is about to exceed limit.
!
   60 nofin = 2*nofin
      levmax = levout
      flag = flag + (b - x0) / (b - a)
      go to 70
!
!    current level is levmax.
!
   62 flag = flag + 1.0_wp
!
!    ***   stage 7   ***   interval converged
!    add contributions into running sums.
!
   70 ans = ans + qnow
      errest = errest + esterr
      cor11  = cor11  + qdiff / 1023.0_wp
!
!    locate next interval.
!
   72 if (nim .eq. 2*(nim/2)) go to 75
      nim = nim/2
      lev = lev-1
      go to 72
   75 nim = nim + 1
      if (lev .le. 0) go to 80
!
!    assemble elements required for the next interval.
!
      qprev = qright(lev)
      x0 = x(16)
      f0 = fff(16)
      do 78 i = 1, 8
         fff(2*i) = fsave(i,lev)
         x(2*i) = xsave(i,lev)
   78 continue
      go to 30
!
!    ***   stage 8   ***   finalize and return
!
   80 ans = ans + cor11
!
!    make sure errest not less than roundoff level.
!
      if (errest .eq. 0.0_wp) return
   82 temp = abs(ans) + errest
      if (temp .ne. abs(ans)) return
      errest = 2.0_wp*errest
      go to 82
end associate
end function


end module
