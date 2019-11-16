Module kinds_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! module setting global working precision for Real numbers and
! double precision for Integer numbers at compile time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Implicit None

  Integer, Parameter :: wp = Selected_Real_Kind(15,307)
  Integer, Parameter :: ip = Selected_Int_Kind(12)
  Real    (kind = wp) :: wpSample
  Integer (kind = ip) :: ipSample
  Integer             :: intSample
! To define a variable:
!        Real( Kind = wp )         ::
! To define a constant:
!        3.14_wp
! Type conversion:
!   a = Real( <integer value>, wp)

End Module kinds_f
