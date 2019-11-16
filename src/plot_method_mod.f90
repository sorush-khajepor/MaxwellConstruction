module plot_method_mod
    use method_mod
    implicit none
!   gp = gnuplot
!   gpd = gnuplot data
!   gpc = gnuplot config
!   res = resolution
    type plot_method_t
        class(method_t),pointer :: method
        Character (len=100):: fname,gpcfname,tit,xl,yl
        real(kind=wp) :: x_st,x_nd
        integer :: funit,res
    contains
        procedure :: init => init_plot_method
        procedure :: f => curve_function
        procedure :: write_gpd_file => write_gpd_file
        procedure :: write_gpc_file=> write_gpc_file
        procedure :: plot => plot_gp
    end type
contains

function curve_function (this,x) result (y)
    class(plot_method_t) :: this
    real(kind=wp) :: x,y
    y = this%method%func(x)
end function

subroutine init_plot_method(this,method,x_st,x_nd,res,fname,tit,xl,yl)
    class(plot_method_t) :: this
    class(method_t),target :: method
    real(kind=wp) :: x_st,x_nd
    integer :: res
    Character (len=*):: fname,tit,xl,yl
    this%method=>method
    this%x_st=x_st
    this%x_nd=x_nd
    this%res=res
    this%fname=fname
    this%tit=tit
    this%xl=xl
    this%yl=yl
    this%gpcfname="plot.cfg"
end subroutine

subroutine write_gpd_file(this)
    class(plot_method_t) :: this
    integer i
    real(kind=wp) :: x,dx
    associate(x_st=>this%x_st,x_nd=>this%x_nd,funit=>this%funit,res=>this%res)

    open(newunit=funit,file= trim(this%fname))
    write(funit,'(A3,2A25)') ' ',this%xl,this%yl
    dx = (x_nd-x_st)/res
    do i = 0, res
        x = x_st + i * dx
        write(funit,*) x, this%f(x)
    end do
    close (funit)
    end associate
end subroutine

subroutine write_gpc_file(this)
    class(plot_method_t) :: this
    integer :: funit
! Gnuplot configuration to produce GIF file
    open(newunit=funit,file=trim(this%gpcfname))
    write(funit,*) "set terminal png "
    write(funit,*) 'set output "'//trim(this%fname)//'.png"'
    write(funit,*) "set nokey"
    !write(funit,*) 'set xrange[0:', xdim,']'
    !write(funit,*) 'set yrange[1:', ydim,']'
    write(funit,*) 'set title  "'//trim(this%tit)//'"'
    write(funit,*) 'set xlabel "'//trim(this%xl)//'"'
    write(funit,*) 'set ylabel "'//trim(this%yl)//'"'

    write(funit,*) "p '"//trim(this%fname)//"' using 1:2 w lp "

    close(funit)
end subroutine

subroutine plot_gp(this)
    class(plot_method_t) :: this

    call this%write_gpd_file()
    call this%write_gpc_file()
    call SYSTEM('gnuplot '//trim(this%gpcfname))
    call SYSTEM('gpicview '//trim(this%fname)//'.png &')

end subroutine

end module plot_method_mod
