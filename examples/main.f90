program main
    use min_MI_method_mod
    use plot_method_mod
    implicit none
    type(brent_t) :: brent_minmi
    type(m4minmi_t),target :: m4minmi
    class (eos_t), allocatable :: eos
    real(kind=wp) :: p_st,p_nd,psat,rhog_sat,rhol_sat,rho_st,rho_nd
    real(kind=wp) :: err1,err2,Tr,dTr,Tr_st,Tr_plot,tmp1,tmp2,a,b
    integer i,Tr_steps,funit_sp
    character(len=50) :: fname_satpoints,tit
    type(plot_method_t) :: pm
    type(m4pres_t) :: m4pres
    real(kind=wp),parameter ::  eps = epsilon(1._wp)

!   ============ User Input ======================================
!   set reduced temperatures, and plot Tr
!   Start reduced temperature
    Tr_st = 0.3_wp
!   Reduced temperature step size    
    dTr = 0.01_wp
!   Number of temperature steps
    Tr_steps = 80
!   Draw p-rho plot at this temperature (needs gnuplot installed)     
    Tr_plot = 0.7_wp

!   Allocating the type of equaiton of state (VW,CS,PR,RK,SRK,...)
!   Uncomment one equation of state
    allocate(srk_t::eos)
    !allocate(rk_t::eos)
    !allocate(cs_t::eos)
    !allocate(pr_t::eos)
    !allocate(vw_t::eos)

!   =============End User Input ===================================


!   Initializing eos based on their type.
    select type (eospt=>eos)
        type is (vw_t)
            a=0.01_wp
            b=0.2_wp
            call eospt%init(a=a,b=b,Tr=Tr_st)
            rho_st = eps
            rho_nd = 1._wp/b - 100*eps
        type is (cs_t)
            a=0.01_wp
            b=0.2_wp
            call eospt%init(a=a,b=b,Tr=Tr_st)
            rho_st= eps
            rho_nd= 4._wp/b - 100*eps
        type is (rk_t)
            a=0.01_wp
            b=0.2_wp
            call eospt%init(a=a,b=b,Tr=Tr_st)
            rho_st = eps
            rho_nd = 1._wp/b- 100*eps
        type is (srk_t)
            a=0.01_wp
            b=0.2_wp
            call eospt%init(a=a,b=b,omega=0.344_wp,Tr=Tr_st)
            rho_st = eps
            rho_nd = 1._wp/b- 100*eps
        type is (pr_t)
            a=0.01_wp
            b=0.2_wp
            call eospt%init(a=a,b=b,omega=0.344_wp,Tr=Tr_st)
            rho_st = eps
            rho_nd = 1._wp/b- 100*eps
        class default
            stop "Error: EOS type is not recognized!"
    end select

!   check critical points
    tmp1 = eos%getdp_drho(eos%rho_c)
    tmp2 = eos%getd2p_drho2(eos%rho_c)
    if (tmp1+tmp2.gt.100*eps) then
        print*, "Epsilon = ", eps
        print*, "dp/drho=",tmp1,"d2p/drho=",tmp2
        stop "Crtitcal rho and T can not make dp/drho=0 and d2p/drho2 = 0!"
    else
        print*, "Error: Critical point is successfully set."
    end if

!   initializing output file of saturation points

    write(fname_satpoints,'(A50)') "satpoints_"//trim(eos%getmyname())//".txt"
    open(newunit=funit_sp,file= trim(adjustl(fname_satpoints)))
    write(funit_sp,'(A1,A6,A5,f6.3,A5,F6.3)') "#", eos%getmyname(),"a =",eos%a,"b =",eos%b
    write(funit_sp,'(A1,A5,f16.10,A11,f16.10,A9,f16.10)') "#", "Tc =",eos%Tc,"    rho_c =",eos%rho_c,"    Pc =", eos%pc
    write(funit_sp,'(A1,A6,7A11,2A16)')"#","Tr ","T  ","rhog_sat","rhol_sat","Psat  ","rhog_sat_r",&
                                "rhol_sat_r","psat_r","Integral_Err","BC_Error"

!   Loop over reduced temperatures Tr
    Tr=Tr_st-dTr
    do i =1,Tr_steps

!       set new temperature
        Tr=Tr+dTr
        call eos%setT(Tr=Tr)

!       initializing the method for minimization of Maxwell integral
        call m4minmi%init(rho_st=rho_st,rho_nd=rho_nd,eos=eos)

!       p is the x of brent for minimization of the Maxwell integral
!       negative p gives no vapor density from vapor branch (look at the graph)
!       zero is not good for a starting point (out of experience), 3*eps is chosen
        p_st = max(100*eps,m4minmi%pmin)
        p_nd = m4minmi%pmax-10*eps

        call brent_minmi%init(method=m4minmi, mbound_st=p_st,mbound_nd=p_nd)

!       the answer of brent is p(or x) = P_saturaion
        Psat=brent_minmi%getroot()

!       finding the intersection of psat and EOS (which are rhog and rhol)
        rhog_sat = m4minmi%getrhog(Psat)
        rhol_sat = m4minmi%getrhol(Psat)

!       calculate Maxwell integral with the aid of results which should be zero
        err1 = m4minmi%getintegral_err(rhog_sat,rhol_sat,psat)
!       boundary conditoins of Max. Integ. i.e. P(rhog)=P(rhol)=P0
        err2 = abs(eos%getp(rhog_sat)-eos%getp(rhol_sat))

!       write saturation points in table format for different reduced Temperatures (Tr)
        write(funit_sp,'(1f7.3,7f15.9,2e16.4e3)') Tr,eos%T,rhog_sat,rhol_sat,psat, &
                            rhog_sat/eos%rho_c, rhol_sat/eos%rho_c, psat/eos%pc,err1,err2

!       plot p vs. rho
        if(abs(Tr-Tr_plot).lt.10*eps) then
                write(tit,'(A5,A4,f5.3,A4,F5.3,A4,F5.3)') eos%getmyname(),"  a=",eos%a,"  b=",eos%b,"Tr=",Tr
                call m4pres%init(eos=eos)
                call pm%init(method=m4pres,x_st=0.1_wp*m4minmi%rhomax, &
                                           x_nd=1.35_wp*m4minmi%rhomin, &
                                           res=100,fname='p_rho_'//trim(eos%getmyname()), &
                                           tit=trim(tit),xl='Density',yl='Pressure')
                call pm%plot()
                write(*,'(A27,F6.3)')"P vs. rho is plotted at Tr=",Tr

        end if
    write(*,'(A3,f5.3,A15)')"Tr=",Tr,"is calculated."
    end do

    close(funit_sp)
end program main
