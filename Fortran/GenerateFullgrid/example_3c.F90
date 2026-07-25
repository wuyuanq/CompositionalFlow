
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com
!!$
!!$ This is a two-phase case. Inject to 1 and output from 2.

program example

    use RST_globalFlashData
    use RST_flashcalculationDriver

    implicit none

    integer :: i

    Nc = 3
    Temp = 480.0

    ! 1 is CO2, 2 is CH4 and 3 is C3H8
    allocate(ct(Nc))
    ct(1) = 304.4
    ct(2) = 190
    ct(3) = 370
    allocate(cp(Nc))
    cp(1) = 7.4*1.D6
    cp(2) = 4.6*1.D6
    cp(3) = 4.2*1.D6
    allocate(af(Nc))
    af(1) = 0.23
    af(2) = 0.01
    af(3) = 0.15
    allocate(mw(Nc))
    mw(1) = 0.044
    mw(2) = 0.016
    mw(3) = 0.044
    allocate(cv(Nc))
    cv(1) = 0.0021
    cv(2) = 0.0062
    cv(3) = 0.0045
    allocate(psatA(Nc))
    psatA(1) = 6.81228
    psatA(2) = 6.69561
    psatA(3) = 6.82973
    allocate(psatB(Nc))
    psatB(1) = 1301.679
    psatB(2) = 405.420
    psatB(3) = 813.2
    allocate(psatC(Nc))
    psatC(1) = 269.506
    psatC(2) = 267.777
    psatC(3) = 248
    allocate(delta(Nc,Nc))
    delta(:,:) = 0.0
    delta(1,2) = 0.15
    delta(2,1) = delta(1,2)
    delta(1,3) = 0.1239
    delta(3,1) = delta(1,3)
    delta(2,3) = 0.036
    delta(3,2) = delta(2,3)

    call flashcalculationDriver()

end program example
