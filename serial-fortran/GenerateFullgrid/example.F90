
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com
!!$
!!$ This is a 2 component example.

program example

    use RST_flashcalculationDriver

    implicit none

    integer :: i

    Nc = 2
    Temp = 2.6D2

    allocate(ct(Nc))
    ct(1) = 190
    ct(2) = 370
    allocate(cp(Nc))
    cp(1) = 4.6*1.D6
    cp(2) = 4.2*1.D6
    allocate(af(Nc))
    af(1) = 0.01
    af(2) = 0.15
    allocate(mw(Nc))
    mw(1) = 0.016
    mw(2) = 0.044
    allocate(cv(Nc))
    cv(1) = 0.0062
    cv(2) = 0.0045
    allocate(psatA(Nc))
    psatA(1) = 6.69561
    psatA(2) = 6.82973
    allocate(psatB(Nc))
    psatB(1) = 405.420
    psatB(2) = 813.2
    allocate(psatC(Nc))
    psatC(1) = 267.777
    psatC(2) = 248
    allocate(delta(Nc,Nc))
    delta(:,:) = 0.0
    delta(1,2) = 0.036
    delta(2,1) = delta(1,2)

    call flashcalculationDriver()

end program example
