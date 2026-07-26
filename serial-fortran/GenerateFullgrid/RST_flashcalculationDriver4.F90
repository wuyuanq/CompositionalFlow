
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_flashcalculationDriver

    use RST_globalFlashData
    use RST_Flashcalculation_fullgrid_3c

contains

    subroutine flashcalculationDriver()

        implicit none

        integer, parameter :: GRIDSIZE_P = 33
        integer, parameter :: GRIDSIZE_Z1 = 513
        integer, parameter :: GRIDSIZE_Z2 = 513
        real(kind=8), parameter :: PMIN = 1.9D6!1.984D6
        real(kind=8), parameter :: PMAX = 1.6D8!2.112D6
        real(kind=8), parameter :: PSTEP = (PMAX-PMIN)/(GRIDSIZE_P-1)

        real(kind=8), dimension(:), pointer :: z, xW, xN, v
        real(kind=8) :: xiW, xiN, densiW, densiN, Sw, Cf
        logical :: isW, isN
        integer :: i, j, m

        open(unit=110, file="dummyfullgrid_3c.txt", status='replace')

        allocate(z(Nc))
        allocate(xW(Nc))
        allocate(xN(Nc))
        allocate(v(Nc))

        allocate(xtable(Nc,TABLESIZE))
        allocate(ytable(Nc,TABLESIZE))
        allocate(xiLtable(TABLESIZE))
        allocate(xiGtable(TABLESIZE))
        allocate(rhoLtable(TABLESIZE))
        allocate(rhoGtable(TABLESIZE))
        allocate(sLtable(TABLESIZE))
        allocate(vtable(Nc,TABLESIZE))
        allocate(Cftable(TABLESIZE))
        allocate(isWtable(TABLESIZE))
        allocate(isNtable(TABLESIZE))

        open(unit=120, file=FULLGRIDPREFIX//"xW1.txt", status='old')
        open(unit=130, file=FULLGRIDPREFIX//"xW2.txt", status='old')
        open(unit=140, file=FULLGRIDPREFIX//"xW3.txt", status='old')
        open(unit=150, file=FULLGRIDPREFIX//"xN1.txt", status='old')
        open(unit=160, file=FULLGRIDPREFIX//"xN2.txt", status='old')
        open(unit=170, file=FULLGRIDPREFIX//"xN3.txt", status='old')
        open(unit=180, file=FULLGRIDPREFIX//"xiW.txt", status='old')
        open(unit=190, file=FULLGRIDPREFIX//"xiN.txt", status='old')
        open(unit=200, file=FULLGRIDPREFIX//"densiW.txt", status='old')
        open(unit=210, file=FULLGRIDPREFIX//"densiN.txt", status='old')
        open(unit=220, file=FULLGRIDPREFIX//"sW.txt", status='old')
        open(unit=230, file=FULLGRIDPREFIX//"v1.txt", status='old')
        open(unit=240, file=FULLGRIDPREFIX//"v2.txt", status='old')
        open(unit=250, file=FULLGRIDPREFIX//"v3.txt", status='old')
        open(unit=260, file=FULLGRIDPREFIX//"Cf.txt", status='old')
        open(unit=270, file=FULLGRIDPREFIX//"liquid.txt", status='old')
        open(unit=280, file=FULLGRIDPREFIX//"gas.txt", status='old')

        read(120,*) xtable(1,1:TABLESIZE)
        read(130,*) xtable(2,1:TABLESIZE)
        read(140,*) xtable(3,1:TABLESIZE)
        read(150,*) ytable(1,1:TABLESIZE)
        read(160,*) ytable(2,1:TABLESIZE)
        read(170,*) ytable(3,1:TABLESIZE)
        read(180,*) xiLtable(1:TABLESIZE)
        read(190,*) xiGtable(1:TABLESIZE)
        read(200,*) rhoLtable(1:TABLESIZE)
        read(210,*) rhoGtable(1:TABLESIZE)
        read(220,*) sLtable(1:TABLESIZE)
        read(230,*) vtable(1,1:TABLESIZE)
        read(240,*) vtable(2,1:TABLESIZE)
        read(250,*) vtable(3,1:TABLESIZE)
        read(260,*) Cftable(1:TABLESIZE)
        read(270,*) isWtable(1:TABLESIZE)
        read(280,*) isNtable(1:TABLESIZE)

        do i = 0, GRIDSIZE_P-1
            do j = 0, GRIDSIZE_Z1-1
                do m = 0, GRIDSIZE_Z2-1

                    print *, i, j, m

                    z(1) = j*1.D0/(GRIDSIZE_Z1-1)
                    z(2) = m*1.D0/(GRIDSIZE_Z2-1)
                    z(3) = 1 - z(1) - z(2)

                    if(z(3) >= 0) then
                        call flashcalculation_fullgrid_3c( PMIN+i*PSTEP, z, xW, xN, xiW, xiN, &!
                            densiW, densiN, Sw, v, Cf, isW, isN )

                        write(110, fmt="(es12.5)") PMIN+i*PSTEP, z, xW, xN, xiW, xiN, densiW, densiN, Sw, v, Cf
                        write(110, fmt="(L1)") isW, isN

                    end if

                end do
            end do
        end do

        deallocate(z)
        deallocate(xW)
        deallocate(xN)
        deallocate(v)

        close(110)
        close(120)
        close(130)
        close(140)
        close(150)
        close(160)
        close(170)
        close(180)
        close(190)
        close(200)
        close(210)
        close(220)
        close(230)
        close(240)
        close(250)
        close(260)
        close(270)
        close(280)

    end subroutine flashcalculationDriver

end module RST_flashcalculationDriver
