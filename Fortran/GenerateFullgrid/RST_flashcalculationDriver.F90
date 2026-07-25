
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_flashcalculationDriver

    use RST_flashcalculation

    character(len=*), parameter :: FULLGRIDPREFIX = '/Users/yuanqingwu/research/CompositionalFlow/Fortran/Fullgrid_2c/'

contains

    subroutine flashcalculationDriver()

        implicit none
        
        integer, parameter :: GRIDSIZE_P = 129
        integer, parameter :: GRIDSIZE_Z1 = 129
        real(kind=8), parameter :: PMIN = 1.9D6!1.800D6!
        real(kind=8), parameter :: PMAX = 2.1D6!3.300D6!
        real(kind=8), parameter :: PSTEP = (PMAX-PMIN)/(GRIDSIZE_P-1)
        real(kind=8), dimension(:), pointer :: z, xW, xN, v
        real(kind=8) :: xiW, xiN, densiW, densiN, Sw, Cf
        logical :: isW, isN, isRea
        integer :: i, j

        open(unit=120, file=FULLGRIDPREFIX//"P.txt", status='replace')
        open(unit=130, file=FULLGRIDPREFIX//"z1.txt", status='replace')
        open(unit=140, file=FULLGRIDPREFIX//"z2.txt", status='replace')
        open(unit=150, file=FULLGRIDPREFIX//"xW1.txt", status='replace')
        open(unit=160, file=FULLGRIDPREFIX//"xW2.txt", status='replace')
        open(unit=170, file=FULLGRIDPREFIX//"xN1.txt", status='replace')
        open(unit=180, file=FULLGRIDPREFIX//"xN2.txt", status='replace')
        open(unit=190, file=FULLGRIDPREFIX//"xiW.txt", status='replace')
        open(unit=200, file=FULLGRIDPREFIX//"xiN.txt", status='replace')
        open(unit=210, file=FULLGRIDPREFIX//"densiW.txt", status='replace')
        open(unit=220, file=FULLGRIDPREFIX//"densiN.txt", status='replace')
        open(unit=230, file=FULLGRIDPREFIX//"sW.txt", status='replace')
        open(unit=240, file=FULLGRIDPREFIX//"v1.txt", status='replace')
        open(unit=250, file=FULLGRIDPREFIX//"v2.txt", status='replace')
        open(unit=260, file=FULLGRIDPREFIX//"Cf.txt", status='replace')
        open(unit=270, file=FULLGRIDPREFIX//"liquid.txt", status='replace')
        open(unit=280, file=FULLGRIDPREFIX//"gas.txt", status='replace')

        allocate(z(Nc))
        allocate(xW(Nc))
        allocate(xN(Nc))
        allocate(v(Nc))

        do i = 0, GRIDSIZE_P-1
            do j = 0, GRIDSIZE_Z1-1

                !print *, i, j

                z(1) = j*1.D0/(GRIDSIZE_Z1-1)
                z(2) = 1 - z(1)

                call flashcalculation( PMIN+i*PSTEP, z, xW, xN, xiW, xiN, &!
                    densiW, densiN, Sw, v, Cf, isW, isN, isRea )

                write(120, fmt="(es24.16)") PMIN+i*PSTEP
                write(130, fmt="(es24.16)") z(1)
                write(140, fmt="(es24.16)") z(2)
                write(150, fmt="(es24.16)") xW(1)
                write(160, fmt="(es24.16)") xW(2)
                write(170, fmt="(es24.16)") xN(1)
                write(180, fmt="(es24.16)") xN(2)
                write(190, fmt="(es24.16)") xiW
                write(200, fmt="(es24.16)") xiN
                write(210, fmt="(es24.16)") densiW
                write(220, fmt="(es24.16)") densiN
                write(230, fmt="(es24.16)") Sw
                write(240, fmt="(es24.16)") v(1)
                write(250, fmt="(es24.16)") v(2)
                write(260, fmt="(es24.16)") Cf
                if(isW) then
                    write(270, fmt="(i1)") 1
                else
                    write(270, fmt="(i1)") 0
                end if
                if(isN) then
                    write(280, fmt="(i1)") 1
                else
                    write(280, fmt="(i1)") 0
                end if

            end do
        end do

        deallocate(z)
        deallocate(xW)
        deallocate(xN)
        deallocate(v)

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
