
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_genPlotfileParaview

    use RST_model
    use RST_globalData

    implicit none

contains

    subroutine genPlotfileParaview()

        character(len=30) :: fpparaviewvtk
        character(len=15) :: charLx, charLy, charLz
        character(len=5) :: charnx, charny, charnz
        character(len=15) :: charsx, charsy, charsz
        character(len=15) :: charox, charoy, charoz
        character(len=15) :: charpn
        character(len=15) :: charPw
        real(kind=8) :: sx, sy, sz, ox, oy, oz
        logical :: alive
        integer :: ierr, i, j, k

        if(myid == 0) then

            fpparaviewvtk = trim(adjustl(soludoc))//"/paraviewplot.vtk"

            open(unit=10, file=trim(adjustl(fpparaviewvtk)), status='replace')

            write(10, fmt="(a)") "# vtk DataFile Version 3.0"
            write(10, fmt="(a)") "vtk output"
            write(10, fmt="(a)") "ASCII"
            write(10, fmt="(a)") "DATASET STRUCTURED_POINTS"
            write(charnx,'(i5)') nx
            write(charny,'(i5)') ny
            write(charnz,'(i5)') nz
            write(10, fmt="(a)") "DIMENSIONS "//charnx//" "//charny//" "//charnz
            write(charLx,'(es15.8)') Lx
            write(charLy,'(es15.8)') Ly
            write(charLz,'(es15.8)') Lz
            sx = Lx/nx
            sy = Ly/ny
            sz = Lz/nz
            write(charsx,'(es15.8)') sx
            write(charsy,'(es15.8)') sy
            write(charsz,'(es15.8)') sz
            ox = sx/2
            oy = sy/2
            oz = sz/2
            write(charox,'(es15.8)') ox
            write(charoy,'(es15.8)') oy
            write(charoz,'(es15.8)') oz
            write(10, fmt="(a)") "ORIGIN "//charox//" "//charoy//" "//charoz
            write(10, fmt="(a)") "SPACING "//charsx//" "//charsy//" "//charsz
            write(charpn,'(i12)') nx*ny*nz
            write(10, fmt="(a)") "POINT_DATA "//charpn
            write(10, fmt="(a)") "SCALARS pre-point float"
            write(10, fmt="(a)") "LOOKUP_TABLE default"
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        write(charPw,'(es15.8)') Pw(i,j,k)
                        write(10, fmt="(a)") charPw//" "
                    end do
                end do
            end do

            close(10)

        end if

    end subroutine genPlotfileParaview

end module RST_genPlotfileParaview
