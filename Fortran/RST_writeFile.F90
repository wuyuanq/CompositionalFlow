
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_writeFile

    use RST_model
    use RST_globalData
    use RST_globalFlashData

    implicit none

contains

    subroutine writeFile()

        logical :: alive
        integer :: i, j, m
        character(len=50) :: fpwtxt, fuxtxt, fuytxt, fswtxt, fxitxt
        character(len=50), dimension(:), pointer :: fmftxt
        character :: charm

        fpwtxt = trim(adjustl(soludoc))//"/soln_cmp2PhFlw_Pw_raw.txt"
        fuxtxt = trim(adjustl(soludoc))//"/soln_cmp2PhFlw_Ux_raw.txt"
        fuytxt = trim(adjustl(soludoc))//"/soln_cmp2PhFlw_Uy_raw.txt"
        fswtxt = trim(adjustl(soludoc))//"/soln_cmp2PhFlw_Sw_raw.txt"
        fxitxt = trim(adjustl(soludoc))//"/soln_cmp2PhFlw_c_raw.txt"
        allocate(fmftxt(Nc))
        do m = 1, Nc
            write(charm,'(i1)') m
            fmftxt(m) = trim(adjustl(soludoc))//"/soln_cmp2PhFlw_x"//charm//"_raw.txt"
        end do

        open(unit=10, file=trim(adjustl(fpwtxt)), status='replace')
        do j = 2, nx+1
            do i = 2, ny+1
                write(10, fmt="(es15.8)") Pw(i,j)
            end do
        end do
        close(10)

        open(unit=10, file=trim(adjustl(fuxtxt)), status='replace')
        do j = 1, nx+1
            do i = 1, ny
                write(10, fmt="(es15.8)") Uwx(i,j)+Unx(i,j)
            end do
        end do
        close(10)

        open(unit=10, file=trim(adjustl(fuytxt)), status='replace')
        do j = 1, nx
            do i = 1, ny+1
                write(10, fmt="(es15.8)") Uwy(i,j)+Uny(i,j)
            end do
        end do
        close(10)

        do m = 1, Nc
            open(unit=10, file=trim(adjustl(fmftxt(m))), status='replace')
            do j = 1, nx
                do i = 1, ny
                    write(10, fmt="(es15.8)") z(m,i,j)
                end do
            end do
            close(10)
        end do

        open(unit=10, file=trim(adjustl(fswtxt)), status='replace')
        do j = 1, nx
            do i = 1, ny
                write(10, fmt="(es15.8)") Sw(i,j)
            end do
        end do
        close(10)

        open(unit=10, file=trim(adjustl(fxitxt)), status='replace')
        do j = 1, nx
            do i = 1, ny
                write(10, fmt="(es15.8)") xiW(i,j)*Sw(i,j)+xiN(i,j)*(1-Sw(i,j))
            end do
        end do
        close(10)

        deallocate(fmftxt)

    end subroutine writeFile

end module RST_writeFile
