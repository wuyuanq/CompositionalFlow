
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com
!!$
!!$ The module provides a subroutine which computes out an optimal scheme to
!!$ allocate the cells to different processes so that the communication cost
!!$ is the lowest. 

module RST_proceAlloc

contains

    subroutine proceAlloc(Np, local_nx, local_ny, local_nz, local_pncols, local_pnrows, local_pnlays)

        implicit none
        integer, intent(in) :: Np ! number of processes
        integer, intent(in) :: local_nx
        integer, intent(in) :: local_ny
        integer, intent(in) :: local_nz
        integer, intent(out) :: local_pncols
        integer, intent(out) :: local_pnrows
        integer, intent(out) :: local_pnlays
        integer(kind=8) :: minlen, curlen
        integer :: i, j, k, l

        minlen = local_nx*local_ny*local_nz*6

        do i = 1, Np
            if(mod(Np,i) == 0) then
                j = Np/i
                do k = 1, j
                    if(mod(j,k) == 0) then
                        l = j/k
                        if((mod(local_nx,i) == 0).and.(mod(local_ny,k) == 0).and.(mod(local_nz,l) == 0)) then
                            curlen = (i-1)*local_ny*local_nz + (k-1)*local_nx*local_nz + (l-1)*local_nx*local_ny
                            if(curlen < minlen) then
                                minlen = curlen
                                local_pncols = i
                                local_pnrows = k
                                local_pnlays = l
                            end if
                        end if
                    end if
                end do
            end if
        end do

        if(minlen == local_nx*local_ny*local_nz*6) then
            print *, 'The number of processes is not right. Please adjust it.'
            print *, 'You must make sure that each process has the same number of cells.'
            stop
        end if

    end subroutine proceAlloc

end module RST_proceAlloc