
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
    include "mpif.h"

contains

    subroutine writeFile()

        implicit none
        
        character(len=50) :: fpwtxt, fuxtxt, fuytxt, fswtxt, fxitxt
        character(len=50), dimension(:), pointer :: fmftxt
        character :: charm
        real(kind=8), dimension(:,:), pointer :: allPw
        real(kind=8), dimension(:,:), pointer :: alluwx
        real(kind=8), dimension(:,:), pointer :: alluwy
        real(kind=8), dimension(:,:), pointer :: allunx
        real(kind=8), dimension(:,:), pointer :: alluny
        real(kind=8), dimension(:,:,:), pointer :: allx
        real(kind=8), dimension(:,:), pointer :: allSw
        real(kind=8), dimension(:,:), pointer :: allxiW
        real(kind=8), dimension(:,:), pointer :: allxiN
        integer :: buffersize
        real(kind=8), dimension(:), pointer :: buffer
        real(kind=8), dimension(:), pointer :: Pwsent, uwxsent, uwysent, unxsent, unysent, &!
            zsent, Swsent, xiWsent, xiNsent
        integer :: p_xlower, p_xupper, p_ylower, p_yupper
        integer :: p_prow, p_pcol
        integer :: position
        integer :: i, j, m, n, prociteration

        integer :: ierr
        integer :: request
        integer :: status(MPI_STATUS_SIZE)
        integer :: commbuffer_size = MAX_COMMBUF
        real(kind=8) :: commbuffer(MAX_COMMBUF)

        buffersize = (Nc+8)*local_size

        allocate(buffer(buffersize))
        call MPI_BUFFER_ATTACH(commbuffer,commbuffer_size,ierr)

        if(myid /= 0) then

            allocate(Pwsent(local_size))
            allocate(uwxsent(local_size))
            allocate(uwysent(local_size))
            allocate(unxsent(local_size))
            allocate(unysent(local_size))
            allocate(zsent(local_size*Nc))
            allocate(Swsent(local_size))
            allocate(xiWsent(local_size))
            allocate(xiNsent(local_size))

            n = 1
            do j = 1, localncols
                do i = 1, localnrows
                    Pwsent(n) = Pw(i,j)
                    n = n + 1
                end do
            end do
            n = 1
            do j = 1, localncols
                do i = 1, localnrows
                    uwxsent(n) = Uwx(i,j)
                    n = n + 1
                end do
            end do
            n = 1
            do j = 1, localncols
                do i = 1, localnrows
                    uwysent(n) = Uwy(i,j)
                    n = n + 1
                end do
            end do
            n = 1
            do j = 1, localncols
                do i = 1, localnrows
                    unxsent(n) = Unx(i,j)
                    n = n + 1
                end do
            end do
            n = 1
            do j = 1, localncols
                do i = 1, localnrows
                    unysent(n) = Uny(i,j)
                    n = n + 1
                end do
            end do
            n = 1
            do j = 1, localncols
                do i = 1, localnrows
                    do m = 1, Nc
                        zsent(n) = z(m,i,j)
                        n = n + 1
                    end do
                end do
            end do
            n = 1
            do j = 1, localncols
                do i = 1, localnrows
                    Swsent(n) = Sw(i,j)
                    n = n + 1
                end do
            end do
            n = 1
            do j = 1, localncols
                do i = 1, localnrows
                    xiWsent(n) = xiW(i,j)
                    n = n + 1
                end do
            end do
            n = 1
            do j = 1, localncols
                do i = 1, localnrows
                    xiNsent(n) = xiN(i,j)
                    n = n + 1
                end do
            end do

            position = 0
            call MPI_PACK(Pwsent, local_size, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, &!
                MPI_COMM_WORLD, ierr)
            call MPI_PACK(uwxsent, local_size, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, &!
                MPI_COMM_WORLD, ierr)
            call MPI_PACK(uwysent, local_size, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, &!
                MPI_COMM_WORLD, ierr)
            call MPI_PACK(unxsent, local_size, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, &!
                MPI_COMM_WORLD, ierr)
            call MPI_PACK(unysent, local_size, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, &!
                MPI_COMM_WORLD, ierr)
            call MPI_PACK(zsent, local_size*Nc, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, &!
                MPI_COMM_WORLD, ierr)
            call MPI_PACK(Swsent, local_size, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, &!
                MPI_COMM_WORLD, ierr)
            call MPI_PACK(xiWsent, local_size, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, &!
                MPI_COMM_WORLD, ierr)
            call MPI_PACK(xiNsent, local_size, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, &!
                MPI_COMM_WORLD, ierr)

            call MPI_IBSEND(buffer,buffersize,MPI_DOUBLE_PRECISION,0,myid+num_procs,MPI_COMM_WORLD,request,ierr)
            call MPI_WAIT(request, status, ierr)

            deallocate(Pwsent)
            deallocate(uwxsent)
            deallocate(uwysent)
            deallocate(unxsent)
            deallocate(unysent)
            deallocate(zsent)
            deallocate(Swsent)
            deallocate(xiWsent)
            deallocate(xiNsent)

        else

            allocate(allPw(ny,nx))
            allocate(alluwx(ny,nx+1))
            allocate(alluwy(ny+1,nx))
            allocate(allunx(ny,nx+1))
            allocate(alluny(ny+1,nx))
            allocate(allx(Nc,ny,nx))
            allocate(allSw(ny,nx))
            allocate(allxiW(ny,nx))
            allocate(allxiN(ny,nx))

            alluwx(1:ny,nx+1) = UwBdryX(2, 1:ny)
            allunx(1:ny,nx+1) = UnBdryX(2, 1:ny)
            alluwy(ny+1,1:nx) = UwBdryY(1:nx, 2)
            alluny(ny+1,1:nx) = UnBdryY(1:nx, 2)

            do j = 1, localncols
                do i = 1, localnrows
                    allPw(i,j) = Pw(i,j)
                    alluwx(i,j) = Uwx(i,j)
                    alluwy(i,j) = Uwy(i,j)
                    allunx(i,j) = Unx(i,j)
                    alluny(i,j) = Uny(i,j)
                    allSw(i,j) = Sw(i,j)
                    allxiW(i,j) = xiW(i,j)
                    allxiN(i,j) = xiN(i,j)
                end do
            end do
            do j = 1, localncols
                do i = 1, localnrows
                    do m = 1, Nc
                        allx(m,i,j) = z(m,i,j)
                    end do
                end do
            end do

            do prociteration = 1, num_procs-1

                p_prow = mod(prociteration,pnrows)+1
                p_pcol = prociteration/pnrows+1

                p_xlower = (p_pcol-1)*localncols+1
                p_xupper = p_pcol*localncols
                p_ylower = (p_prow-1)*localnrows+1
                p_yupper = p_prow*localnrows

                call MPI_RECV(buffer, buffersize, MPI_DOUBLE_PRECISION, prociteration, &!
                    prociteration+num_procs, MPI_COMM_WORLD, status, ierr)

                n = 1
                do j = p_xlower, p_xupper
                    do i = p_ylower, p_yupper
                        allPw(i,j) = buffer(n)
                        n = n + 1
                    end do
                end do
                do j = p_xlower, p_xupper
                    do i = p_ylower, p_yupper
                        alluwx(i,j) = buffer(n)
                        n = n + 1
                    end do
                end do
                do j = p_xlower, p_xupper
                    do i = p_ylower, p_yupper
                        alluwy(i,j) = buffer(n)
                        n = n + 1
                    end do
                end do
                do j = p_xlower, p_xupper
                    do i = p_ylower, p_yupper
                        allunx(i,j) = buffer(n)
                        n = n + 1
                    end do
                end do
                do j = p_xlower, p_xupper
                    do i = p_ylower, p_yupper
                        alluny(i,j) = buffer(n)
                        n = n + 1
                    end do
                end do
                do j = p_xlower, p_xupper
                    do i = p_ylower, p_yupper
                        do m = 1, Nc
                            allx(m,i,j) = buffer(n)
                            n = n + 1
                        end do
                    end do
                end do
                do j = p_xlower, p_xupper
                    do i = p_ylower, p_yupper
                        allSw(i,j) = buffer(n)
                        n = n + 1
                    end do
                end do
                do j = p_xlower, p_xupper
                    do i = p_ylower, p_yupper
                        allxiW(i,j) = buffer(n)
                        n = n + 1
                    end do
                end do
                do j = p_xlower, p_xupper
                    do i = p_ylower, p_yupper
                        allxiN(i,j) = buffer(n)
                        n = n + 1
                    end do
                end do

            end do

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
            do j = 1, nx
                do i = 1, ny
                    write(10, fmt="(es15.8)") allPw(i,j)
                end do
            end do
            close(10)

            open(unit=10, file=trim(adjustl(fuxtxt)), status='replace')
            do j = 1, nx+1
                do i = 1, ny
                    write(10, fmt="(es15.8)") alluwx(i,j)+allunx(i,j)
                end do
            end do
            close(10)

            open(unit=10, file=trim(adjustl(fuytxt)), status='replace')
            do j = 1, nx
                do i = 1, ny+1
                    write(10, fmt="(es15.8)") alluwy(i,j)+alluny(i,j)
                end do
            end do
            close(10)

            do m = 1, Nc
                open(unit=10, file=trim(adjustl(fmftxt(m))), status='replace')
                do j = 1, nx
                    do i = 1, ny
                        write(10, fmt="(es15.8)") allx(m,i,j)
                    end do
                end do
                close(10)
            end do

            open(unit=10, file=trim(adjustl(fswtxt)), status='replace')
            do j = 1, nx    
                do i = 1, ny
                    write(10, fmt="(es15.8)") allSw(i,j)
                end do
            end do
            close(10)

            open(unit=10, file=trim(adjustl(fxitxt)), status='replace')
            do j = 1, nx
                do i = 1, ny
                    write(10, fmt="(es15.8)") allxiW(i,j)*allSw(i,j)+allxiN(i,j)*(1-allSw(i,j))
                end do
            end do
            close(10)

            deallocate(fmftxt)

            deallocate(allPw)
            deallocate(alluwx)
            deallocate(alluwy)
            deallocate(allunx)
            deallocate(alluny)
            deallocate(allx)
            deallocate(allSw)
            deallocate(allxiW)
            deallocate(allxiN)

        end if

        deallocate(buffer)
        call MPI_BUFFER_DETACH(commbuffer,commbuffer_size,ierr)

    end subroutine writeFile

end module RST_writeFile
