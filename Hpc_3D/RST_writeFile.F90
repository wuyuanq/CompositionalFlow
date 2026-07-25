
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
       
        character(len=50) :: fpwtxt, fuxtxt, fuytxt, fuztxt, fswtxt, fxitxt
        character(len=50), dimension(:), pointer :: fmftxt
        character :: charm
        real(kind=8), dimension(:,:,:), pointer :: allPw
        real(kind=8), dimension(:,:,:), pointer :: alluwx
        real(kind=8), dimension(:,:,:), pointer :: alluwy
        real(kind=8), dimension(:,:,:), pointer :: alluwz
        real(kind=8), dimension(:,:,:), pointer :: allunx
        real(kind=8), dimension(:,:,:), pointer :: alluny
        real(kind=8), dimension(:,:,:), pointer :: allunz
        real(kind=8), dimension(:,:,:,:), pointer :: allx
        real(kind=8), dimension(:,:,:), pointer :: allSw
        real(kind=8), dimension(:,:,:), pointer :: allxiW
        real(kind=8), dimension(:,:,:), pointer :: allxiN
        integer :: buffersize
        real(kind=8), dimension(:), pointer :: buffer
        real(kind=8), dimension(:), pointer :: Pwsent, uwxsent, uwysent, uwzsent, unxsent, unysent, &!
            unzsent, zsent, Swsent, xiWsent, xiNsent
        integer :: p_xlower, p_xupper, p_ylower, p_yupper, p_zlower, p_zupper
        integer :: p_pcol, p_prow, p_play
        integer :: position
        integer :: i, j, k, m, n, prociteration

        integer :: ierr
        integer :: request
        integer :: status(MPI_STATUS_SIZE)
        integer :: commbuffer_size = MAX_COMMBUF
        real(kind=8) :: commbuffer(MAX_COMMBUF)

        buffersize = (Nc+10)*local_size

        allocate(buffer(buffersize))
        call MPI_BUFFER_ATTACH(commbuffer,commbuffer_size,ierr)

        if(myid /= 0) then

            allocate(Pwsent(local_size))
            allocate(uwxsent(local_size))
            allocate(uwysent(local_size))
            allocate(uwzsent(local_size))
            allocate(unxsent(local_size))
            allocate(unysent(local_size))
            allocate(unzsent(local_size))
            allocate(zsent(local_size*Nc))
            allocate(Swsent(local_size))
            allocate(xiWsent(local_size))
            allocate(xiNsent(local_size))

            n = 1
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        Pwsent(n) = Pw(i,j,k)
                        n = n + 1
                    end do
                end do
            end do

            n = 1
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        uwxsent(n) = Uwx(i,j,k)
                        n = n + 1
                    end do
                end do
            end do

            n = 1
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        uwysent(n) = Uwy(i,j,k)
                        n = n + 1
                    end do
                end do
            end do

            n = 1
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        uwzsent(n) = Uwz(i,j,k)
                        n = n + 1
                    end do
                end do
            end do

            n = 1
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        unxsent(n) = Unx(i,j,k)
                        n = n + 1
                    end do
                end do
            end do

            n = 1
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        unysent(n) = Uny(i,j,k)
                        n = n + 1
                    end do
                end do
            end do

            n = 1
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        unzsent(n) = Unz(i,j,k)
                        n = n + 1
                    end do
                end do
            end do

            n = 1
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        do m = 1, Nc
                            zsent(n) = z(m,i,j,k)
                            n = n + 1
                        end do
                    end do
                end do
            end do

            n = 1
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        Swsent(n) = Sw(i,j,k)
                        n = n + 1
                    end do
                end do
            end do

            n = 1
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        xiWsent(n) = xiW(i,j,k)
                        n = n + 1
                    end do
                end do
            end do

            n = 1
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        xiNsent(n) = xiN(i,j,k)
                        n = n + 1
                    end do
                end do
            end do

            position = 0
            call MPI_PACK(Pwsent, local_size, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, MPI_COMM_WORLD, ierr)
            call MPI_PACK(uwxsent, local_size, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, MPI_COMM_WORLD, ierr)
            call MPI_PACK(uwysent, local_size, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, MPI_COMM_WORLD, ierr)
            call MPI_PACK(uwzsent, local_size, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, MPI_COMM_WORLD, ierr)
            call MPI_PACK(unxsent, local_size, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, MPI_COMM_WORLD, ierr)
            call MPI_PACK(unysent, local_size, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, MPI_COMM_WORLD, ierr)
            call MPI_PACK(unzsent, local_size, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, MPI_COMM_WORLD, ierr)
            call MPI_PACK(zsent, local_size*Nc, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, MPI_COMM_WORLD, ierr)
            call MPI_PACK(Swsent, local_size, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, MPI_COMM_WORLD, ierr)
            call MPI_PACK(xiWsent, local_size, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, MPI_COMM_WORLD, ierr)
            call MPI_PACK(xiNsent, local_size, MPI_DOUBLE_PRECISION, buffer, buffersize*8, position, MPI_COMM_WORLD, ierr)

            call MPI_IBSEND(buffer,buffersize,MPI_DOUBLE_PRECISION,0,myid+num_procs,MPI_COMM_WORLD,request,ierr)
            call MPI_WAIT(request, status, ierr)

            deallocate(Pwsent)
            deallocate(uwxsent)
            deallocate(uwysent)
            deallocate(uwzsent)
            deallocate(unxsent)
            deallocate(unysent)
            deallocate(unzsent)
            deallocate(zsent)
            deallocate(Swsent)
            deallocate(xiWsent)
            deallocate(xiNsent)

        else

            allocate(allPw(nx,ny,nz))
            allocate(alluwx(nx+1,ny,nz))
            allocate(alluwy(nx,ny+1,nz))
            allocate(alluwz(nx,ny,nz+1))
            allocate(allunx(nx+1,ny,nz))
            allocate(alluny(nx,ny+1,nz))
            allocate(allunz(nx,ny,nz+1))
            allocate(allx(Nc,nx,ny,nz))
            allocate(allSw(nx,ny,nz))
            allocate(allxiW(nx,ny,nz))
            allocate(allxiN(nx,ny,nz))

            alluwx(nx+1, 1:ny, 1:nz) = UwBdryX(2, 1:ny, 1:nz)
            allunx(nx+1, 1:ny, 1:nz) = UnBdryX(2, 1:ny, 1:nz)
            alluwy(1:nx, ny+1, 1:nz) = UwBdryY(1:nx, 2, 1:nz)
            alluny(1:nx, ny+1, 1:nz) = UnBdryY(1:nx, 2, 1:nz)
            alluwz(1:nx, 1:ny, nz+1) = UwBdryZ(1:nx, 1:ny, 2)
            allunz(1:nx, 1:ny, nz+1) = UnBdryZ(1:nx, 1:ny, 2)

            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        allPw(i,j,k) = Pw(i,j,k)
                        alluwx(i,j,k) = Uwx(i,j,k)
                        alluwy(i,j,k) = Uwy(i,j,k)
                        alluwz(i,j,k) = Uwz(i,j,k)
                        allunx(i,j,k) = Unx(i,j,k)
                        alluny(i,j,k) = Uny(i,j,k)
                        allunz(i,j,k) = Unz(i,j,k)
                        allSw(i,j,k) = Sw(i,j,k)
                        allxiW(i,j,k) = xiW(i,j,k)
                        allxiN(i,j,k) = xiN(i,j,k)
                    end do
                end do
            end do
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        do m = 1, Nc
                            allx(m,i,j,k) = z(m,i,j,k)
                        end do
                    end do
                end do
            end do

            do prociteration = 1, num_procs-1

                p_play = prociteration/(pnrows*pncols)+1
                p_prow = (prociteration-(p_play-1)*pnrows*pncols)/pncols+1
                p_pcol = (prociteration-(p_play-1)*pnrows*pncols)-(p_prow-1)*pncols+1

                p_xlower = (p_pcol-1)*localncols+1
                p_xupper = p_pcol*localncols
                p_ylower = (p_prow-1)*localnrows+1
                p_yupper = p_prow*localnrows
                p_zlower = (p_play-1)*localnlays+1
                p_zupper = p_play*localnlays

                call MPI_RECV(buffer, buffersize, MPI_DOUBLE_PRECISION, prociteration, &!
                    prociteration+num_procs, MPI_COMM_WORLD, status, ierr)

                n = 1
                do k = p_zlower, p_zupper
                    do j = p_ylower, p_yupper
                        do i = p_xlower, p_xupper
                            allPw(i,j,k) = buffer(n)
                            n = n + 1
                        end do
                    end do
                end do

                do k = p_zlower, p_zupper
                    do j = p_ylower, p_yupper
                        do i = p_xlower, p_xupper
                            alluwx(i,j,k) = buffer(n)
                            n = n + 1
                        end do
                    end do
                end do

                do k = p_zlower, p_zupper
                    do j = p_ylower, p_yupper
                        do i = p_xlower, p_xupper
                            alluwy(i,j,k) = buffer(n)
                            n = n + 1
                        end do
                    end do
                end do

                do k = p_zlower, p_zupper
                    do j = p_ylower, p_yupper
                        do i = p_xlower, p_xupper
                            alluwz(i,j,k) = buffer(n)
                            n = n + 1
                        end do
                    end do
                end do

                do k = p_zlower, p_zupper
                    do j = p_ylower, p_yupper
                        do i = p_xlower, p_xupper
                            allunx(i,j,k) = buffer(n)
                            n = n + 1
                        end do
                    end do
                end do

                do k = p_zlower, p_zupper
                    do j = p_ylower, p_yupper
                        do i = p_xlower, p_xupper
                            alluny(i,j,k) = buffer(n)
                            n = n + 1
                        end do
                    end do
                end do

                do k = p_zlower, p_zupper
                    do j = p_ylower, p_yupper
                        do i = p_xlower, p_xupper
                            allunz(i,j,k) = buffer(n)
                            n = n + 1
                        end do
                    end do
                end do

                do k = p_zlower, p_zupper
                    do j = p_ylower, p_yupper
                        do i = p_xlower, p_xupper
                            do m = 1, Nc
                                allx(m,i,j,k) = buffer(n)
                                n = n + 1
                            end do
                        end do
                    end do
                end do

                do k = p_zlower, p_zupper
                    do j = p_ylower, p_yupper
                        do i = p_xlower, p_xupper
                            allSw(i,j,k) = buffer(n)
                            n = n + 1
                        end do
                    end do
                end do

                do k = p_zlower, p_zupper
                    do j = p_ylower, p_yupper
                        do i = p_xlower, p_xupper
                            allxiW(i,j,k) = buffer(n)
                            n = n + 1
                        end do
                    end do
                end do

                do k = p_zlower, p_zupper
                    do j = p_ylower, p_yupper
                        do i = p_xlower, p_xupper
                            allxiN(i,j,k) = buffer(n)
                            n = n + 1
                        end do
                    end do
                end do

            end do

            fpwtxt = trim(adjustl(soludoc))//"/soln_cmp2PhFlw_Pw_raw.txt"
            fuxtxt = trim(adjustl(soludoc))//"/soln_cmp2PhFlw_Ux_raw.txt"
            fuytxt = trim(adjustl(soludoc))//"/soln_cmp2PhFlw_Uy_raw.txt"
            fuztxt = trim(adjustl(soludoc))//"/soln_cmp2PhFlw_Uz_raw.txt"
            fswtxt = trim(adjustl(soludoc))//"/soln_cmp2PhFlw_Sw_raw.txt"
            fxitxt = trim(adjustl(soludoc))//"/soln_cmp2PhFlw_c_raw.txt"
            
            allocate(fmftxt(Nc))
            do m = 1, Nc
                write(charm,'(i1)') m
                fmftxt(m) = trim(adjustl(soludoc))//"/soln_cmp2PhFlw_x"//charm//"_raw.txt"
            end do

            open(unit=10, file=trim(adjustl(fpwtxt)), status='replace')
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        write(10, fmt="(es15.8)") allPw(i,j,k)
                    end do
                end do
            end do
            close(10)

            open(unit=10, file=trim(adjustl(fuxtxt)), status='replace')
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx+1
                        write(10, fmt="(es15.8)") alluwx(i,j,k)+allunx(i,j,k)
                    end do
                end do
            end do
            close(10)

            open(unit=10, file=trim(adjustl(fuytxt)), status='replace')
            do k = 1, nz
                do j = 1, ny+1
                    do i = 1, nx
                        write(10, fmt="(es15.8)") alluwy(i,j,k)+alluny(i,j,k)
                    end do
                end do
            end do
            close(10)

            open(unit=10, file=trim(adjustl(fuztxt)), status='replace')
            do k = 1, nz+1
                do j = 1, ny
                    do i = 1, nx
                        write(10, fmt="(es15.8)") alluwz(i,j,k)+allunz(i,j,k)
                    end do
                end do
            end do
            close(10)

            do m = 1, Nc
                open(unit=10, file=trim(adjustl(fmftxt(m))), status='replace')
                do k = 1, nz
                    do j = 1, ny
                        do i = 1, nx
                            write(10, fmt="(es15.8)") allx(m,i,j,k)
                        end do
                    end do
                end do
                close(10)
            end do

            open(unit=10, file=trim(adjustl(fswtxt)), status='replace')
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        write(10, fmt="(es15.8)") allSw(i,j,k)
                    end do
                end do
            end do
            close(10)

            open(unit=10, file=trim(adjustl(fxitxt)), status='replace')
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        write(10, fmt="(es15.8)") allxiW(i,j,k)*allSw(i,j,k)+allxiN(i,j,k)*(1-allSw(i,j,k))
                    end do
                end do
            end do
            close(10)

            deallocate(fmftxt)

            deallocate(allPw)
            deallocate(alluwx)
            deallocate(alluwy)
            deallocate(alluwz)
            deallocate(allunx)
            deallocate(alluny)
            deallocate(allunz)
            deallocate(allx)
            deallocate(allSw)
            deallocate(allxiW)
            deallocate(allxiN)

        end if

        deallocate(buffer)
        call MPI_BUFFER_DETACH(commbuffer,commbuffer_size,ierr)

    end subroutine writeFile

end module RST_writeFile
