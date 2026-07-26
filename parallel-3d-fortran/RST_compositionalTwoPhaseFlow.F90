
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_compositionalTwoPhaseFlow

    use RST_model
    use RST_globalData
#ifdef SPARSE
    use RST_flashcalculation_fullgrid
#else
    use RST_flashcalculation
#endif
    use RST_viscosity
    include "mpif.h"

contains

    subroutine initialize(modelCase)

        implicit none

        type(model), intent(in out) :: modelCase
       
        character(len=50) :: fmhtxt, fmrtxt
        integer :: indexl, indexr, indexu, indexd, indexf, indexb
        integer :: ierr, errorcode
        integer :: i, j, k, m, n 

        call MPI_INIT(ierr)

        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        timestart = MPI_Wtime()
        flashtime = 0.0
        solvertime = 0.0
        commtime = 0.0

        pncols = modelCase%pncols
        pnrows = modelCase%pnrows
        pnlays = modelCase%pnlays
        Nc = modelCase%Nc
        Temp = modelCase%Temp
        Lx = modelCase%Lx
        Ly = modelCase%Ly
        Lz = modelCase%Lz
        timeEnd = modelCase%timeEnd
        nx = modelCase%nx
        ny = modelCase%ny
        nz = modelCase%nz
        nt = modelCase%nt
        gravX = modelCase%gravX
        gravY = modelCase%gravY
        gravZ = modelCase%gravZ

        allocate(xs(nx+1))
        allocate(ys(ny+1))
        allocate(zs(nz+1))
        allocate(ts(nt+1))
        allocate(Kxx(nx,ny,nz))
        allocate(Kyy(nx,ny,nz))
        allocate(Kzz(nx,ny,nz))
        allocate(poro(nx,ny,nz))
        allocate(src(Nc,nx,ny,nz))
        allocate(isDiriX(2,ny,nz))
        allocate(isDiriY(nx,2,nz))
        allocate(isDiriZ(nx,ny,2))
        allocate(PwBdryX(2,ny,nz))
        allocate(PwBdryY(nx,2,nz))
        allocate(PwBdryZ(nx,ny,2))
        allocate(PwInit(nx,ny,nz))
        allocate(zBdryX(Nc,2,ny,nz))
        allocate(zBdryY(Nc,nx,2,nz))
        allocate(zBdryZ(Nc,nx,ny,2))
        allocate(zInit(Nc,nx,ny,nz))
        allocate(UwBdryX(2,ny,nz))
        allocate(UwBdryY(nx,2,nz))
        allocate(UwBdryZ(nx,ny,2))
        allocate(UnBdryX(2,ny,nz))
        allocate(UnBdryY(nx,2,nz))
        allocate(UnBdryZ(nx,ny,2))
        allocate(ct(Nc))
        allocate(cp(Nc))
        allocate(af(Nc))
        allocate(mw(Nc))
        allocate(cv(Nc))
        allocate(psatA(Nc))
        allocate(psatB(Nc))
        allocate(psatC(Nc))
        allocate(delta(Nc,Nc))

        xs = modelCase%xs
        ys = modelCase%ys
        zs = modelCase%zs
        ts = modelCase%ts
        Kxx = modelCase%Kxx
        Kyy = modelCase%Kyy
        Kzz = modelCase%Kzz
        poro = modelCase%poro
        src = modelCase%src
        isDiriX = modelCase%isDiriX
        isDiriY = modelCase%isDiriY
        isDiriZ = modelCase%isDiriZ
        PwBdryX = modelCase%PwBdryX
        PwBdryY = modelCase%PwBdryY
        PwBdryZ = modelCase%PwBdryZ
        PwInit = modelCase%PwInit
        zBdryX = modelCase%zBdryX
        zBdryY = modelCase%zBdryY
        zBdryZ = modelCase%zBdryZ
        zInit = modelCase%zInit
        UwBdryX = modelCase%UwBdryX
        UwBdryY = modelCase%UwBdryY
        UwBdryZ = modelCase%UwBdryZ
        UnBdryX = modelCase%UnBdryX
        UnBdryY = modelCase%UnBdryY
        UnBdryZ = modelCase%UnBdryZ
        ct = modelCase%ct
        cp = modelCase%cp
        af = modelCase%af
        mw = modelCase%mw
        cv = modelCase%cv
        psatA = modelCase%psatA
        psatB = modelCase%psatB
        psatC = modelCase%psatC
        delta = modelCase%delta
        soludoc = modelCase%soludoc

        UwBdryX(1, 1:ny, 1:nz) = -UwBdryX(1, 1:ny, 1:nz)
        UwBdryY(1:nx, 1, 1:nz) = -UwBdryY(1:nx, 1, 1:nz)
        UwBdryZ(1:nx, 1:ny, 1) = -UwBdryZ(1:nx, 1:ny, 1)
        UnBdryX(1, 1:ny, 1:nz) = -UnBdryX(1, 1:ny, 1:nz)
        UnBdryY(1:nx, 1, 1:nz) = -UnBdryY(1:nx, 1, 1:nz)
        UnBdryZ(1:nx, 1:ny, 1) = -UnBdryZ(1:nx, 1:ny, 1)

        call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

        if(num_procs /= pnrows*pncols*pnlays) then
            print *, 'Please set the number of processes as ', pnrows*pncols*pnlays
            call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
        end if

        if(mod(nx*ny*nz, num_procs) /= 0) then
            print *, 'Please adjust the number of processes in order to make sure every'
            print *, 'process has the same number of cells.'
            call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
        end if

        localncols = nx/pncols
        localnrows = ny/pnrows
        localnlays = nz/pnlays

        play = myid/(pnrows*pncols)+1
        prow = (myid-(play-1)*pnrows*pncols)/pncols+1
        pcol = (myid-(play-1)*pnrows*pncols)-(prow-1)*pncols+1

        xlower = (pcol-1)*localncols+1
        xupper = pcol*localncols
        ylower = (prow-1)*localnrows+1
        yupper = prow*localnrows
        zlower = (play-1)*localnlays+1
        zupper = play*localnlays

        allocate(Pw(0:localncols+1, 0:localnrows+1, 0:localnlays+1))

        allocate(Uwx(localncols+1, localnrows, localnlays))
        Uwx = 0
        allocate(Uwy(localncols, localnrows+1, localnlays))
        Uwy = 0
        allocate(Uwz(localncols, localnrows, localnlays+1))
        Uwz = 0
        allocate(Unx(localncols+1, localnrows, localnlays))
        Unx = 0
        allocate(Uny(localncols, localnrows+1, localnlays))
        Uny = 0
        allocate(Unz(localncols, localnrows, localnlays+1))
        Unz = 0

        allocate(Sw(0:localncols+1, 0:localnrows+1, 0:localnlays+1))

        allocate(lambdawx(localncols+1, localnrows, localnlays))
        allocate(lambdawy(localncols, localnrows+1, localnlays))
        allocate(lambdawz(localncols, localnrows, localnlays+1))
        allocate(lambdanx(localncols+1, localnrows, localnlays))
        allocate(lambdany(localncols, localnrows+1, localnlays))
        allocate(lambdanz(localncols, localnrows, localnlays+1))

        allocate(Kxxbar(localncols+1, localnrows, localnlays))
        allocate(Kyybar(localncols, localnrows+1, localnlays))
        allocate(Kzzbar(localncols, localnrows, localnlays+1))

        allocate(z(Nc, 0:localncols+1, 0:localnrows+1, 0:localnlays+1))

        allocate(densiW(0:localncols+1, 0:localnrows+1, 0:localnlays+1))
        allocate(densiN(0:localncols+1, 0:localnrows+1, 0:localnlays+1))
        allocate(densiWbarx(localncols+1, localnrows, localnlays))
        allocate(densiWbary(localncols, localnrows+1, localnlays))
        allocate(densiWbarz(localncols, localnrows, localnlays+1))
        allocate(densiNbarx(localncols+1, localnrows, localnlays))
        allocate(densiNbary(localncols, localnrows+1, localnlays))
        allocate(densiNbarz(localncols, localnrows, localnlays+1))

        allocate(xW(Nc, 0:localncols+1, 0:localnrows+1, 0:localnlays+1))
        allocate(xN(Nc, 0:localncols+1, 0:localnrows+1, 0:localnlays+1))
        allocate(xWbarx(Nc, localncols+1, localnrows, localnlays))
        allocate(xWbary(Nc, localncols, localnrows+1, localnlays))
        allocate(xWbarz(Nc, localncols, localnrows, localnlays+1))
        allocate(xNbarx(Nc, localncols+1, localnrows, localnlays))
        allocate(xNbary(Nc, localncols, localnrows+1, localnlays))
        allocate(xNbarz(Nc, localncols, localnrows, localnlays+1))

        allocate(xiW(0:localncols+1, 0:localnrows+1, 0:localnlays+1))
        allocate(xiN(0:localncols+1, 0:localnrows+1, 0:localnlays+1))
        allocate(xiWbarx(localncols+1, localnrows, localnlays))
        allocate(xiWbary(localncols, localnrows+1, localnlays))
        allocate(xiWbarz(localncols, localnrows, localnlays+1))
        allocate(xiNbarx(localncols+1, localnrows, localnlays))
        allocate(xiNbary(localncols, localnrows+1, localnlays))
        allocate(xiNbarz(localncols, localnrows, localnlays+1))

        allocate(v(Nc, 0:localncols+1, 0:localnrows+1, 0:localnlays+1))

        allocate(Cf(0:localncols+1, 0:localnrows+1, 0:localnlays+1))

        allocate(viscW(0:localncols+1, 0:localnrows+1, 0:localnlays+1))
        allocate(viscN(0:localncols+1, 0:localnrows+1, 0:localnlays+1))

        if(pcol == 1) then
            Uwx(1, 1:localnrows, 1:localnlays) = UwBdryX(1, ylower:ylower+localnrows-1, zlower:zlower+localnlays-1)
            Unx(1, 1:localnrows, 1:localnlays) = UnBdryX(1, ylower:ylower+localnrows-1, zlower:zlower+localnlays-1)
        end if
        if(pcol == pncols) then
            Uwx(localncols+1, 1:localnrows, 1:localnlays) = UwBdryX(2, ylower:ylower+localnrows-1, zlower:zlower+localnlays-1)
            Unx(localncols+1, 1:localnrows, 1:localnlays) = UnBdryX(2, ylower:ylower+localnrows-1, zlower:zlower+localnlays-1)
        end if
        if(prow == 1) then
            Uwy(1:localncols, 1, 1:localnlays) = UwBdryY(xlower:xlower+localncols-1, 1, zlower:zlower+localnlays-1)
            Uny(1:localncols, 1, 1:localnlays) = UnBdryY(xlower:xlower+localncols-1, 1, zlower:zlower+localnlays-1)
        end if
        if(prow == pnrows) then
            Uwy(1:localncols, localnrows+1, 1:localnlays) = UwBdryY(xlower:xlower+localncols-1, 2, zlower:zlower+localnlays-1)
            Uny(1:localncols, localnrows+1, 1:localnlays) = UnBdryY(xlower:xlower+localncols-1, 2, zlower:zlower+localnlays-1)
        end if
        if(play == 1) then
            Uwz(1:localncols, 1:localnrows, 1) = UwBdryZ(xlower:xlower+localncols-1, ylower:ylower+localnrows-1, 1)
            Unz(1:localncols, 1:localnrows, 1) = UnBdryZ(xlower:xlower+localncols-1, ylower:ylower+localnrows-1, 1)
        end if
        if(play == pnlays) then
            Uwz(1:localncols, 1:localnrows, localnlays+1) = UwBdryZ(xlower:xlower+localncols-1, ylower:ylower+localnrows-1, 2)
            Unz(1:localncols, 1:localnrows, localnlays+1) = UnBdryZ(xlower:xlower+localncols-1, ylower:ylower+localnrows-1, 2)
        end if

        if(pncols == 1) then ! only one process column
            do k = 1, localnlays
                do j = 1, localnrows
                    Kxxbar(1,j,k) = Kxx(1,ylower+j-1,zlower+k-1)
                end do
            end do
            do k = 1, localnlays
                do j = 1, localnrows
                    Kxxbar(localncols+1,j,k) = Kxx(nx,ylower+j-1,zlower+k-1)
                end do
            end do
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 2, localncols
                        Kxxbar(i,j,k) = (xs(xlower+i)-xs(xlower+i-2)) / ((xs(xlower+i-1)-xs(xlower+i-2))/&!
                            Kxx(xlower+i-2,ylower+j-1,zlower+k-1)+(xs(xlower+i)-xs(xlower+i-1))&!
                            /Kxx(xlower+i-1,ylower+j-1,zlower+k-1))
                    end do
                end do
            end do
        elseif(pcol == 1) then
            do k = 1, localnlays
                do j = 1, localnrows
                    Kxxbar(1,j,k) = Kxx(1,ylower+j-1,zlower+k-1)
                end do
            end do
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 2, localncols+1
                        Kxxbar(i,j,k) = (xs(xlower+i)-xs(xlower+i-2)) / ((xs(xlower+i-1)-xs(xlower+i-2))/&!
                            Kxx(xlower+i-2,ylower+j-1,zlower+k-1)+(xs(xlower+i)-xs(xlower+i-1))&!
                            /Kxx(xlower+i-1,ylower+j-1,zlower+k-1))
                    end do
                end do
            end do
        elseif(pcol == pncols) then
            do k = 1, localnlays
                do j = 1, localnrows
                    Kxxbar(localncols+1,j,k) = Kxx(nx,ylower+j-1,zlower+k-1)
                end do
            end do
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        Kxxbar(i,j,k) = (xs(xlower+i)-xs(xlower+i-2)) / ((xs(xlower+i-1)-xs(xlower+i-2))/&!
                            Kxx(xlower+i-2,ylower+j-1,zlower+k-1)+(xs(xlower+i)-xs(xlower+i-1)) &!
                            /Kxx(xlower+i-1,ylower+j-1,zlower+k-1))
                    end do
                end do
            end do
        else
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols+1
                        Kxxbar(i,j,k) = (xs(xlower+i)-xs(xlower+i-2)) / ((xs(xlower+i-1)-xs(xlower+i-2))/&!
                            Kxx(xlower+i-2,ylower+j-1,zlower+k-1)+(xs(xlower+i)-xs(xlower+i-1)) &!
                            /Kxx(xlower+i-1,ylower+j-1,zlower+k-1))
                    end do
                end do
            end do
        end if

        if(pnrows == 1) then ! only one process row
            do k = 1, localnlays
                do i = 1, localncols
                    Kyybar(i,1,k) = Kyy(xlower+i-1,1,zlower+k-1)
                end do
            end do
            do k = 1, localnlays
                do i = 1, localncols
                    Kyybar(i,localnrows+1,k) = Kyy(xlower+i-1,ny,zlower+k-1)
                end do
            end do
            do k = 1, localnlays
                do j = 2, localnrows
                    do i = 1, localncols
                        Kyybar(i,j,k) = (ys(ylower+j)-ys(ylower+j-2)) / ((ys(ylower+j-1)-ys(ylower+j-2))/&!
                            Kyy(xlower+i-1,ylower+j-2,zlower+k-1)+(ys(ylower+j)-ys(ylower+j-1)) &!
                            /Kyy(xlower+i-1,ylower+j-1,zlower+k-1))
                    end do
                end do
            end do
        elseif(prow == 1) then
            do k = 1, localnlays
                do i = 1, localncols
                    Kyybar(i,1,k) = Kyy(xlower+i-1,1,zlower+k-1)
                end do
            end do
            do k = 1, localnlays
                do j = 2, localnrows+1
                    do i = 1, localncols
                        Kyybar(i,j,k) = (ys(ylower+j)-ys(ylower+j-2)) / ((ys(ylower+j-1)-ys(ylower+j-2))/&!
                            Kyy(xlower+i-1,ylower+j-2,zlower+k-1)+(ys(ylower+j)-ys(ylower+j-1)) &!
                            /Kyy(xlower+i-1,ylower+j-1,zlower+k-1))
                    end do
                end do
            end do
        elseif(prow == pnrows) then
            do k = 1, localnlays
                do i = 1, localncols
                    Kyybar(i,localnrows+1,k) = Kyy(xlower+i-1, ny,zlower+k-1)
                end do
            end do
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        Kyybar(i,j,k) = (ys(ylower+j)-ys(ylower+j-2)) / ((ys(ylower+j-1)-ys(ylower+j-2))/&!
                            Kyy(xlower+i-1,ylower+j-2,zlower+k-1)+(ys(ylower+j)-ys(ylower+j-1)) &!
                            /Kyy(xlower+i-1,ylower+j-1,zlower+k-1))
                    end do
                end do
            end do
        else
            do k = 1, localnlays
                do j = 1, localnrows+1
                    do i = 1, localncols
                        Kyybar(i,j,k) = (ys(ylower+j)-ys(ylower+j-2)) / ((ys(ylower+j-1)-ys(ylower+j-2))/&!
                            Kyy(xlower+i-1,ylower+j-2,zlower+k-1)+(ys(ylower+j)-ys(ylower+j-1)) &!
                            /Kyy(xlower+i-1,ylower+j-1,zlower+k-1))
                    end do
                end do
            end do
        end if

        if(pnlays == 1) then ! only one process lay
            do j = 1, localnrows
                do i = 1, localncols
                    Kzzbar(i,j,1) = Kzz(xlower+i-1,ylower+j-1,1)
                end do
            end do
            do j = 1, localnrows
                do i = 1, localncols
                    Kzzbar(i,j,localnlays+1) = Kzz(xlower+i-1,ylower+j-1,nz)
                end do
            end do
            do k = 2, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        Kzzbar(i,j,k) = (zs(zlower+k)-zs(zlower+k-2)) / ((zs(zlower+k-1)-zs(zlower+k-2))/&!
                            Kzz(xlower+i-1,ylower+j-1,zlower+k-2)+(zs(zlower+k)-zs(zlower+k-1)) &!
                            /Kzz(xlower+i-1,ylower+j-1,zlower+k-1))
                    end do
                end do
            end do
        elseif(play == 1) then
            do j = 1, localnrows
                do i = 1, localncols
                    Kzzbar(i,j,1) = Kzz(xlower+i-1,ylower+j-1,1)
                end do
            end do
            do k = 2, localnlays+1
                do j = 1, localnrows
                    do i = 1, localncols
                        Kzzbar(i,j,k) = (zs(zlower+k)-zs(zlower+k-2)) / ((zs(zlower+k-1)-zs(zlower+k-2))/&!
                            Kzz(xlower+i-1,ylower+j-1,zlower+k-2)+(zs(zlower+k)-zs(zlower+k-1)) &!
                            /Kzz(xlower+i-1,ylower+j-1,zlower+k-1))
                    end do
                end do
            end do
        elseif(play == pnlays) then
            do j = 1, localnrows
                do i = 1, localncols
                    Kzzbar(i,j,localnlays+1) = Kzz(xlower+i-1,ylower+j-1,nz)
                end do
            end do
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        Kzzbar(i,j,k) = (zs(zlower+k)-zs(zlower+k-2)) / ((zs(zlower+k-1)-zs(zlower+k-2))/&!
                            Kzz(xlower+i-1,ylower+j-1,zlower+k-2)+(zs(zlower+k)-zs(zlower+k-1)) &!
                            /Kzz(xlower+i-1,ylower+j-1,zlower+k-1))
                    end do
                end do
            end do
        else
            do k = 1, localnlays+1
                do j = 1, localnrows
                    do i = 1, localncols
                        Kzzbar(i,j,k) = (zs(zlower+k)-zs(zlower+k-2)) / ((zs(zlower+k-1)-zs(zlower+k-2))/&!
                            Kzz(xlower+i-1,ylower+j-1,zlower+k-2)+(zs(zlower+k)-zs(zlower+k-1)) &!
                            /Kzz(xlower+i-1,ylower+j-1,zlower+k-1))
                    end do
                end do
            end do
        end if

        indexl = 0
        indexr = localncols+1
        indexd = 0
        indexu = localnrows+1
        indexf = 0
        indexb = localnlays+1

        if(pcol == 1) then
            indexl = 1
        end if
        if(pcol == pncols) then
            indexr = localncols
        end if
        if(prow == 1) then
            indexd = 1
        end if
        if(prow == pnrows) then
            indexu = localnrows
        end if
        if(play == 1) then
            indexf = 1
        end if
        if(play == pnlays) then
            indexb = localnlays
        end if

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr
                    Pw(i,j,k) = PwInit(xlower+i-1,ylower+j-1,zlower+k-1)
                end do
            end do
        end do
        if(pcol==1) then
            do k = 1, localnlays
                do j = 1, localnrows
                    Pw(0,j,k) = PwBdryX(1,ylower+j-1,zlower+k-1)
                end do
            end do
        end if
        if(pcol==pncols) then
            do k = 1, localnlays
                do j = 1, localnrows
                    Pw(localncols+1,j,k) = PwBdryX(2,ylower+j-1,zlower+k-1)
                end do
            end do
        end if
        if(prow==1) then
            do k = 1, localnlays
                do i = 1, localncols
                    Pw(i,0,k) = PwBdryY(xlower+i-1,1,zlower+k-1)
                end do
            end do
        end if
        if(prow==pnrows) then
            do k = 1, localnlays
                do i = 1, localncols
                    Pw(i,localnrows+1,k) = PwBdryY(xlower+i-1,2,zlower+k-1)
                end do
            end do
        end if
        if(play==1) then
            do j = 1, localnrows
                do i = 1, localncols
                    Pw(i,j,0) = PwBdryZ(xlower+i-1,ylower+j-1,1)
                end do
            end do
        end if
        if(play==pnlays) then
            do j = 1, localnrows
                do i = 1, localncols
                    Pw(i,j,localnlays+1) = PwBdryZ(xlower+i-1,ylower+j-1,2)
                end do
            end do
        end if

        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr
                    do m = 1, Nc
                        z(m,i,j,k) = zInit(m,xlower+i-1,ylower+j-1,zlower+k-1)
                    end do
                end do
            end do
        end do

#ifdef SPARSE
        allocate(xtable(1:Nc,1:TABLESIZE))
        allocate(ytable(1:Nc,1:TABLESIZE))
        allocate(xiLtable(1:TABLESIZE))
        allocate(xiGtable(1:TABLESIZE))
        allocate(rhoLtable(1:TABLESIZE))
        allocate(rhoGtable(1:TABLESIZE))
        allocate(sLtable(1:TABLESIZE))
        allocate(vtable(1:Nc,1:TABLESIZE))
        allocate(Cftable(1:TABLESIZE))
        allocate(isWtable(1:TABLESIZE))
        allocate(isNtable(1:TABLESIZE))

        open(unit=120, file=FULLGRIDPREFIX//"xW1.txt", status='old')
        open(unit=130, file=FULLGRIDPREFIX//"xW2.txt", status='old')
        open(unit=150, file=FULLGRIDPREFIX//"xN1.txt", status='old')
        open(unit=160, file=FULLGRIDPREFIX//"xN2.txt", status='old')
        open(unit=180, file=FULLGRIDPREFIX//"xiW.txt", status='old')
        open(unit=190, file=FULLGRIDPREFIX//"xiN.txt", status='old')
        open(unit=200, file=FULLGRIDPREFIX//"densiW.txt", status='old')
        open(unit=210, file=FULLGRIDPREFIX//"densiN.txt", status='old')
        open(unit=220, file=FULLGRIDPREFIX//"sW.txt", status='old')
        open(unit=230, file=FULLGRIDPREFIX//"v1.txt", status='old')
        open(unit=240, file=FULLGRIDPREFIX//"v2.txt", status='old')
        open(unit=260, file=FULLGRIDPREFIX//"Cf.txt", status='old')
        open(unit=270, file=FULLGRIDPREFIX//"liquid.txt", status='old')
        open(unit=280, file=FULLGRIDPREFIX//"gas.txt", status='old')
        if(Nc == 3) then
            open(unit=140, file=FULLGRIDPREFIX//"xW3.txt", status='old')
            open(unit=170, file=FULLGRIDPREFIX//"xN3.txt", status='old')
            open(unit=250, file=FULLGRIDPREFIX//"v3.txt", status='old')
        end if

        read(120,*) xtable(1,1:TABLESIZE)
        read(130,*) xtable(2,1:TABLESIZE)
        read(150,*) ytable(1,1:TABLESIZE)
        read(160,*) ytable(2,1:TABLESIZE)
        read(180,*) xiLtable(1:TABLESIZE)
        read(190,*) xiGtable(1:TABLESIZE)
        read(200,*) rhoLtable(1:TABLESIZE)
        read(210,*) rhoGtable(1:TABLESIZE)
        read(220,*) sLtable(1:TABLESIZE)
        read(230,*) vtable(1,1:TABLESIZE)
        read(240,*) vtable(2,1:TABLESIZE)
        read(260,*) Cftable(1:TABLESIZE)
        read(270,*) isWtable(1:TABLESIZE)
        read(280,*) isNtable(1:TABLESIZE)
        if(Nc == 3) then
            read(140,*) xtable(3,1:TABLESIZE) 
            read(170,*) ytable(3,1:TABLESIZE) 
            read(250,*) vtable(3,1:TABLESIZE) 
        end if
#endif

        totalmole = 0.0
        t = 2

        local_size = nx*ny*nz/num_procs
        allocate(initial_x_guess(local_size))
        n = 1
        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols
                    initial_x_guess(n) = Pw(i,j,k)
                    n = n + 1
                end do
            end do
        end do

        stencil_indices(1:7) = (/0, 1, 2, 3, 4, 5, 6/)
        offsets(1:7,1:3) = reshape([0,-1,1,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,-1,1], [7,3])

        ilower(1) = xlower
        iupper(1) = xupper
        ilower(2) = ylower
        iupper(2) = yupper
        ilower(3) = zlower
        iupper(3) = zupper

        call HYPRE_StructGridCreate(MPI_COMM_WORLD, 3, grid, ierr)
        call HYPRE_StructGridSetExtents(grid, ilower, iupper, ierr)
        call HYPRE_StructGridAssemble(grid, ierr)

        call HYPRE_StructStencilCreate(3, 7, stencil, ierr)
        do n = 1, 7
            call HYPRE_StructStencilSetElement(stencil, n-1, offsets(n,1:3), ierr)
        end do

        call HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, global_A, ierr)
        call HYPRE_StructMatrixInitialize(global_A, ierr)

        call HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, global_b, ierr)
        call HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, global_x, ierr)
        call HYPRE_StructVectorInitialize(global_b, ierr)
        call HYPRE_StructVectorInitialize(global_x, ierr)

        call HYPRE_StructSMGCreate(MPI_COMM_WORLD, solver, ierr)
        call HYPRE_StructSMGSetMemoryUse(solver, 0, ierr)
        call HYPRE_StructSMGSetTol(solver, 1.0e-07, ierr)
        call HYPRE_StructSMGSetRelChange(solver, 0, ierr)
        call HYPRE_StructSMGSetNumPreRelax(solver, 1, ierr)
        call HYPRE_StructSMGSetNumPostRelax(solver, 1, ierr)

        if(myid == 0) then
            fmhtxt = trim(adjustl(soludoc))//"/soln_cmp2PhFlw_moleHistory.txt"
            fmrtxt = trim(adjustl(soludoc))//"/soln_cmp2PhFlw_moleRatioHistory.txt"

            open(unit=20, file=trim(adjustl(fmhtxt)), status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file error. ', ierr
                call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
            end if
            open(unit=30, file=trim(adjustl(fmrtxt)), status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file error. ', ierr
                call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
            end if
        end if

        deallocate(modelCase%xs)
        deallocate(modelCase%ys)
        deallocate(modelCase%zs)
        deallocate(modelCase%ts)
        deallocate(modelCase%Kxx)
        deallocate(modelCase%Kyy)
        deallocate(modelCase%Kzz)
        deallocate(modelCase%poro)
        deallocate(modelCase%src)
        deallocate(modelCase%isDiriX)
        deallocate(modelCase%isDiriY)
        deallocate(modelCase%isDiriZ)
        deallocate(modelCase%PwBdryX)
        deallocate(modelCase%PwBdryY)
        deallocate(modelCase%PwBdryZ)
        deallocate(modelCase%PwInit)
        deallocate(modelCase%zBdryX)
        deallocate(modelCase%zBdryY)
        deallocate(modelCase%zBdryZ)
        deallocate(modelCase%zInit)
        deallocate(modelCase%UwBdryX)
        deallocate(modelCase%UwBdryY)
        deallocate(modelCase%UwBdryZ)
        deallocate(modelCase%UnBdryX)
        deallocate(modelCase%UnBdryY)
        deallocate(modelCase%UnBdryZ)
        deallocate(modelCase%ct)
        deallocate(modelCase%cp)
        deallocate(modelCase%af)
        deallocate(modelCase%mw)
        deallocate(modelCase%cv)
        deallocate(modelCase%delta)

    end subroutine initialize

    subroutine computeParameters()

        implicit none
        real(kind=8), dimension(:), pointer :: xWtemp, xNtemp, ztemp, vtemp
        real(kind=8) :: ignore1, Swtemp, viscWtemp, viscNtemp
        logical :: isW, isN, isRea
        real(kind=8), dimension(:,:,:,:), pointer :: moleincell
        real(kind=8), dimension(:), pointer :: leftmole
        real(kind=8) :: totaldesiredleftmole, sum, leftmole1sum, sentbuffer(2)
        integer :: indexl, indexr, indexd, indexu, indexf, indexb
        integer :: i, j, k, m, prociteration

        integer :: ierr, errorcode
        integer :: request, position
        integer :: status(MPI_STATUS_SIZE)
        integer :: buffer_size = MAX_BUF
        real(kind=8) :: buffer(MAX_BUF)
        real(kind=8) :: flashtimestart, flashtimefinish

        allocate(ztemp(Nc))
        allocate(xWtemp(Nc))
        allocate(xNtemp(Nc))
        allocate(vtemp(Nc))

        indexl = 0
        indexr = localncols+1
        indexd = 0
        indexu = localnrows+1
        indexf = 0
        indexb = localnlays+1

        if(pcol == 1) then
            indexl = 1
        end if
        if(pcol == pncols) then
            indexr = localncols
        end if
        if(prow == 1) then
            indexd = 1
        end if
        if(prow == pnrows) then
            indexu = localnrows
        end if
        if(play == 1) then
            indexf = 1
        end if
        if(play == pnlays) then
            indexb = localnlays
        end if

        ! the loop order should be the first index of the array be the innerest loop, and the third index should be the most outside.
        ! Such order makes the read of the array faster.
        do k = indexf, indexb
            do j = indexd, indexu
                do i = indexl, indexr

                    if((.not.((i==0).and.(j==0))).and. &!
                        (.not.((i==0).and.(j==localnrows+1))).and. &!
                        (.not.((i==localncols+1).and.(j==0))).and. &!
                        (.not.((i==localncols+1).and.(j==localnrows+1))).and. &!
                        (.not.((i==0).and.(k==0))).and. &!
                        (.not.((i==0).and.(k==localnlays+1))).and. &!
                        (.not.((i==localncols+1).and.(k==0))).and. &!
                        (.not.((i==localncols+1).and.(k==localnlays+1))).and. &!
                        (.not.((j==0).and.(k==0))).and. &!
                        (.not.((j==0).and.(k==localnlays+1))).and. &!
                        (.not.((j==localnrows+1).and.(k==0))).and. &!
                        (.not.((j==localnrows+1).and.(k==localnlays+1)))) then

                        ztemp(1:Nc) = z(1:Nc,i,j,k)

                        flashtimestart = MPI_Wtime()

#ifdef SPARSE
                        call flashcalculation_fullgrid( Pw(i,j,k), ztemp, xWtemp, xNtemp, xiW(i,j,k), xiN(i,j,k), &!
                            densiW(i,j,k), densiN(i,j,k), Sw(i,j,k), vtemp, Cf(i,j,k), isW, isN )
#else
                        call flashcalculation( Pw(i,j,k), ztemp, &!
                            xWtemp, xNtemp, xiW(i,j,k), xiN(i,j,k), densiW(i,j,k), densiN(i,j,k), &!
                            Sw(i,j,k), vtemp, Cf(i,j,k), isW, isN, isRea )
#endif

                        flashtimefinish = MPI_Wtime()
                        flashtime = flashtime + flashtimefinish - flashtimestart

                        do m = 1, Nc
                            xW(m,i,j,k) = xWtemp(m)
                            xN(m,i,j,k) = xNtemp(m)
                            v(m,i,j,k) = vtemp(m)
                        end do

                        if(isW) then
                            viscW(i,j,k) = viscosity( xWtemp, xiW(i,j,k), Pw(i,j,k), 'l' )
                        else
                            viscW(i,j,k) = 1.D12
                        end if
                
                        if(isN) then
                            viscN(i,j,k) = viscosity( xNtemp, xiN(i,j,k), Pw(i,j,k), 'g' )
                        else
                            viscN(i,j,k) = 1.D12
                        end if
                    end if

                end do
            end do
        end do

        allocate(moleincell(Nc, localncols, localnrows, localnlays))
        
        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols
                    do m = 1, Nc
                        moleincell(m,i,j,k) = z(m,i,j,k)*(xiW(i,j,k)*Sw(i,j,k) + xiN(i,j,k)*(1-Sw(i,j,k)))&!
                            *(xs(xlower+i)-xs(xlower+i-1))*(ys(ylower+j)-ys(ylower+j-1))*(zs(zlower+k)-zs(zlower+k-1)) &!
                            *poro(xlower+i-1,ylower+j-1,zlower+k-1)
                    end do
                end do
            end do
        end do

        allocate(leftmole(Nc))
        leftmole = 0
        do m = 1, Nc
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        leftmole(m) = leftmole(m) + moleincell(m,i,j,k)
                    end do
                end do
            end do
        end do

        totaldesiredleftmole = 0.0
        do m = 2, Nc
            totaldesiredleftmole = totaldesiredleftmole + leftmole(m)
        end do

        call MPI_BUFFER_ATTACH(buffer,buffer_size,ierr)

        if(myid /= 0) then

            position = 0
            call MPI_PACK(totaldesiredleftmole, 1, MPI_DOUBLE_PRECISION, sentbuffer, 2*8, position, MPI_COMM_WORLD, ierr)
            call MPI_PACK(leftmole(1), 1, MPI_DOUBLE_PRECISION, sentbuffer, 2*8, position, MPI_COMM_WORLD, ierr)
            call MPI_IBSEND(sentbuffer,2,MPI_DOUBLE_PRECISION,0,myid+num_procs,MPI_COMM_WORLD,request,ierr)
            call MPI_WAIT(request, status, ierr)

        else

            sum = totaldesiredleftmole
            leftmole1sum = leftmole(1)

            do prociteration = 1, num_procs-1

                call MPI_RECV(sentbuffer, 2, MPI_DOUBLE_PRECISION, prociteration, prociteration+num_procs, &!
                    MPI_COMM_WORLD, status, ierr)

                sum = sum + sentbuffer(1)
                leftmole1sum = leftmole1sum + sentbuffer(2)

            end do

            if(t == 2) then
                totalmole = sum
            end if

            write(20, fmt="(es12.5)") (totalmole-sum)/totalmole
print *, (totalmole-sum)/totalmole
            write(30, fmt="(es12.5)") sum/leftmole1sum

        end if

        deallocate(leftmole)
        deallocate(moleincell)

        call MPI_BUFFER_DETACH(buffer,buffer_size,ierr)

        if(pncols == 1) then
            do k = 1, localnlays
                do j = 1, localnrows
                    if((Uwx(1,j,k) > 0).or.(Unx(1,j,k) > 0)) then
                        ztemp(1:Nc) = zBdryX(1:Nc,1,ylower+j-1,zlower+k-1)
#ifdef SPARSE
                        call flashcalculation_fullgrid( Pw(1,j,k), ztemp, xWtemp, xNtemp, xiWbarx(1,j,k), xiNbarx(1,j,k), &!
                            densiWbarx(1,j,k), densiNbarx(1,j,k), Swtemp, vtemp, ignore1, isW, isN )
#else
                        call flashcalculation(Pw(1,j,k), ztemp, xWtemp, xNtemp, xiWbarx(1,j,k), &!
                            xiNbarx(1,j,k), densiWbarx(1,j,k), densiNbarx(1,j,k), Swtemp, vtemp, ignore1, isW, isN, isRea)
#endif
                        xWbarx(1:Nc,1,j,k) = xWtemp(1:Nc)
                        xNbarx(1:Nc,1,j,k) = xNtemp(1:Nc)
                        if(isW) then
                            viscWtemp = viscosity( xWtemp, xiWbarx(1,j,k), Pw(1,j,k), 'l' )
                        end if
                        if(isN) then
                            viscNtemp = viscosity( xNtemp, xiNbarx(1,j,k), Pw(1,j,k), 'g' )
                        end if
                    end if
                    if(Uwx(1,j,k) > 0) then
                        lambdawx(1,j,k) = Kxxbar(1,j,k)*computekr_W(Swtemp)/viscWtemp
                    elseif(Uwx(1,j,k) < 0) then
                        lambdawx(1,j,k) = Kxxbar(1,j,k)*computekr_W(Sw(1,j,k))/viscW(1,j,k)
                        xiWbarx(1,j,k) = xiW(1,j,k)
                        densiWbarx(1,j,k) = densiW(1,j,k)
                        xWbarx(1:Nc,1,j,k) = xW(1:Nc,1,j,k)
                    elseif(isDiriX(1,j,k) == 1) then
                        lambdawx(1,j,k) = Kxxbar(1,j,k)*computekr_W(Sw(1,j,k))/viscW(1,j,k)
                        xiWbarx(1,j,k) = xiW(1,j,k)
                        densiWbarx(1,j,k) = densiW(1,j,k)
                        xWbarx(1:Nc,1,j,k) = xW(1:Nc,1,j,k)
                    end if
                    if(Unx(1,j,k) > 0) then
                        lambdanx(1,j,k) = Kxxbar(1,j,k)*computekr_N(Swtemp)/viscNtemp
                    elseif(Unx(1,j,k) < 0) then
                        lambdanx(1,j,k) = Kxxbar(1,j,k)*computekr_N(Sw(1,j,k))/viscN(1,j,k)
                        xiNbarx(1,j,k) = xiN(1,j,k)
                        densiNbarx(1,j,k) = densiN(1,j,k)
                        xNbarx(1:Nc,1,j,k) = xN(1:Nc,1,j,k)
                    elseif(isDiriX(1,j,k) == 1) then
                        lambdanx(1,j,k) = Kxxbar(1,j,k)*computekr_N(Sw(1,j,k))/viscN(1,j,k)
                        xiNbarx(1,j,k) = xiN(1,j,k)
                        densiNbarx(1,j,k) = densiN(1,j,k)
                        xNbarx(1:Nc,1,j,k) = xN(1:Nc,1,j,k)
                    end if
                end do
            end do
            do k = 1, localnlays
                do j = 1, localnrows
                    if((Uwx(localncols+1,j,k) < 0).or.(Unx(localncols+1,j,k) < 0)) then
                        ztemp(1:Nc) = zBdryX(1:Nc,2,ylower+j-1,zlower+k-1)
#ifdef SPARSE
                        call flashcalculation_fullgrid( Pw(localncols,j,k), ztemp, xWtemp, &!
                            xNtemp, xiWbarx(localncols+1,j,k), xiNbarx(localncols+1,j,k), densiWbarx(localncols+1,j,k), &!
                            densiNbarx(localncols+1,j,k), Swtemp, vtemp, ignore1, isW, isN )
#else
                        call flashcalculation(Pw(localncols,j,k), ztemp, xWtemp, &!
                            xNtemp, xiWbarx(localncols+1,j,k), xiNbarx(localncols+1,j,k), densiWbarx(localncols+1,j,k), &!
                            densiNbarx(localncols+1,j,k), Swtemp, vtemp, ignore1, isW, isN, isRea)
#endif
                        xWbarx(1:Nc,localncols+1,j,k) = xWtemp(1:Nc)
                        xNbarx(1:Nc,localncols+1,j,k) = xNtemp(1:Nc)
                        if(isW) then
                            viscWtemp = viscosity( xWtemp, xiWbarx(localncols+1,j,k), Pw(localncols,j,k), 'l' )
                        end if
                        if(isN) then
                            viscNtemp = viscosity( xNtemp, xiNbarx(localncols+1,j,k), Pw(localncols,j,k), 'g' )
                        end if
                    end if
                    if(Uwx(localncols+1,j,k) < 0) then
                        lambdawx(localncols+1,j,k) = Kxxbar(localncols+1,j,k)*computekr_W(Swtemp)/viscWtemp
                    elseif(Uwx(localncols+1,j,k) > 0) then
                        lambdawx(localncols+1,j,k) = Kxxbar(localncols+1,j,k)*computekr_W(Sw(localncols,j,k))/viscW(localncols,j,k)
                        xiWbarx(localncols+1,j,k) = xiW(localncols,j,k)
                        densiWbarx(localncols+1,j,k) = densiW(localncols,j,k)
                        xWbarx(1:Nc,localncols+1,j,k) = xW(1:Nc,localncols,j,k)
                    elseif(isDiriX(2,ylower+j-1,zlower+k-1) == 1) then
                        lambdawx(localncols+1,j,k) = Kxxbar(localncols+1,j,k)*computekr_W(Sw(localncols,j,k))/viscW(localncols,j,k)
                        xiWbarx(localncols+1,j,k) = xiW(localncols,j,k)
                        densiWbarx(localncols+1,j,k) = densiW(localncols,j,k)
                        xWbarx(1:Nc,localncols+1,j,k) = xW(1:Nc,localncols,j,k)
                    end if
                    if(Unx(localncols+1,j,k) < 0) then
                        lambdanx(localncols+1,j,k) = Kxxbar(localncols+1,j,k)*computekr_N(Swtemp)/viscNtemp
                    elseif(Unx(localncols+1,j,k) > 0) then
                        lambdanx(localncols+1,j,k) = Kxxbar(localncols+1,j,k)*computekr_N(Sw(localncols,j,k))/viscN(localncols,j,k)
                        xiNbarx(localncols+1,j,k) = xiN(localncols,j,k)
                        densiNbarx(localncols+1,j,k) = densiN(localncols,j,k)
                        xNbarx(1:Nc,localncols+1,j,k) = xN(1:Nc,localncols,j,k)
                    elseif(isDiriX(2,ylower+j-1,zlower+k-1) == 1) then
                        lambdanx(localncols+1,j,k) = Kxxbar(localncols+1,j,k)*computekr_N(Sw(localncols,j,k))/viscN(localncols,j,k)
                        xiNbarx(localncols+1,j,k) = xiN(localncols,j,k)
                        densiNbarx(localncols+1,j,k) = densiN(localncols,j,k)
                        xNbarx(1:Nc,localncols+1,j,k) = xN(1:Nc,localncols,j,k)
                    end if
                end do
            end do
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 2, localncols
                        if(Uwx(i,j,k) > 0) then
                            lambdawx(i,j,k) = Kxxbar(i,j,k)*computekr_W(Sw(i-1,j,k))/viscW(i-1,j,k)
                            xiWbarx(i,j,k) = xiW(i-1,j,k)
                            densiWbarx(i,j,k) = densiW(i-1,j,k)
                            xWbarx(1:Nc,i,j,k) = xW(1:Nc,i-1,j,k)
                        else
                            lambdawx(i,j,k) = Kxxbar(i,j,k)*computekr_W(Sw(i,j,k))/viscW(i,j,k)
                            xiWbarx(i,j,k) = xiW(i,j,k)
                            densiWbarx(i,j,k) = densiW(i,j,k)
                            xWbarx(1:Nc,i,j,k) = xW(1:Nc,i,j,k)
                        end if
                        if(Unx(i,j,k) > 0) then
                            lambdanx(i,j,k) = Kxxbar(i,j,k)*computekr_N(Sw(i-1,j,k))/viscN(i-1,j,k)
                            xiNbarx(i,j,k) = xiN(i-1,j,k)
                            densiNbarx(i,j,k) = densiN(i-1,j,k)
                            xNbarx(1:Nc,i,j,k) = xN(1:Nc,i-1,j,k)
                        else
                            lambdanx(i,j,k) = Kxxbar(i,j,k)*computekr_N(Sw(i,j,k))/viscN(i,j,k)
                            xiNbarx(i,j,k) = xiN(i,j,k)
                            densiNbarx(i,j,k) = densiN(i,j,k)
                            xNbarx(1:Nc,i,j,k) = xN(1:Nc,i,j,k)
                        end if
                    end do
                end do
            end do
        elseif(pcol==1) then
            do k = 1, localnlays
                do j = 1, localnrows
                    if((Uwx(1,j,k) > 0).or.(Unx(1,j,k) > 0)) then
                        ztemp(1:Nc) = zBdryX(1:Nc,1,ylower+j-1,zlower+k-1)
#ifdef SPARSE
                        call flashcalculation_fullgrid( Pw(1,j,k), ztemp, xWtemp, xNtemp, xiWbarx(1,j,k), xiNbarx(1,j,k), &!
                            densiWbarx(1,j,k), densiNbarx(1,j,k), Swtemp, vtemp, ignore1, isW, isN )
#else
                        call flashcalculation(Pw(1,j,k), ztemp, xWtemp, xNtemp, xiWbarx(1,j,k), &!
                            xiNbarx(1,j,k), densiWbarx(1,j,k), densiNbarx(1,j,k), Swtemp, vtemp, ignore1, isW, isN, isRea)
#endif
                        xWbarx(1:Nc,1,j,k) = xWtemp(1:Nc)
                        xNbarx(1:Nc,1,j,k) = xNtemp(1:Nc)
                        if(isW) then
                            viscWtemp = viscosity( xWtemp, xiWbarx(1,j,k), Pw(1,j,k), 'l' )
                        end if
                        if(isN) then
                            viscNtemp = viscosity( xNtemp, xiNbarx(1,j,k), Pw(1,j,k), 'g' )
                        end if
                    end if
                    if(Uwx(1,j,k) > 0) then
                        lambdawx(1,j,k) = Kxxbar(1,j,k)*computekr_W(Swtemp)/viscWtemp
                    elseif(Uwx(1,j,k) < 0) then
                        lambdawx(1,j,k) = Kxxbar(1,j,k)*computekr_W(Sw(1,j,k))/viscW(1,j,k)
                        xiWbarx(1,j,k) = xiW(1,j,k)
                        densiWbarx(1,j,k) = densiW(1,j,k)
                        xWbarx(1:Nc,1,j,k) = xW(1:Nc,1,j,k)
                    elseif(isDiriX(1,ylower+j-1,zlower+k-1) == 1) then
                        lambdawx(1,j,k) = Kxxbar(1,j,k)*computekr_W(Sw(1,j,k))/viscW(1,j,k)
                        xiWbarx(1,j,k) = xiW(1,j,k)
                        densiWbarx(1,j,k) = densiW(1,j,k)
                        xWbarx(1:Nc,1,j,k) = xW(1:Nc,1,j,k)
                    end if
                    if(Unx(1,j,k) > 0) then
                        lambdanx(1,j,k) = Kxxbar(1,j,k)*computekr_N(Swtemp)/viscNtemp
                    elseif(Unx(1,j,k) < 0) then
                        lambdanx(1,j,k) = Kxxbar(1,j,k)*computekr_N(Sw(1,j,k))/viscN(1,j,k)
                        xiNbarx(1,j,k) = xiN(1,j,k)
                        densiNbarx(1,j,k) = densiN(1,j,k)
                        xNbarx(1:Nc,1,j,k) = xN(1:Nc,1,j,k)
                    elseif(isDiriX(1,ylower+j-1,zlower+k-1) == 1) then
                        lambdanx(1,j,k) = Kxxbar(1,j,k)*computekr_N(Sw(1,j,k))/viscN(1,j,k)
                        xiNbarx(1,j,k) = xiN(1,j,k)
                        densiNbarx(1,j,k) = densiN(1,j,k)
                        xNbarx(1:Nc,1,j,k) = xN(1:Nc,1,j,k)
                    end if
                end do
            end do
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 2, localncols+1
                        if(Uwx(i,j,k) > 0) then
                            lambdawx(i,j,k) = Kxxbar(i,j,k)*computekr_W(Sw(i-1,j,k))/viscW(i-1,j,k)
                            xiWbarx(i,j,k) = xiW(i-1,j,k)
                            densiWbarx(i,j,k) = densiW(i-1,j,k)
                            xWbarx(1:Nc,i,j,k) = xW(1:Nc,i-1,j,k)
                        else
                            lambdawx(i,j,k) = Kxxbar(i,j,k)*computekr_W(Sw(i,j,k))/viscW(i,j,k)
                            xiWbarx(i,j,k) = xiW(i,j,k)
                            densiWbarx(i,j,k) = densiW(i,j,k)
                            xWbarx(1:Nc,i,j,k) = xW(1:Nc,i,j,k)
                        end if
                        if(Unx(i,j,k) > 0) then
                            lambdanx(i,j,k) = Kxxbar(i,j,k)*computekr_N(Sw(i-1,j,k))/viscN(i-1,j,k)
                            xiNbarx(i,j,k) = xiN(i-1,j,k)
                            densiNbarx(i,j,k) = densiN(i-1,j,k)
                            xNbarx(1:Nc,i,j,k) = xN(1:Nc,i-1,j,k)
                        else
                            lambdanx(i,j,k) = Kxxbar(i,j,k)*computekr_N(Sw(i,j,k))/viscN(i,j,k)
                            xiNbarx(i,j,k) = xiN(i,j,k)
                            densiNbarx(i,j,k) = densiN(i,j,k)
                            xNbarx(1:Nc,i,j,k) = xN(1:Nc,i,j,k)
                        end if
                    end do
                end do
            end do
        elseif(pcol==pncols) then
            do k = 1, localnlays
                do j = 1, localnrows
                    if((Uwx(localncols+1,j,k) < 0).or.(Unx(localncols+1,j,k) < 0)) then
                        ztemp(1:Nc) = zBdryX(1:Nc,2,ylower+j-1,zlower+k-1)
#ifdef SPARSE
                        call flashcalculation_fullgrid( Pw(localncols,j,k), ztemp, xWtemp, &!
                            xNtemp, xiWbarx(localncols+1,j,k), xiNbarx(localncols+1,j,k), densiWbarx(localncols+1,j,k), &!
                            densiNbarx(localncols+1,j,k), Swtemp, vtemp, ignore1, isW, isN )
#else
                        call flashcalculation(Pw(localncols,j,k), ztemp, xWtemp, &!
                            xNtemp, xiWbarx(localncols+1,j,k), xiNbarx(localncols+1,j,k), densiWbarx(localncols+1,j,k), &!
                            densiNbarx(localncols+1,j,k), Swtemp, vtemp, ignore1, isW, isN, isRea)
#endif
                        xWbarx(1:Nc,localncols+1,j,k) = xWtemp(1:Nc)
                        xNbarx(1:Nc,localncols+1,j,k) = xNtemp(1:Nc)
                        if(isW) then
                            viscWtemp = viscosity( xWtemp, xiWbarx(localncols+1,j,k), Pw(localncols,j,k), 'l' )
                        end if
                        if(isN) then
                            viscNtemp = viscosity( xNtemp, xiNbarx(localncols+1,j,k), Pw(localncols,j,k), 'g' )
                        end if
                    end if
                    if(Uwx(localncols+1,j,k) < 0) then
                        lambdawx(localncols+1,j,k) = Kxxbar(localncols+1,j,k)*computekr_W(Swtemp)/viscWtemp
                    elseif(Uwx(localncols+1,j,k) > 0) then
                        lambdawx(localncols+1,j,k) = Kxxbar(localncols+1,j,k)*computekr_W(Sw(localncols,j,k))/viscW(localncols,j,k)
                        xiWbarx(localncols+1,j,k) = xiW(localncols,j,k)
                        densiWbarx(localncols+1,j,k) = densiW(localncols,j,k)
                        xWbarx(1:Nc,localncols+1,j,k) = xW(1:Nc,localncols,j,k)
                    elseif(isDiriX(2,ylower+j-1,zlower+k-1) == 1) then
                        lambdawx(localncols+1,j,k) = Kxxbar(localncols+1,j,k)*computekr_W(Sw(localncols,j,k))/viscW(localncols,j,k)
                        xiWbarx(localncols+1,j,k) = xiW(localncols,j,k)
                        densiWbarx(localncols+1,j,k) = densiW(localncols,j,k)
                        xWbarx(1:Nc,localncols+1,j,k) = xW(1:Nc,localncols,j,k)
                    end if
                    if(Unx(localncols+1,j,k) < 0) then
                        lambdanx(localncols+1,j,k) = Kxxbar(localncols+1,j,k)*computekr_N(Swtemp)/viscNtemp
                    elseif(Unx(localncols+1,j,k) > 0) then
                        lambdanx(localncols+1,j,k) = Kxxbar(localncols+1,j,k)*computekr_N(Sw(localncols,j,k))/viscN(localncols,j,k)
                        xiNbarx(localncols+1,j,k) = xiN(localncols,j,k)
                        densiNbarx(localncols+1,j,k) = densiN(localncols,j,k)
                        xNbarx(1:Nc,localncols+1,j,k) = xN(1:Nc,localncols,j,k)
                    elseif(isDiriX(2,ylower+j-1,zlower+k-1) == 1) then
                        lambdanx(localncols+1,j,k) = Kxxbar(localncols+1,j,k)*computekr_N(Sw(localncols,j,k))/viscN(localncols,j,k)
                        xiNbarx(localncols+1,j,k) = xiN(localncols,j,k)
                        densiNbarx(localncols+1,j,k) = densiN(localncols,j,k)
                        xNbarx(1:Nc,localncols+1,j,k) = xN(1:Nc,localncols,j,k)
                    end if
                end do
            end do
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        if(Uwx(i,j,k) > 0) then
                            lambdawx(i,j,k) = Kxxbar(i,j,k)*computekr_W(Sw(i-1,j,k))/viscW(i-1,j,k)
                            xiWbarx(i,j,k) = xiW(i-1,j,k)
                            densiWbarx(i,j,k) = densiW(i-1,j,k)
                            xWbarx(1:Nc,i,j,k) = xW(1:Nc,i-1,j,k)
                        else
                            lambdawx(i,j,k) = Kxxbar(i,j,k)*computekr_W(Sw(i,j,k))/viscW(i,j,k)
                            xiWbarx(i,j,k) = xiW(i,j,k)
                            densiWbarx(i,j,k) = densiW(i,j,k)
                            xWbarx(1:Nc,i,j,k) = xW(1:Nc,i,j,k)
                        end if
                        if(Unx(i,j,k) > 0) then
                            lambdanx(i,j,k) = Kxxbar(i,j,k)*computekr_N(Sw(i-1,j,k))/viscN(i-1,j,k)
                            xiNbarx(i,j,k) = xiN(i-1,j,k)
                            densiNbarx(i,j,k) = densiN(i-1,j,k)
                            xNbarx(1:Nc,i,j,k) = xN(1:Nc,i-1,j,k)
                        else
                            lambdanx(i,j,k) = Kxxbar(i,j,k)*computekr_N(Sw(i,j,k))/viscN(i,j,k)
                            xiNbarx(i,j,k) = xiN(i,j,k)
                            densiNbarx(i,j,k) = densiN(i,j,k)
                            xNbarx(1:Nc,i,j,k) = xN(1:Nc,i,j,k)
                        end if
                    end do
                end do
            end do
        else
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols+1
                        if(Uwx(i,j,k) > 0) then
                            lambdawx(i,j,k) = Kxxbar(i,j,k)*computekr_W(Sw(i-1,j,k))/viscW(i-1,j,k)
                            xiWbarx(i,j,k) = xiW(i-1,j,k)
                            densiWbarx(i,j,k) = densiW(i-1,j,k)
                            xWbarx(1:Nc,i,j,k) = xW(1:Nc,i-1,j,k)
                        else
                            lambdawx(i,j,k) = Kxxbar(i,j,k)*computekr_W(Sw(i,j,k))/viscW(i,j,k)
                            xiWbarx(i,j,k) = xiW(i,j,k)
                            densiWbarx(i,j,k) = densiW(i,j,k)
                            xWbarx(1:Nc,i,j,k) = xW(1:Nc,i,j,k)
                        end if
                        if(Unx(i,j,k) > 0) then
                            lambdanx(i,j,k) = Kxxbar(i,j,k)*computekr_N(Sw(i-1,j,k))/viscN(i-1,j,k)
                            xiNbarx(i,j,k) = xiN(i-1,j,k)
                            densiNbarx(i,j,k) = densiN(i-1,j,k)
                            xNbarx(1:Nc,i,j,k) = xN(1:Nc,i-1,j,k)
                        else
                            lambdanx(i,j,k) = Kxxbar(i,j,k)*computekr_N(Sw(i,j,k))/viscN(i,j,k)
                            xiNbarx(i,j,k) = xiN(i,j,k)
                            densiNbarx(i,j,k) = densiN(i,j,k)
                            xNbarx(1:Nc,i,j,k) = xN(1:Nc,i,j,k)
                        end if
                    end do
                end do
            end do
        end if

        if(pnrows == 1) then
            do k = 1, localnlays
                do i = 1, localncols
                    if((Uwy(i,1,k) > 0).or.(Uny(i,1,k) > 0)) then
                        ztemp(1:Nc) = zBdryY(1:Nc,xlower+i-1,1,zlower+k-1)
#ifdef SPARSE
                        call flashcalculation_fullgrid( Pw(i,1,k), ztemp, xWtemp, xNtemp, xiWbary(i,1,k), xiNbary(i,1,k), &!
                            densiWbary(i,1,k), densiNbary(i,1,k), Swtemp, vtemp, ignore1, isW, isN )
#else
                        call flashcalculation(Pw(i,1,k), ztemp, xWtemp, xNtemp, xiWbary(i,1,k), xiNbary(i,1,k), &!
                            densiWbary(i,1,k), densiNbary(i,1,k), Swtemp, vtemp, ignore1, isW, isN, isRea)
#endif
                        xWbary(1:Nc,i,1,k) = xWtemp(1:Nc)
                        xNbary(1:Nc,i,1,k) = xNtemp(1:Nc)
                        if(isW) then
                            viscWtemp = viscosity( xWtemp, xiWbary(i,1,k), Pw(i,1,k), 'l' )
                        end if
                        if(isN) then
                            viscNtemp = viscosity( xNtemp, xiNbary(i,1,k), Pw(i,1,k), 'g' )
                        end if
                    end if
                    if(Uwy(i,1,k) > 0) then
                        lambdawy(i,1,k) = Kyybar(i,1,k)*computekr_W(Swtemp)/viscWtemp
                    elseif(Uwy(i,1,k) < 0) then
                        lambdawy(i,1,k) = Kyybar(i,1,k)*computekr_W(Sw(i,1,k))/viscW(i,1,k)
                        xiWbary(i,1,k) = xiW(i,1,k)
                        densiWbary(i,1,k) = densiW(i,1,k)
                        xWbary(1:Nc,i,1,k) = xW(1:Nc,i,1,k)
                    elseif(isDiriY(xlower+i-1,1,zlower+k-1) == 1) then
                        lambdawy(i,1,k) = Kyybar(i,1,k)*computekr_W(Sw(i,1,k))/viscW(i,1,k)
                        xiWbary(i,1,k) = xiW(i,1,k)
                        densiWbary(i,1,k) = densiW(i,1,k)
                        xWbary(1:Nc,i,1,k) = xW(1:Nc,i,1,k)
                    end if
                    if(Uny(i,1,k) > 0) then
                        lambdany(i,1,k) = Kyybar(i,1,k)*computekr_N(Swtemp)/viscNtemp
                    elseif(Uny(i,1,k) < 0) then
                        lambdany(i,1,k) = Kyybar(i,1,k)*computekr_N(Sw(i,1,k))/viscN(i,1,k)
                        xiNbary(i,1,k) = xiN(i,1,k)
                        densiNbary(i,1,k) = densiN(i,1,k)
                        xNbary(1:Nc,i,1,k) = xN(1:Nc,i,1,k)
                    elseif(isDiriY(xlower+i-1,1,zlower+k-1) == 1) then
                        lambdany(i,1,k) = Kyybar(i,1,k)*computekr_N(Sw(i,1,k))/viscN(i,1,k)
                        xiNbary(i,1,k) = xiN(i,1,k)
                        densiNbary(i,1,k) = densiN(i,1,k)
                        xNbary(1:Nc,i,1,k) = xN(1:Nc,i,1,k)
                    end if
                end do
            end do
            do k = 1, localnlays
                do i = 1, localncols
                    if((Uwy(i,localnrows+1,k) < 0).or.(Uny(i,localnrows+1,k) < 0)) then
                        ztemp(1:Nc) = zBdryY(1:Nc,xlower+i-1,2,zlower+k-1)
#ifdef SPARSE
                        call flashcalculation_fullgrid( Pw(i,localnrows,k), ztemp, xWtemp, &!
                            xNtemp, xiWbary(i,localnrows+1,k), xiNbary(i,localnrows+1,k), densiWbary(i,localnrows+1,k), &!
                            densiNbary(i,localnrows+1,k), Swtemp, vtemp, ignore1, isW, isN )
#else
                        call flashcalculation(Pw(i,localnrows,k), ztemp, xWtemp, &!
                            xNtemp, xiWbary(i,localnrows+1,k), xiNbary(i,localnrows+1,k), densiWbary(i,localnrows+1,k), &!
                            densiNbary(i,localnrows+1,k), Swtemp, vtemp, ignore1, isW, isN, isRea)
#endif
                        xWbary(1:Nc,i,localnrows+1,k) = xWtemp(1:Nc)
                        xNbary(1:Nc,i,localnrows+1,k) = xNtemp(1:Nc)
                        if(isW) then
                            viscWtemp = viscosity( xWtemp, xiWbary(i,localnrows+1,k), Pw(i,localnrows,k), 'l' )
                        end if
                        if(isN) then
                            viscNtemp = viscosity( xNtemp, xiNbary(i,localnrows+1,k), Pw(i,localnrows,k), 'g' )
                        end if
                    end if
                    if(Uwy(i,localnrows+1,k) < 0) then
                        lambdawy(i,localnrows+1,k) = Kyybar(i,localnrows+1,k)*computekr_W(Swtemp)/viscWtemp
                    elseif(Uwy(i,localnrows+1,k) > 0) then
                        lambdawy(i,localnrows+1,k) = Kyybar(i,localnrows+1,k)*computekr_W(Sw(i,localnrows,k))/viscW(i,localnrows,k)
                        xiWbary(i,localnrows+1,k) = xiW(i,localnrows,k)
                        densiWbary(i,localnrows+1,k) = densiW(i,localnrows,k)
                        xWbary(1:Nc,i,localnrows+1,k) = xW(1:Nc,i,localnrows,k)
                    elseif(isDiriY(xlower+i-1,2,zlower+k-1) == 1) then
                        lambdawy(i,localnrows+1,k) = Kyybar(i,localnrows+1,k)*computekr_W(Sw(i,localnrows,k))/viscW(i,localnrows,k)
                        xiWbary(i,localnrows+1,k) = xiW(i,localnrows,k)
                        densiWbary(i,localnrows+1,k) = densiW(i,localnrows,k)
                        xWbary(1:Nc,i,localnrows+1,k) = xW(1:Nc,i,localnrows,k)
                    end if
                    if(Uny(i,localnrows+1,k) < 0) then
                        lambdany(i,localnrows+1,k) = Kyybar(i,localnrows+1,k)*computekr_N(Swtemp)/viscNtemp
                    elseif(Uny(i,localnrows+1,k) > 0) then
                        lambdany(i,localnrows+1,k) = Kyybar(i,localnrows+1,k)*computekr_N(Sw(i,localnrows,k))/viscN(i,localnrows,k)
                        xiNbary(i,localnrows+1,k) = xiN(i,localnrows,k)
                        densiNbary(i,localnrows+1,k) = densiN(i,localnrows,k)
                        xNbary(1:Nc,i,localnrows+1,k) = xN(1:Nc,i,localnrows,k)
                    elseif(isDiriY(xlower+i-1,2,zlower+k-1) == 1) then
                        lambdany(i,localnrows+1,k) = Kyybar(i,localnrows+1,k)*computekr_N(Sw(i,localnrows,k))/viscN(i,localnrows,k)
                        xiNbary(i,localnrows+1,k) = xiN(i,localnrows,k)
                        densiNbary(i,localnrows+1,k) = densiN(i,localnrows,k)
                        xNbary(1:Nc,i,localnrows+1,k) = xN(1:Nc,i,localnrows,k)
                    end if
                end do
            end do
            do k = 1, localnlays
                do j = 2, localnrows
                    do i = 1, localncols
                        if(Uwy(i,j,k) > 0) then
                            lambdawy(i,j,k) = Kyybar(i,j,k)*computekr_W(Sw(i,j-1,k))/viscW(i,j-1,k)
                            xiWbary(i,j,k) = xiW(i,j-1,k)
                            densiWbary(i,j,k) = densiW(i,j-1,k)
                            xWbary(1:Nc,i,j,k) = xW(1:Nc,i,j-1,k)
                        else
                            lambdawy(i,j,k) = Kyybar(i,j,k)*computekr_W(Sw(i,j,k))/viscW(i,j,k)
                            xiWbary(i,j,k) = xiW(i,j,k)
                            densiWbary(i,j,k) = densiW(i,j,k)
                            xWbary(1:Nc,i,j,k) = xW(1:Nc,i,j,k)
                        end if
                        if(Uny(i,j,k) > 0) then
                            lambdany(i,j,k) = Kyybar(i,j,k)*computekr_N(Sw(i,j-1,k))/viscN(i,j-1,k)
                            xiNbary(i,j,k) = xiN(i,j-1,k)
                            densiNbary(i,j,k) = densiN(i,j-1,k)
                            xNbary(1:Nc,i,j,k) = xN(1:Nc,i,j-1,k)
                        else
                            lambdany(i,j,k) = Kyybar(i,j,k)*computekr_N(Sw(i,j,k))/viscN(i,j,k)
                            xiNbary(i,j,k) = xiN(i,j,k)
                            densiNbary(i,j,k) = densiN(i,j,k)
                            xNbary(1:Nc,i,j,k) = xN(1:Nc,i,j,k)
                        end if
                    end do
                end do
            end do
        elseif(prow==1) then
            do k = 1, localnlays
                do i = 1, localncols
                    if((Uwy(i,1,k) > 0).or.(Uny(i,1,k) > 0)) then
                        ztemp(1:Nc) = zBdryY(1:Nc,xlower+i-1,1,zlower+k-1)
#ifdef SPARSE
                        call flashcalculation_fullgrid( Pw(i,1,k), ztemp, xWtemp, xNtemp, xiWbary(i,1,k), xiNbary(i,1,k), &!
                            densiWbary(i,1,k), densiNbary(i,1,k), Swtemp, vtemp, ignore1, isW, isN )
#else
                        call flashcalculation(Pw(i,1,k), ztemp, xWtemp, xNtemp, xiWbary(i,1,k), &!
                            xiNbary(i,1,k), densiWbary(i,1,k), densiNbary(i,1,k), Swtemp, vtemp, ignore1, isW, isN, isRea)
#endif
                        xWbary(1:Nc,i,1,k) = xWtemp(1:Nc)
                        xNbary(1:Nc,i,1,k) = xNtemp(1:Nc)
                        if(isW) then
                            viscWtemp = viscosity( xWtemp, xiWbary(i,1,k), Pw(i,1,k), 'l' )
                        end if
                        if(isN) then
                            viscNtemp = viscosity( xNtemp, xiNbary(i,1,k), Pw(i,1,k), 'g' )
                        end if
                    end if
                    if(Uwy(i,1,k) > 0) then
                        lambdawy(i,1,k) = Kyybar(i,1,k)*computekr_W(Swtemp)/viscWtemp
                    elseif(Uwy(i,1,k) < 0) then
                        lambdawy(i,1,k) = Kyybar(i,1,k)*computekr_W(Sw(i,1,k))/viscW(i,1,k)
                        xiWbary(i,1,k) = xiW(i,1,k)
                        densiWbary(i,1,k) = densiW(i,1,k)
                        xWbary(1:Nc,i,1,k) = xW(1:Nc,i,1,k)
                    elseif(isDiriY(xlower+i-1,1,zlower+k-1) == 1) then
                        lambdawy(i,1,k) = Kyybar(i,1,k)*computekr_W(Sw(i,1,k))/viscW(i,1,k)
                        xiWbary(i,1,k) = xiW(i,1,k)
                        densiWbary(i,1,k) = densiW(i,1,k)
                        xWbary(1:Nc,i,1,k) = xW(1:Nc,i,1,k)
                    end if
                    if(Uny(i,1,k) > 0) then
                        lambdany(i,1,k) = Kyybar(i,1,k)*computekr_N(Swtemp)/viscNtemp
                    elseif(Uny(i,1,k) < 0) then
                        lambdany(i,1,k) = Kyybar(i,1,k)*computekr_N(Sw(i,1,k))/viscN(i,1,k)
                        xiNbary(i,1,k) = xiN(i,1,k)
                        densiNbary(i,1,k) = densiN(i,1,k)
                        xNbary(1:Nc,i,1,k) = xN(1:Nc,i,1,k)
                    elseif(isDiriY(xlower+i-1,1,zlower+k-1) == 1) then
                        lambdany(i,1,k) = Kyybar(i,1,k)*computekr_N(Sw(i,1,k))/viscN(i,1,k)
                        xiNbary(i,1,k) = xiN(i,1,k)
                        densiNbary(i,1,k) = densiN(i,1,k)
                        xNbary(1:Nc,i,1,k) = xN(1:Nc,i,1,k)
                    end if
                end do
            end do
            do k = 1, localnlays
                do j = 2, localnrows+1
                    do i = 1, localncols
                        if(Uwy(i,j,k) > 0) then
                            lambdawy(i,j,k) = Kyybar(i,j,k)*computekr_W(Sw(i,j-1,k))/viscW(i,j-1,k)
                            xiWbary(i,j,k) = xiW(i,j-1,k)
                            densiWbary(i,j,k) = densiW(i,j-1,k)
                            xWbary(1:Nc,i,j,k) = xW(1:Nc,i,j-1,k)
                        else
                            lambdawy(i,j,k) = Kyybar(i,j,k)*computekr_W(Sw(i,j,k))/viscW(i,j,k)
                            xiWbary(i,j,k) = xiW(i,j,k)
                            densiWbary(i,j,k) = densiW(i,j,k)
                            xWbary(1:Nc,i,j,k) = xW(1:Nc,i,j,k)
                        end if
                        if(Uny(i,j,k) > 0) then
                            lambdany(i,j,k) = Kyybar(i,j,k)*computekr_N(Sw(i,j-1,k))/viscN(i,j-1,k)
                            xiNbary(i,j,k) = xiN(i,j-1,k)
                            densiNbary(i,j,k) = densiN(i,j-1,k)
                            xNbary(1:Nc,i,j,k) = xN(1:Nc,i,j-1,k)
                        else
                            lambdany(i,j,k) = Kyybar(i,j,k)*computekr_N(Sw(i,j,k))/viscN(i,j,k)
                            xiNbary(i,j,k) = xiN(i,j,k)
                            densiNbary(i,j,k) = densiN(i,j,k)
                            xNbary(1:Nc,i,j,k) = xN(1:Nc,i,j,k)
                        end if
                    end do
                end do
            end do
        elseif(prow==pnrows) then
            do k = 1, localnlays
                do i = 1, localncols
                    if((Uwy(i,localnrows+1,k) < 0).or.(Uny(i,localnrows+1,k) < 0)) then
                        ztemp(1:Nc) = zBdryY(1:Nc,xlower+i-1,2,zlower+k-1)
#ifdef SPARSE
                        call flashcalculation_fullgrid( Pw(i,localnrows,k), ztemp, xWtemp, &!
                            xNtemp, xiWbary(i,localnrows+1,k), xiNbary(i,localnrows+1,k), densiWbary(i,localnrows+1,k), &!
                            densiNbary(i,localnrows+1,k), Swtemp, vtemp, ignore1, isW, isN )
#else
                        call flashcalculation(Pw(i,localnrows,k), ztemp, xWtemp, &!
                            xNtemp, xiWbary(i,localnrows+1,k), xiNbary(i,localnrows+1,k), densiWbary(i,localnrows+1,k), &!
                            densiNbary(i,localnrows+1,k), Swtemp, vtemp, ignore1, isW, isN, isRea)
#endif
                        xWbary(1:Nc,i,localnrows+1,k) = xWtemp(1:Nc)
                        xNbary(1:Nc,i,localnrows+1,k) = xNtemp(1:Nc)
                        if(isW) then
                            viscWtemp = viscosity( xWtemp, xiWbary(i,localnrows+1,k), Pw(i,localnrows,k), 'l' )
                        end if
                        if(isN) then
                            viscNtemp = viscosity( xNtemp, xiNbary(i,localnrows+1,k), Pw(i,localnrows,k), 'g' )
                        end if
                    end if
                    if(Uwy(i,localnrows+1,k) < 0) then
                        lambdawy(i,localnrows+1,k) = Kyybar(i,localnrows+1,k)*computekr_W(Swtemp)/viscWtemp
                    elseif(Uwy(i,localnrows+1,k) > 0) then
                        lambdawy(i,localnrows+1,k) = Kyybar(i,localnrows+1,k)*computekr_W(Sw(i,localnrows,k))/viscW(i,localnrows,k)
                        xiWbary(i,localnrows+1,k) = xiW(i,localnrows,k)
                        densiWbary(i,localnrows+1,k) = densiW(i,localnrows,k)
                        xWbary(1:Nc,i,localnrows+1,k) = xW(1:Nc,i,localnrows,k)
                    elseif(isDiriY(xlower+i-1,2,zlower+k-1) == 1) then
                        lambdawy(i,localnrows+1,k) = Kyybar(i,localnrows+1,k)*computekr_W(Sw(i,localnrows,k))/viscW(i,localnrows,k)
                        xiWbary(i,localnrows+1,k) = xiW(i,localnrows,k)
                        densiWbary(i,localnrows+1,k) = densiW(i,localnrows,k)
                        xWbary(1:Nc,i,localnrows+1,k) = xW(1:Nc,i,localnrows,k)
                    end if
                    if(Uny(i,localnrows+1,k) < 0) then
                        lambdany(i,localnrows+1,k) = Kyybar(i,localnrows+1,k)*computekr_N(Swtemp)/viscNtemp
                    elseif(Uny(i,localnrows+1,k) > 0) then
                        lambdany(i,localnrows+1,k) = Kyybar(i,localnrows+1,k)*computekr_N(Sw(i,localnrows,k))/viscN(i,localnrows,k)
                        xiNbary(i,localnrows+1,k) = xiN(i,localnrows,k)
                        densiNbary(i,localnrows+1,k) = densiN(i,localnrows,k)
                        xNbary(1:Nc,i,localnrows+1,k) = xN(1:Nc,i,localnrows,k)
                    elseif(isDiriY(xlower+i-1,2,zlower+k-1) == 1) then
                        lambdany(i,localnrows+1,k) = Kyybar(i,localnrows+1,k)*computekr_N(Sw(i,localnrows,k))/viscN(i,localnrows,k)
                        xiNbary(i,localnrows+1,k) = xiN(i,localnrows,k)
                        densiNbary(i,localnrows+1,k) = densiN(i,localnrows,k)
                        xNbary(1:Nc,i,localnrows+1,k) = xN(1:Nc,i,localnrows,k)
                    end if
                end do
            end do
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        if(Uwy(i,j,k) > 0) then
                            lambdawy(i,j,k) = Kyybar(i,j,k)*computekr_W(Sw(i,j-1,k))/viscW(i,j-1,k)
                            xiWbary(i,j,k) = xiW(i,j-1,k)
                            densiWbary(i,j,k) = densiW(i,j-1,k)
                            xWbary(1:Nc,i,j,k) = xW(1:Nc,i,j-1,k)
                        else
                            lambdawy(i,j,k) = Kyybar(i,j,k)*computekr_W(Sw(i,j,k))/viscW(i,j,k)
                            xiWbary(i,j,k) = xiW(i,j,k)
                            densiWbary(i,j,k) = densiW(i,j,k)
                            xWbary(1:Nc,i,j,k) = xW(1:Nc,i,j,k)
                        end if
                        if(Uny(i,j,k) > 0) then
                            lambdany(i,j,k) = Kyybar(i,j,k)*computekr_N(Sw(i,j-1,k))/viscN(i,j-1,k)
                            xiNbary(i,j,k) = xiN(i,j-1,k)
                            densiNbary(i,j,k) = densiN(i,j-1,k)
                            xNbary(1:Nc,i,j,k) = xN(1:Nc,i,j-1,k)
                        else
                            lambdany(i,j,k) = Kyybar(i,j,k)*computekr_N(Sw(i,j,k))/viscN(i,j,k)
                            xiNbary(i,j,k) = xiN(i,j,k)
                            densiNbary(i,j,k) = densiN(i,j,k)
                            xNbary(1:Nc,i,j,k) = xN(1:Nc,i,j,k)
                        end if
                    end do
                end do
            end do
        else
            do k = 1, localnlays
                do j = 1, localnrows+1
                    do i = 1, localncols
                        if(Uwy(i,j,k) > 0) then
                            lambdawy(i,j,k) = Kyybar(i,j,k)*computekr_W(Sw(i,j-1,k))/viscW(i,j-1,k)
                            xiWbary(i,j,k) = xiW(i,j-1,k)
                            densiWbary(i,j,k) = densiW(i,j-1,k)
                            xWbary(1:Nc,i,j,k) = xW(1:Nc,i,j-1,k)
                        else
                            lambdawy(i,j,k) = Kyybar(i,j,k)*computekr_W(Sw(i,j,k))/viscW(i,j,k)
                            xiWbary(i,j,k) = xiW(i,j,k)
                            densiWbary(i,j,k) = densiW(i,j,k)
                            xWbary(1:Nc,i,j,k) = xW(1:Nc,i,j,k)
                        end if
                        if(Uny(i,j,k) > 0) then
                            lambdany(i,j,k) = Kyybar(i,j,k)*computekr_N(Sw(i,j-1,k))/viscN(i,j-1,k)
                            xiNbary(i,j,k) = xiN(i,j-1,k)
                            densiNbary(i,j,k) = densiN(i,j-1,k)
                            xNbary(1:Nc,i,j,k) = xN(1:Nc,i,j-1,k)
                        else
                            lambdany(i,j,k) = Kyybar(i,j,k)*computekr_N(Sw(i,j,k))/viscN(i,j,k)
                            xiNbary(i,j,k) = xiN(i,j,k)
                            densiNbary(i,j,k) = densiN(i,j,k)
                            xNbary(1:Nc,i,j,k) = xN(1:Nc,i,j,k)
                        end if
                    end do
                end do
            end do
        end if

        if(pnlays == 1) then
            do j = 1, localnrows
                do i = 1, localncols
                    if((Uwz(i,j,1) > 0).or.(Unz(i,j,1) > 0)) then
                        ztemp(1:Nc) = zBdryZ(1:Nc,xlower+i-1,ylower+j-1,1)
#ifdef SPARSE
                        call flashcalculation_fullgrid( Pw(i,j,1), ztemp, xWtemp, &!
                            xNtemp, xiWbarz(i,j,1), xiNbarz(i,j,1), densiWbarz(i,j,1), densiNbarz(i,j,1), &!
                            Swtemp, vtemp, ignore1, isW, isN )
#else
                        call flashcalculation(Pw(i,j,1), ztemp, xWtemp, &!
                            xNtemp, xiWbarz(i,j,1), xiNbarz(i,j,1), densiWbarz(i,j,1), densiNbarz(i,j,1), &!
                            Swtemp, vtemp, ignore1, isW, isN, isRea)
#endif
                        xWbarz(1:Nc,i,j,1) = xWtemp(1:Nc)
                        xNbarz(1:Nc,i,j,1) = xNtemp(1:Nc)
                        if(isW) then
                            viscWtemp = viscosity( xWtemp, xiWbarz(i,j,1), Pw(i,j,1), 'l' )
                        end if
                        if(isN) then
                            viscNtemp = viscosity( xNtemp, xiNbarz(i,j,1), Pw(i,j,1), 'g' )
                        end if
                    end if
                    if(Uwz(i,j,1) > 0) then
                        lambdawz(i,j,1) = Kzzbar(i,j,1)*computekr_W(Swtemp)/viscWtemp
                    elseif(Uwz(i,j,1) < 0) then
                        lambdawz(i,j,1) = Kzzbar(i,j,1)*computekr_W(Sw(i,j,1))/viscW(i,j,1)
                        xiWbarz(i,j,1) = xiW(i,j,1)
                        densiWbarz(i,j,1) = densiW(i,j,1)
                        xWbarz(1:Nc,i,j,1) = xW(1:Nc,i,j,1)
                    elseif(isDiriZ(xlower+i-1,ylower+j-1,1) == 1) then
                        lambdawz(i,j,1) = Kzzbar(i,j,1)*computekr_W(Sw(i,j,1))/viscW(i,j,1)
                        xiWbarz(i,j,1) = xiW(i,j,1)
                        densiWbarz(i,j,1) = densiW(i,j,1)
                        xWbarz(1:Nc,i,j,1) = xW(1:Nc,i,j,1)
                    end if
                    if(Unz(i,j,1) > 0) then
                        lambdanz(i,j,1) = Kzzbar(i,j,1)*computekr_N(Swtemp)/viscNtemp
                    elseif(Unz(i,j,1) < 0) then
                        lambdanz(i,j,1) = Kzzbar(i,j,1)*computekr_N(Sw(i,j,1))/viscN(i,j,1)
                        xiNbarz(i,j,1) = xiN(i,j,1)
                        densiNbarz(i,j,1) = densiN(i,j,1)
                        xNbarz(1:Nc,i,j,1) = xN(1:Nc,i,j,1)
                    elseif(isDiriZ(xlower+i-1,ylower+j-1,1) == 1) then
                        lambdanz(i,j,1) = Kzzbar(i,j,1)*computekr_N(Sw(i,j,1))/viscN(i,j,1)
                        xiNbarz(i,j,1) = xiN(i,j,1)
                        densiNbarz(i,j,1) = densiN(i,j,1)
                        xNbarz(1:Nc,i,j,1) = xN(1:Nc,i,j,1)
                    end if
                end do
            end do
            do j = 1, localnrows
                do i = 1, localncols
                    if((Uwz(i,j,localnlays+1) < 0).or.(Unz(i,j,localnlays+1) < 0)) then
                        ztemp(1:Nc) = zBdryZ(1:Nc,xlower+i-1,ylower+j-1,2)
#ifdef SPARSE
                        call flashcalculation_fullgrid( Pw(i,j,localnlays), ztemp, xWtemp, &!
                            xNtemp, xiWbarz(i,j,localnlays+1), xiNbarz(i,j,localnlays+1), densiWbarz(i,j,localnlays+1), &!
                            densiNbarz(i,j,localnlays+1), Swtemp, vtemp, ignore1, isW, isN )
#else
                        call flashcalculation(Pw(i,j,localnlays), ztemp, xWtemp, &!
                            xNtemp, xiWbarz(i,j,localnlays+1), xiNbarz(i,j,localnlays+1), densiWbarz(i,j,localnlays+1), &!
                            densiNbarz(i,j,localnlays+1), Swtemp, vtemp, ignore1, isW, isN, isRea)
#endif
                        xWbarz(1:Nc,i,j,localnlays+1) = xWtemp(1:Nc)
                        xNbarz(1:Nc,i,j,localnlays+1) = xNtemp(1:Nc)
                        if(isW) then
                            viscWtemp = viscosity( xWtemp, xiWbarz(i,j,localnlays+1), Pw(i,j,localnlays), 'l' )
                        end if
                        if(isN) then
                            viscNtemp = viscosity( xNtemp, xiNbarz(i,j,localnlays+1), Pw(i,j,localnlays), 'g' )
                        end if
                    end if
                    if(Uwz(i,j,localnlays+1) < 0) then
                        lambdawz(i,j,localnlays+1) = Kzzbar(i,j,localnlays+1)*computekr_W(Swtemp)/viscWtemp
                    elseif(Uwz(i,j,localnlays+1) > 0) then
                        lambdawz(i,j,localnlays+1) = Kzzbar(i,j,localnlays+1)*computekr_W(Sw(i,j,localnlays))/viscW(i,j,localnlays)
                        xiWbarz(i,j,localnlays+1) = xiW(i,j,localnlays)
                        densiWbarz(i,j,localnlays+1) = densiW(i,j,localnlays)
                        xWbarz(1:Nc,i,j,localnlays+1) = xW(1:Nc,i,j,localnlays)
                    elseif(isDiriZ(xlower+i-1,ylower+j-1,2) == 1) then
                        lambdawz(i,j,localnlays+1) = Kzzbar(i,j,localnlays+1)*computekr_W(Sw(i,j,localnlays))/viscW(i,j,localnlays)
                        xiWbarz(i,j,localnlays+1) = xiW(i,j,localnlays)
                        densiWbarz(i,j,localnlays+1) = densiW(i,j,localnlays)
                        xWbarz(1:Nc,i,j,localnlays+1) = xW(1:Nc,i,j,localnlays)
                    end if
                    if(Unz(i,j,localnlays+1) < 0) then
                        lambdanz(i,j,localnlays+1) = Kzzbar(i,j,localnlays+1)*computekr_N(Swtemp)/viscNtemp
                    elseif(Unz(i,j,localnlays+1) > 0) then
                        lambdanz(i,j,localnlays+1) = Kzzbar(i,j,localnlays+1)*computekr_N(Sw(i,j,localnlays))/viscN(i,j,localnlays)
                        xiNbarz(i,j,localnlays+1) = xiN(i,j,localnlays)
                        densiNbarz(i,j,localnlays+1) = densiN(i,j,localnlays)
                        xNbarz(1:Nc,i,j,localnlays+1) = xN(1:Nc,i,j,localnlays)
                    elseif(isDiriZ(xlower+i-1,ylower+j-1,2) == 1) then
                        lambdanz(i,j,localnlays+1) = Kzzbar(i,j,localnlays+1)*computekr_N(Sw(i,j,localnlays))/viscN(i,j,localnlays)
                        xiNbarz(i,j,localnlays+1) = xiN(i,j,localnlays)
                        densiNbarz(i,j,localnlays+1) = densiN(i,j,localnlays)
                        xNbarz(1:Nc,i,j,localnlays+1) = xN(1:Nc,i,j,localnlays)
                    end if
                end do
            end do
            do k = 2, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        if(Uwz(i,j,k) > 0) then
                            lambdawz(i,j,k) = Kzzbar(i,j,k)*computekr_W(Sw(i,j,k-1))/viscW(i,j,k-1)
                            xiWbarz(i,j,k) = xiW(i,j,k-1)
                            densiWbarz(i,j,k) = densiW(i,j,k-1)
                            xWbarz(1:Nc,i,j,k) = xW(1:Nc,i,j,k-1)
                        else
                            lambdawz(i,j,k) = Kzzbar(i,j,k)*computekr_W(Sw(i,j,k))/viscW(i,j,k)
                            xiWbarz(i,j,k) = xiW(i,j,k)
                            densiWbarz(i,j,k) = densiW(i,j,k)
                            xWbarz(1:Nc,i,j,k) = xW(1:Nc,i,j,k)
                        end if
                        if(Unz(i,j,k) > 0) then
                            lambdanz(i,j,k) = Kzzbar(i,j,k)*computekr_N(Sw(i,j,k-1))/viscN(i,j,k-1)
                            xiNbarz(i,j,k) = xiN(i,j,k-1)
                            densiNbarz(i,j,k) = densiN(i,j,k-1)
                            xNbarz(1:Nc,i,j,k) = xN(1:Nc,i,j,k-1)
                        else
                            lambdanz(i,j,k) = Kzzbar(i,j,k)*computekr_N(Sw(i,j,k))/viscN(i,j,k)
                            xiNbarz(i,j,k) = xiN(i,j,k)
                            densiNbarz(i,j,k) = densiN(i,j,k)
                            xNbarz(1:Nc,i,j,k) = xN(1:Nc,i,j,k)
                        end if
                    end do
                end do
            end do
        elseif(play == 1) then
            do j = 1, localnrows
                do i = 1, localncols
                    if((Uwz(i,j,1) > 0).or.(Unz(i,j,1) > 0)) then
                        ztemp(1:Nc) = zBdryZ(1:Nc,xlower+i-1,ylower+j-1,1)
#ifdef SPARSE
                        call flashcalculation_fullgrid( Pw(i,j,1), ztemp, xWtemp, xNtemp, xiWbarz(i,j,1), xiNbarz(i,j,1), &!
                            densiWbarz(i,j,1), densiNbarz(i,j,1), Swtemp, vtemp, ignore1, isW, isN )
#else
                        call flashcalculation(Pw(i,j,1), ztemp, xWtemp, xNtemp, xiWbarz(i,j,1), &!
                            xiNbarz(i,j,1), densiWbarz(i,j,1), densiNbarz(i,j,1), Swtemp, vtemp, ignore1, isW, isN, isRea)
#endif
                        xWbarz(1:Nc,i,j,1) = xWtemp(1:Nc)
                        xNbarz(1:Nc,i,j,1) = xNtemp(1:Nc)
                        if(isW) then
                            viscWtemp = viscosity( xWtemp, xiWbarz(i,j,1), Pw(i,j,1), 'l' )
                        end if
                        if(isN) then
                            viscNtemp = viscosity( xNtemp, xiNbarz(i,j,1), Pw(i,j,1), 'g' )
                        end if
                    end if
                    if(Uwz(i,j,1) > 0) then
                        lambdawz(i,j,1) = Kzzbar(i,j,1)*computekr_W(Swtemp)/viscWtemp
                    elseif(Uwz(i,j,1) < 0) then
                        lambdawz(i,j,1) = Kzzbar(i,j,1)*computekr_W(Sw(i,j,1))/viscW(i,j,1)
                        xiWbarz(i,j,1) = xiW(i,j,1)
                        densiWbarz(i,j,1) = densiW(i,j,1)
                        xWbarz(1:Nc,i,j,1) = xW(1:Nc,i,j,1)
                    elseif(isDiriZ(xlower+i-1,ylower+j-1,1) == 1) then
                        lambdawz(i,j,1) = Kzzbar(i,j,1)*computekr_W(Sw(i,j,1))/viscW(i,j,1)
                        xiWbarz(i,j,1) = xiW(i,j,1)
                        densiWbarz(i,j,1) = densiW(i,j,1)
                        xWbarz(1:Nc,i,j,1) = xW(1:Nc,i,j,1)
                    end if
                    if(Unz(i,j,1) > 0) then
                        lambdanz(i,j,1) = Kzzbar(i,j,1)*computekr_N(Swtemp)/viscNtemp
                    elseif(Unz(i,j,1) < 0) then
                        lambdanz(i,j,1) = Kzzbar(i,j,1)*computekr_N(Sw(i,j,1))/viscN(i,j,1)
                        xiNbarz(i,j,1) = xiN(i,j,1)
                        densiNbarz(i,j,1) = densiN(i,j,1)
                        xNbarz(1:Nc,i,j,1) = xN(1:Nc,i,j,1)
                    elseif(isDiriZ(xlower+i-1,ylower+j-1,1) == 1) then
                        lambdanz(i,j,1) = Kzzbar(i,j,1)*computekr_N(Sw(i,j,1))/viscN(i,j,1)
                        xiNbarz(i,j,1) = xiN(i,j,1)
                        densiNbarz(i,j,1) = densiN(i,j,1)
                        xNbarz(1:Nc,i,j,1) = xN(1:Nc,i,j,1)
                    end if
                end do
            end do
            do k = 2, localnlays+1
                do j = 1, localnrows
                    do i = 1, localncols
                        if(Uwz(i,j,k) > 0) then
                            lambdawz(i,j,k) = Kzzbar(i,j,k)*computekr_W(Sw(i,j,k-1))/viscW(i,j,k-1)
                            xiWbarz(i,j,k) = xiW(i,j,k-1)
                            densiWbarz(i,j,k) = densiW(i,j,k-1)
                            xWbarz(1:Nc,i,j,k) = xW(1:Nc,i,j,k-1)
                        else
                            lambdawz(i,j,k) = Kzzbar(i,j,k)*computekr_W(Sw(i,j,k))/viscW(i,j,k)
                            xiWbarz(i,j,k) = xiW(i,j,k)
                            densiWbarz(i,j,k) = densiW(i,j,k)
                            xWbarz(1:Nc,i,j,k) = xW(1:Nc,i,j,k)
                        end if
                        if(Unz(i,j,k) > 0) then
                            lambdanz(i,j,k) = Kzzbar(i,j,k)*computekr_N(Sw(i,j,k-1))/viscN(i,j,k-1)
                            xiNbarz(i,j,k) = xiN(i,j,k-1)
                            densiNbarz(i,j,k) = densiN(i,j,k-1)
                            xNbarz(1:Nc,i,j,k) = xN(1:Nc,i,j,k-1)
                        else
                            lambdanz(i,j,k) = Kzzbar(i,j,k)*computekr_N(Sw(i,j,k))/viscN(i,j,k)
                            xiNbarz(i,j,k) = xiN(i,j,k)
                            densiNbarz(i,j,k) = densiN(i,j,k)
                            xNbarz(1:Nc,i,j,k) = xN(1:Nc,i,j,k)
                        end if
                    end do
                end do
            end do
        elseif(play == pnlays) then
            do j = 1, localnrows
                do i = 1, localncols
                    if((Uwz(i,j,localnlays+1) < 0).or.(Unz(i,j,localnlays+1) < 0)) then
                        ztemp(1:Nc) = zBdryZ(1:Nc,xlower+i-1,ylower+j-1,2)
#ifdef SPARSE
                        call flashcalculation_fullgrid( Pw(i,j,localnlays), ztemp, xWtemp, &!
                            xNtemp, xiWbarz(i,j,localnlays+1), xiNbarz(i,j,localnlays+1), densiWbarz(i,j,localnlays+1), &!
                            densiNbarz(i,j,localnlays+1), Swtemp, vtemp, ignore1, isW, isN )
#else
                        call flashcalculation(Pw(i,j,localnlays), ztemp, xWtemp, &!
                            xNtemp, xiWbarz(i,j,localnlays+1), xiNbarz(i,j,localnlays+1), densiWbarz(i,j,localnlays+1), &!
                            densiNbarz(i,j,localnlays+1), Swtemp, vtemp, ignore1, isW, isN, isRea)
#endif
                        xWbarz(1:Nc,i,j,localnlays+1) = xWtemp(1:Nc)
                        xNbarz(1:Nc,i,j,localnlays+1) = xNtemp(1:Nc)
                        if(isW) then
                            viscWtemp = viscosity( xWtemp, xiWbarz(i,j,localnlays+1), Pw(i,j,localnlays), 'l' )
                        end if
                        if(isN) then
                            viscNtemp = viscosity( xNtemp, xiNbarz(i,j,localnlays+1), Pw(i,j,localnlays), 'g' )
                        end if
                    end if
                    if(Uwz(i,j,localnlays+1) < 0) then
                        lambdawz(i,j,localnlays+1) = Kzzbar(i,j,localnlays+1)*computekr_W(Swtemp)/viscWtemp
                    elseif(Uwz(i,j,localnlays+1) > 0) then
                        lambdawz(i,j,localnlays+1) = Kzzbar(i,j,localnlays+1)*computekr_W(Sw(i,j,localnlays))/viscW(i,j,localnlays)
                        xiWbarz(i,j,localnlays+1) = xiW(i,j,localnlays)
                        densiWbarz(i,j,localnlays+1) = densiW(i,j,localnlays)
                        xWbarz(1:Nc,i,j,localnlays+1) = xW(1:Nc,i,j,localnlays)
                    elseif(isDiriZ(xlower+i-1,ylower+j-1,2) == 1) then
                        lambdawz(i,j,localnlays+1) = Kzzbar(i,j,localnlays+1)*computekr_W(Sw(i,j,localnlays))/viscW(i,j,localnlays)
                        xiWbarz(i,j,localnlays+1) = xiW(i,j,localnlays)
                        densiWbarz(i,j,localnlays+1) = densiW(i,j,localnlays)
                        xWbarz(1:Nc,i,j,localnlays+1) = xW(1:Nc,i,j,localnlays)
                    end if
                    if(Unz(i,j,localnlays+1) < 0) then
                        lambdanz(i,j,localnlays+1) = Kzzbar(i,j,localnlays+1)*computekr_N(Swtemp)/viscNtemp
                    elseif(Unz(i,j,localnlays+1) > 0) then
                        lambdanz(i,j,localnlays+1) = Kzzbar(i,j,localnlays+1)*computekr_N(Sw(i,j,localnlays))/viscN(i,j,localnlays)
                        xiNbarz(i,j,localnlays+1) = xiN(i,j,localnlays)
                        densiNbarz(i,j,localnlays+1) = densiN(i,j,localnlays)
                        xNbarz(1:Nc,i,j,localnlays+1) = xN(1:Nc,i,j,localnlays)
                    elseif(isDiriZ(xlower+i-1,ylower+j-1,2) == 1) then
                        lambdanz(i,j,localnlays+1) = Kzzbar(i,j,localnlays+1)*computekr_N(Sw(i,j,localnlays))/viscN(i,j,localnlays)
                        xiNbarz(i,j,localnlays+1) = xiN(i,j,localnlays)
                        densiNbarz(i,j,localnlays+1) = densiN(i,j,localnlays)
                        xNbarz(1:Nc,i,j,localnlays+1) = xN(1:Nc,i,j,localnlays)
                    end if
                end do
            end do
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        if(Uwz(i,j,k) > 0) then
                            lambdawz(i,j,k) = Kzzbar(i,j,k)*computekr_W(Sw(i,j,k-1))/viscW(i,j,k-1)
                            xiWbarz(i,j,k) = xiW(i,j,k-1)
                            densiWbarz(i,j,k) = densiW(i,j,k-1)
                            xWbarz(1:Nc,i,j,k) = xW(1:Nc,i,j,k-1)
                        else
                            lambdawz(i,j,k) = Kzzbar(i,j,k)*computekr_W(Sw(i,j,k))/viscW(i,j,k)
                            xiWbarz(i,j,k) = xiW(i,j,k)
                            densiWbarz(i,j,k) = densiW(i,j,k)
                            xWbarz(1:Nc,i,j,k) = xW(1:Nc,i,j,k)
                        end if
                        if(Unz(i,j,k) > 0) then
                            lambdanz(i,j,k) = Kzzbar(i,j,k)*computekr_N(Sw(i,j,k-1))/viscN(i,j,k-1)
                            xiNbarz(i,j,k) = xiN(i,j,k-1)
                            densiNbarz(i,j,k) = densiN(i,j,k-1)
                            xNbarz(1:Nc,i,j,k) = xN(1:Nc,i,j,k-1)
                        else
                            lambdanz(i,j,k) = Kzzbar(i,j,k)*computekr_N(Sw(i,j,k))/viscN(i,j,k)
                            xiNbarz(i,j,k) = xiN(i,j,k)
                            densiNbarz(i,j,k) = densiN(i,j,k)
                            xNbarz(1:Nc,i,j,k) = xN(1:Nc,i,j,k)
                        end if
                    end do
                end do
            end do
        else
            do k = 1, localnlays+1
                do j = 1, localnrows
                    do i = 1, localncols
                        if(Uwz(i,j,k) > 0) then
                            lambdawz(i,j,k) = Kzzbar(i,j,k)*computekr_W(Sw(i,j,k-1))/viscW(i,j,k-1)
                            xiWbarz(i,j,k) = xiW(i,j,k-1)
                            densiWbarz(i,j,k) = densiW(i,j,k-1)
                            xWbarz(1:Nc,i,j,k) = xW(1:Nc,i,j,k-1)
                        else
                            lambdawz(i,j,k) = Kzzbar(i,j,k)*computekr_W(Sw(i,j,k))/viscW(i,j,k)
                            xiWbarz(i,j,k) = xiW(i,j,k)
                            densiWbarz(i,j,k) = densiW(i,j,k)
                            xWbarz(1:Nc,i,j,k) = xW(1:Nc,i,j,k)
                        end if
                        if(Unz(i,j,k) > 0) then
                            lambdanz(i,j,k) = Kzzbar(i,j,k)*computekr_N(Sw(i,j,k-1))/viscN(i,j,k-1)
                            xiNbarz(i,j,k) = xiN(i,j,k-1)
                            densiNbarz(i,j,k) = densiN(i,j,k-1)
                            xNbarz(1:Nc,i,j,k) = xN(1:Nc,i,j,k-1)
                        else
                            lambdanz(i,j,k) = Kzzbar(i,j,k)*computekr_N(Sw(i,j,k))/viscN(i,j,k)
                            xiNbarz(i,j,k) = xiN(i,j,k)
                            densiNbarz(i,j,k) = densiN(i,j,k)
                            xNbarz(1:Nc,i,j,k) = xN(1:Nc,i,j,k)
                        end if
                    end do
                end do
            end do
        end if

        deallocate(ztemp)
        deallocate(xWtemp)
        deallocate(xNtemp)
        deallocate(vtemp)

    end subroutine computeParameters

    subroutine computePres()

        implicit none
       
        real(kind=8) :: xedge, yedge, zedge, ledge, redge, dedge, uedge, fedge, bedge
        real(kind=8) :: cotwx1, cotnx1, cotwx2, cotnx2, cotwy1, cotny1, cotwy2, cotny2, cotwz1, cotnz1, cotwz2, cotnz2
        real(kind=8) :: left, right, down, up, front, back
        real(kind=8) :: sum1, sum2
        real(kind=8), dimension(:), pointer :: Pwsent
        real(kind=8), dimension(:), pointer :: sentbuffer, recvbuffer
        integer :: i, j, k, m, n, r

        integer :: ierr, errorcode, num_iter
        integer :: status(MPI_STATUS_SIZE)
        integer :: requestl, requestr, requestd, requestu, requestf, requestb
        integer :: buffer_size = MAX_BUF
        real(kind=8) :: buffer(MAX_BUF)
        real(kind=8) :: solvertimestart, solvertimefinish
        real(kind=8) :: commtimestart, commtimefinish
        real(kind=8), dimension(:), pointer :: values
        real(kind=8), dimension(:), pointer :: rhs_values, x_values

        allocate(rhs_values(local_size))
        allocate(x_values(local_size))
        allocate(values(local_size*7))
        values = 0

        n = 1
        r = 1
        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols

                    xedge = xs(xlower+i) - xs(xlower+i-1)
                    yedge = ys(ylower+j) - ys(ylower+j-1)
                    zedge = zs(zlower+k) - zs(zlower+k-1)

                    if((pcol==1).and.(i==1)) then
                        ledge = 0
                    else
                        ledge = xs(xlower+i-1) - xs(xlower+i-2)
                    end if
                    if((pcol==pncols).and.(i==localncols)) then
                        redge = 0
                    else
                        redge = xs(xlower+i+1) - xs(xlower+i)
                    end if

                    if((prow==1).and.(j==1)) then
                        dedge = 0
                    else
                        dedge = ys(ylower+j-1) - ys(ylower+j-2)
                    end if
                    if((prow==pnrows).and.(j==localnrows)) then
                        uedge = 0
                    else
                        uedge = ys(ylower+j+1) - ys(ylower+j)
                    end if

                    if((play==1).and.(k==1)) then
                        fedge = 0
                    else
                        fedge = zs(zlower+k-1) - zs(zlower+k-2)
                    end if
                    if((play==pnlays).and.(k==localnlays)) then
                        bedge = 0
                    else
                        bedge = zs(zlower+k+1) - zs(zlower+k)
                    end if

                    values(r) = poro(xlower+i-1,ylower+j-1,zlower+k-1)*Cf(i,j,k)/(timeEnd/nt)

                    rhs_values(n) = poro(xlower+i-1,ylower+j-1,zlower+k-1)*Cf(i,j,k)/(timeEnd/nt)*Pw(i,j,k)
                    do m = 1, Nc
                        rhs_values(n) = rhs_values(n) + v(m,i,j,k)*src(m,xlower+i-1,ylower+j-1,zlower+k-1)
                    end do

                    cotwx1 = 0
                    cotnx1 = 0
                    cotwx2 = 0
                    cotnx2 = 0
                    cotwy1 = 0
                    cotny1 = 0
                    cotwy2 = 0
                    cotny2 = 0
                    cotwz1 = 0
                    cotnz1 = 0
                    cotwz2 = 0
                    cotnz2 = 0
                    do m = 1, Nc
                        cotwx1 = cotwx1 + lambdawx(i+1,j,k)*xiWbarx(i+1,j,k)*xWbarx(m,i+1,j,k)*v(m,i,j,k)
                        cotnx1 = cotnx1 + lambdanx(i+1,j,k)*xiNbarx(i+1,j,k)*xNbarx(m,i+1,j,k)*v(m,i,j,k)
                        cotwx2 = cotwx2 + lambdawx(i,j,k)*xiWbarx(i,j,k)*xWbarx(m,i,j,k)*v(m,i,j,k)
                        cotnx2 = cotnx2 + lambdanx(i,j,k)*xiNbarx(i,j,k)*xNbarx(m,i,j,k)*v(m,i,j,k)
                        
                        cotwy1 = cotwy1 + lambdawy(i,j+1,k)*xiWbary(i,j+1,k)*xWbary(m,i,j+1,k)*v(m,i,j,k)
                        cotny1 = cotny1 + lambdany(i,j+1,k)*xiNbary(i,j+1,k)*xNbary(m,i,j+1,k)*v(m,i,j,k)
                        cotwy2 = cotwy2 + lambdawy(i,j,k)*xiWbary(i,j,k)*xWbary(m,i,j,k)*v(m,i,j,k)
                        cotny2 = cotny2 + lambdany(i,j,k)*xiNbary(i,j,k)*xNbary(m,i,j,k)*v(m,i,j,k)

                        cotwz1 = cotwz1 + lambdawz(i,j,k+1)*xiWbarz(i,j,k+1)*xWbarz(m,i,j,k+1)*v(m,i,j,k)
                        cotnz1 = cotnz1 + lambdanz(i,j,k+1)*xiNbarz(i,j,k+1)*xNbarz(m,i,j,k+1)*v(m,i,j,k)
                        cotwz2 = cotwz2 + lambdawz(i,j,k)*xiWbarz(i,j,k)*xWbarz(m,i,j,k)*v(m,i,j,k)
                        cotnz2 = cotnz2 + lambdanz(i,j,k)*xiNbarz(i,j,k)*xNbarz(m,i,j,k)*v(m,i,j,k)
                    end do

                    left = -2*(cotwx2+cotnx2)/xedge/(xedge+ledge)
                    right = -2*(cotwx1+cotnx1)/xedge/(xedge+redge)
                    down = -2*(cotwy2+cotny2)/yedge/(yedge+dedge)
                    up = -2*(cotwy1+cotny1)/yedge/(yedge+uedge)
                    front = -2*(cotwz2+cotnz2)/zedge/(zedge+fedge)
                    back = -2*(cotwz1+cotnz1)/zedge/(zedge+bedge)

                    if((pcol == 1).and.(i == 1).and.(isDiriX(1,ylower+j-1,zlower+k-1) == 0)) then
                        sum1 = 0
                        sum2 = 0
                        do m = 1, Nc
                            sum1 = sum1 + xWbarx(m,i,j,k)*v(m,i,j,k)
                            sum2 = sum2 + xNbarx(m,i,j,k)*v(m,i,j,k)
                        end do
                        rhs_values(n) = rhs_values(n) + UwBdryX(1,ylower+j-1,zlower+k-1)*xiWbarx(i,j,k)*sum1/xedge + &!
                            UnBdryX(1,ylower+j-1,zlower+k-1)*xiNbarx(i,j,k)*sum2/xedge
                    else if((pcol == 1).and.(i == 1).and.(isDiriX(1,ylower+j-1,zlower+k-1) == 1)) then
                        rhs_values(n) = rhs_values(n) - left*Pw(i-1,j,k) + cotwx2*densiWbarx(i,j,k)* &!
                            gravX/xedge + cotnx2*densiNbarx(i,j,k)*gravX/xedge
                        values(r) = values(r) - left
                    else
                        values(r+1) = left
                        values(r) = values(r) - left
                        rhs_values(n) = rhs_values(n) + cotwx2*densiWbarx(i,j,k)*gravX/xedge + &!
                            cotnx2*densiNbarx(i,j,k)*gravX/xedge
                    end if

                    if((pcol == pncols).and.(i == localncols).and.(isDiriX(2,ylower+j-1,zlower+k-1) == 0)) then
                        sum1 = 0
                        sum2 = 0
                        do m = 1, Nc
                            sum1 = sum1 + xWbarx(m,i+1,j,k)*v(m,i,j,k)
                            sum2 = sum2 + xNbarx(m,i+1,j,k)*v(m,i,j,k)
                        end do
                        rhs_values(n) = rhs_values(n) - UwBdryX(2,ylower+j-1,zlower+k-1)*xiWbarx(i+1,j,k)*sum1/xedge - &!
                            UnBdryX(2,ylower+j-1,zlower+k-1)*xiNbarx(i+1,j,k)*sum2/xedge
                    else if((pcol == pncols).and.(i == localncols).and.(isDiriX(2,ylower+j-1,zlower+k-1) == 1)) then
                        rhs_values(n) = rhs_values(n) - right*Pw(i+1,j,k)-cotwx1*densiWbarx(i+1,j,k)* &!
                            gravX/xedge - cotnx1*densiNbarx(i+1,j,k)*gravX/xedge
                        values(r) = values(r) - right
                    else
                        values(r+2) = right
                        values(r) = values(r) - right
                        rhs_values(n) = rhs_values(n) - cotwx1*densiWbarx(i+1,j,k)*gravX/xedge - &!
                            cotnx1*densiNbarx(i+1,j,k)*gravX/xedge
                    end if

                    if((prow == 1).and.(j == 1).and.(isDiriY(xlower+i-1,1,zlower+k-1) == 0)) then
                        sum1 = 0
                        sum2 = 0
                        do m = 1, Nc
                            sum1 = sum1 + xWbary(m,i,j,k)*v(m,i,j,k)
                            sum2 = sum2 + xNbary(m,i,j,k)*v(m,i,j,k)
                        end do
                        rhs_values(n) = rhs_values(n) + UwBdryY(xlower+i-1,1,zlower+k-1)*xiWbary(i,j,k)*sum1/yedge + &!
                            UnBdryY(xlower+i-1,1,zlower+k-1)*xiNbary(i,j,k)*sum2/yedge
                    elseif((prow == 1).and.(j == 1).and.(isDiriY(xlower+i-1,1,zlower+k-1) == 1)) then
                        rhs_values(n) = rhs_values(n) - down*Pw(i,j-1,k) + cotwy2*densiWbary(i,j,k)* &!
                            gravY/yedge + cotny2*densiNbary(i,j,k)*gravY/yedge
                        values(r) = values(r) - down
                    else
                        values(r+3) = down
                        values(r) = values(r) - down
                        rhs_values(n) = rhs_values(n) + cotwy2*densiWbary(i,j,k)*gravY/yedge + &!
                            cotny2*densiNbary(i,j,k)*gravY/yedge
                    end if

                    if((prow == pnrows).and.(j == localnrows).and.(isDiriY(xlower+i-1,2,zlower+k-1) == 0)) then
                        sum1 = 0
                        sum2 = 0
                        do m = 1, Nc
                            sum1 = sum1 + xWbary(m,i,j+1,k)*v(m,i,j,k)
                            sum2 = sum2 + xNbary(m,i,j+1,k)*v(m,i,j,k)
                        end do
                        rhs_values(n) = rhs_values(n) - UwBdryY(xlower+i-1,2,zlower+k-1)*xiWbary(i,j+1,k)*sum1/yedge - &!
                            UnBdryY(xlower+i-1,2,zlower+k-1)*xiNbary(i,j+1,k)*sum2/yedge
                    else if((prow == pnrows).and.(j == localnrows).and.(isDiriY(xlower+i-1,2,zlower+k-1) == 1)) then
                        rhs_values(n) = rhs_values(n) - up*Pw(i,j+1,k) - cotwy1*densiWbary(i,j+1,k)*gravY&!
                            /yedge - cotny1*densiNbary(i,j+1,k)*gravY/yedge
                        values(r) = values(r) - up
                    else
                        values(r+4) = up
                        values(r) = values(r) - up
                        rhs_values(n) = rhs_values(n) - cotwy1*densiWbary(i,j+1,k)*gravY/yedge - &!
                            cotny1*densiNbary(i,j+1,k)*gravY/yedge
                    end if

                    if((play == 1).and.(k == 1).and.(isDiriZ(xlower+i-1,ylower+j-1,1) == 0)) then
                        sum1 = 0
                        sum2 = 0
                        do m = 1, Nc
                            sum1 = sum1 + xWbarz(m,i,j,k)*v(m,i,j,k)
                            sum2 = sum2 + xNbarz(m,i,j,k)*v(m,i,j,k)
                        end do
                        rhs_values(n) = rhs_values(n) + UwBdryZ(xlower+i-1,ylower+j-1,1)*xiWbarz(i,j,k)*sum1/zedge + &!
                            UnBdryZ(xlower+i-1,ylower+j-1,1)*xiNbarz(i,j,k)*sum2/zedge
                    elseif((play == 1).and.(k == 1).and.(isDiriZ(xlower+i-1,ylower+j-1,1) == 1)) then
                        rhs_values(n) = rhs_values(n) - front*Pw(i,j,k-1) + cotwz2*densiWbarz(i,j,k)* &!
                            gravZ/zedge + cotnz2*densiNbarz(i,j,k)*gravZ/zedge
                        values(r) = values(r) - front
                    else
                        values(r+5) = front
                        values(r) = values(r) - front
                        rhs_values(n) = rhs_values(n) + cotwz2*densiWbarz(i,j,k)*gravZ/zedge + &!
                            cotnz2*densiNbarz(i,j,k)*gravZ/zedge
                    end if

                    if((play == pnlays).and.(k == localnlays).and.(isDiriZ(xlower+i-1,ylower+j-1,2) == 0)) then
                        sum1 = 0
                        sum2 = 0
                        do m = 1, Nc
                            sum1 = sum1 + xWbarz(m,i,j,k+1)*v(m,i,j,k)
                            sum2 = sum2 + xNbarz(m,i,j,k+1)*v(m,i,j,k)
                        end do
                        rhs_values(n) = rhs_values(n) - UwBdryZ(xlower+i-1,ylower+j-1,2)*xiWbarz(i,j,k+1)*sum1/zedge - &!
                            UnBdryZ(xlower+i-1,ylower+j-1,2)*xiNbarz(i,j,k+1)*sum2/zedge
                    else if((play == pnlays).and.(k == localnlays).and.(isDiriZ(xlower+i-1,ylower+j-1,2) == 1)) then
                        rhs_values(n) = rhs_values(n) - back*Pw(i,j,k+1) - cotwz1*densiWbarz(i,j,k+1)*gravZ&!
                            /zedge - cotnz1*densiNbarz(i,j,k+1)*gravZ/zedge
                        values(r) = values(r) - back
                    else
                        values(r+6) = back
                        values(r) = values(r) - back
                        rhs_values(n) = rhs_values(n) - cotwz1*densiWbarz(i,j,k+1)*gravZ/zedge - &!
                            cotnz1*densiNbarz(i,j,k+1)*gravZ/zedge
                    end if

                    r = r + 7
                    n = n + 1

                end do
            end do
        end do

        call HYPRE_StructMatrixSetBoxValues(global_A, ilower, iupper, 7, stencil_indices, values, ierr)
        call HYPRE_StructMatrixAssemble(global_A, ierr)

        call HYPRE_StructVectorSetBoxValues(global_b, ilower, iupper, rhs_values, ierr)
        x_values(:) = initial_x_guess(:)
        call HYPRE_StructVectorSetBoxValues(global_x, ilower, iupper, x_values, ierr)
        call HYPRE_StructVectorAssemble(global_b, ierr)
        call HYPRE_StructVectorAssemble(global_x, ierr)

        solvertimestart = MPI_Wtime()
        call HYPRE_StructSMGSetup(solver, global_A, global_b, global_x, ierr)
        call HYPRE_StructSMGSolve(solver, global_A, global_b, global_x, ierr)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart
        call HYPRE_StructSMGGetNumIterations(solver, num_iter, ierr)
        if((ierr /= 0).or.(num_iter == 0)) then
            print *, 'The solver error.', ierr, num_iter
            call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
        end if

        call HYPRE_StructVectorGetBoxValues(global_x, ilower, iupper, x_values, ierr)

        n = 1
        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols
                    Pw(i,j,k) = x_values(n)
                    n = n + 1
                end do
            end do
        end do

        initial_x_guess(:) = x_values(:)

        deallocate(rhs_values)
        deallocate(x_values)
        deallocate(values)

        commtimestart = MPI_Wtime()

        call MPI_BUFFER_ATTACH(buffer,buffer_size,ierr)

        if(pcol /= 1) then
            allocate(Pwsent(localnrows*localnlays))
            n = 1
            do k = 1, localnlays
                do j = 1, localnrows
                    Pwsent(n) = Pw(1,j,k)
                    n = n+1
                end do
            end do
            call MPI_IBSEND(Pwsent, localnrows*localnlays, MPI_DOUBLE_PRECISION, myid-1, myid, &!
                MPI_COMM_WORLD, requestl, ierr)
            deallocate(Pwsent)
        end if

        if(pcol /= pncols) then
            allocate(Pwsent(localnrows*localnlays))
            n = 1
            do k = 1, localnlays
                do j = 1, localnrows
                    Pwsent(n) = Pw(localncols,j,k)
                    n = n+1
                end do
            end do
            call MPI_IBSEND(Pwsent, localnrows*localnlays, MPI_DOUBLE_PRECISION, myid+1, myid, &!
                MPI_COMM_WORLD, requestr, ierr)
            deallocate(Pwsent)
        end if

        if(prow /= 1) then
            allocate(Pwsent(localncols*localnlays))
            n = 1
            do k = 1, localnlays
                do i = 1, localncols
                    Pwsent(n) = Pw(i,1,k)
                    n = n+1
                end do
            end do
            call MPI_IBSEND(Pwsent, localncols*localnlays, MPI_DOUBLE_PRECISION, myid-pncols, &!
                myid, MPI_COMM_WORLD, requestd, ierr)
            deallocate(Pwsent)
        end if

        if(prow /= pnrows) then
            allocate(Pwsent(localncols*localnlays))
            n = 1
            do k = 1, localnlays
                do i = 1, localncols
                    Pwsent(n) = Pw(i,localnrows,k)
                    n = n+1
                end do
            end do
            call MPI_IBSEND(Pwsent, localncols*localnlays, MPI_DOUBLE_PRECISION, myid+pncols, &!
                myid, MPI_COMM_WORLD, requestu, ierr)
            deallocate(Pwsent)
        end if

        if(play /= 1) then
            allocate(Pwsent(localncols*localnrows))
            n = 1
            do j = 1, localnrows
                do i = 1, localncols
                    Pwsent(n) = Pw(i,j,1)
                    n = n+1
                end do
            end do
            call MPI_IBSEND(Pwsent, localncols*localnrows, MPI_DOUBLE_PRECISION, myid-pncols*pnrows, &!
                myid, MPI_COMM_WORLD, requestf, ierr)
            deallocate(Pwsent)
        end if

        if(play /= pnlays) then
            allocate(Pwsent(localncols*localnrows))
            n = 1
            do j = 1, localnrows
                do i = 1, localncols
                    Pwsent(n) = Pw(i,j,localnlays)
                    n = n+1
                end do
            end do
            call MPI_IBSEND(Pwsent, localncols*localnrows, MPI_DOUBLE_PRECISION, myid+pncols*pnrows, &!
                myid, MPI_COMM_WORLD, requestb, ierr)
            deallocate(Pwsent)
        end if

        if(pcol /= 1) then
            allocate(recvbuffer(localnrows*localnlays))
            call MPI_RECV(recvbuffer, localnrows*localnlays, MPI_DOUBLE_PRECISION, myid-1, myid-1, &!
                MPI_COMM_WORLD, status, ierr)
            n = 1
            do k = 1, localnlays
                do j = 1, localnrows
                    Pw(0,j,k) = recvbuffer(n)
                    n = n+1
                end do
            end do
            deallocate(recvbuffer)
        end if

        if(pcol /= pncols) then
            allocate(recvbuffer(localnrows*localnlays))
            call MPI_RECV(recvbuffer, localnrows*localnlays, MPI_DOUBLE_PRECISION, myid+1, myid+1, &!
                MPI_COMM_WORLD, status, ierr)
            n = 1
            do k = 1, localnlays
                do j = 1, localnrows
                    Pw(localncols+1,j,k) = recvbuffer(n)
                    n = n+1
                end do
            end do
            deallocate(recvbuffer)
        end if

        if(prow /= 1) then
            allocate(recvbuffer(localncols*localnlays))
            call MPI_RECV(recvbuffer, localncols*localnlays, MPI_DOUBLE_PRECISION, myid-pncols, myid-pncols, &!
                MPI_COMM_WORLD, status, ierr)
            n = 1
            do k = 1, localnlays
                do i = 1, localncols
                    Pw(i,0,k) = recvbuffer(n)
                    n = n+1
                end do
            end do
            deallocate(recvbuffer)
        end if

        if(prow /= pnrows) then
            allocate(recvbuffer(localncols*localnlays))
            call MPI_RECV(recvbuffer, localncols*localnlays, MPI_DOUBLE_PRECISION, myid+pncols, myid+pncols, &!
                MPI_COMM_WORLD, status, ierr)
            n = 1
            do k = 1, localnlays
                do i = 1, localncols
                    Pw(i,localnrows+1,k) = recvbuffer(n)
                    n = n+1
                end do
            end do
            deallocate(recvbuffer)
        end if

        if(play /= 1) then
            allocate(recvbuffer(localncols*localnrows))
            call MPI_RECV(recvbuffer, localncols*localnrows, MPI_DOUBLE_PRECISION, myid-pncols*pnrows, myid-pncols*pnrows, &!
                MPI_COMM_WORLD, status, ierr)
            n = 1
            do j = 1, localnrows
                do i = 1, localncols
                    Pw(i,j,0) = recvbuffer(n)
                    n = n+1
                end do
            end do
            deallocate(recvbuffer)
        end if

        if(play /= pnlays) then
            allocate(recvbuffer(localncols*localnrows))
            call MPI_RECV(recvbuffer, localncols*localnrows, MPI_DOUBLE_PRECISION, myid+pncols*pnrows, myid+pncols*pnrows, &!
                MPI_COMM_WORLD, status, ierr)
            n = 1
            do j = 1, localnrows
                do i = 1, localncols
                    Pw(i,j,localnlays+1) = recvbuffer(n)
                    n = n+1
                end do
            end do
            deallocate(recvbuffer)
        end if

        if(pcol /= 1) then
            call MPI_WAIT(requestl, status, ierr)
        end if
        if(pcol /= pncols) then
            call MPI_WAIT(requestr, status, ierr)
        end if
        if(prow /= 1) then
            call MPI_WAIT(requestd, status, ierr)
        end if
        if(prow /= pnrows) then
            call MPI_WAIT(requestu, status, ierr)
        end if
        if(play /= 1) then
            call MPI_WAIT(requestf, status, ierr)
        end if
        if(play /= pnlays) then
            call MPI_WAIT(requestb, status, ierr)
        end if

        call MPI_BUFFER_DETACH(buffer,buffer_size,ierr)

        commtimefinish = MPI_Wtime()
        commtime = commtime + commtimefinish - commtimestart

    end subroutine computePres

    subroutine computeMoleFrac()

        implicit none
       
        real(kind=8), dimension(:,:,:,:), pointer :: xic
        real(kind=8) :: left, right, down, up, front, back, div
        real(kind=8) :: xinew
        real(kind=8), dimension(:), pointer :: zsent, recvbuffer
        integer :: requestl, requestr, requestd, requestu, requestf, requestb
        integer :: status(MPI_STATUS_SIZE)
        integer :: i, j, k, m, n

        integer :: ierr, errorcode
        integer :: buffer_size = MAX_BUF
        real(kind=8) :: buffer(MAX_BUF)
        real(kind=8) :: commtimestart, commtimefinish

        allocate(xic(Nc, localncols, localnrows, localnlays))

        do k = 1, localnlays
            do j = 1, localnrows
                do i = 1, localncols
                    do m = 1, Nc
                        left = Uwx(i,j,k)*xiWbarx(i,j,k)*xWbarx(m,i,j,k) + Unx(i,j,k)*xiNbarx(i,j,k)*xNbarx(m,i,j,k)
                        right = Uwx(i+1,j,k)*xiWbarx(i+1,j,k)*xWbarx(m,i+1,j,k) + Unx(i+1,j,k)*xiNbarx(i+1,j,k) &!
                            *xNbarx(m,i+1,j,k)
                        down = Uwy(i,j,k)*xiWbary(i,j,k)*xWbary(m,i,j,k) + Uny(i,j,k)*xiNbary(i,j,k)*xNbary(m,i,j,k)
                        up = Uwy(i,j+1,k)*xiWbary(i,j+1,k)*xWbary(m,i,j+1,k) + Uny(i,j+1,k)*xiNbary(i,j+1,k) &!
                            *xNbary(m,i,j+1,k)
                        front = Uwz(i,j,k)*xiWbarz(i,j,k)*xWbarz(m,i,j,k) + Unz(i,j,k)*xiNbarz(i,j,k)*xNbarz(m,i,j,k)
                        back = Uwz(i,j,k+1)*xiWbarz(i,j,k+1)*xWbarz(m,i,j,k+1) + Unz(i,j,k+1)*xiNbarz(i,j,k+1) &!
                            *xNbarz(m,i,j,k+1)
                        div = (right-left)/(xs(xlower+i)-xs(xlower+i-1)) + (up-down)/(ys(ylower+j)-ys(ylower+j-1)) &!
                            + (back-front)/(zs(zlower+k)-zs(zlower+k-1))
                        xic(m,i,j,k) = (src(m,xlower+i-1,ylower+j-1,zlower+k-1)-div)*(timeEnd/nt)/ &!
                            poro(xlower+i-1,ylower+j-1,zlower+k-1) + z(m,i,j,k)*(Sw(i,j,k)*xiW(i,j,k)+(1-Sw(i,j,k))*xiN(i,j,k))
                        if(xic(m,i,j,k)<0) then
                            print *, 'Please tune the time step.', xic(m,i,j,k)
                            call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
                        end if
                    end do
                    xinew = 0
                    do m = 1, Nc
                        xinew = xinew + xic(m,i,j,k)
                    end do
                    z(1:Nc,i,j,k) = xic(1:Nc,i,j,k)/xinew
                    do m = 1, Nc
                        if(z(m,i,j,k) < 1.D-99) then
                            z(m,i,j,k) = 0.0
                        end if
                    end do
                end do
            end do
        end do

        commtimestart = MPI_Wtime()

        call MPI_BUFFER_ATTACH(buffer,buffer_size,ierr)

        if(pcol /= 1) then
            allocate(zsent(Nc*localnrows*localnlays))
            n = 1
            do k = 1, localnlays
                do j = 1, localnrows
                    do m = 1, Nc
                        zsent(n) = z(m,1,j,k)
                        n = n + 1
                    end do
                end do
            end do
            call MPI_IBSEND(zsent, Nc*localnrows*localnlays, MPI_DOUBLE_PRECISION, myid-1, myid, &!
                MPI_COMM_WORLD, requestl, ierr)
            deallocate(zsent)
        end if

        if(pcol /= pncols) then
            allocate(zsent(Nc*localnrows*localnlays))
            n = 1
            do k = 1, localnlays
                do j = 1, localnrows
                    do m = 1, Nc
                        zsent(n) = z(m,localncols,j,k)
                        n = n + 1
                    end do
                end do
            end do
            call MPI_IBSEND(zsent, Nc*localnrows*localnlays, MPI_DOUBLE_PRECISION, myid+1, myid, &!
                MPI_COMM_WORLD, requestr, ierr)
            deallocate(zsent)
        end if

        if(prow /= 1) then
            allocate(zsent(Nc*localncols*localnlays))
            n = 1
            do k = 1, localnlays
                do i = 1, localncols
                    do m = 1, Nc
                        zsent(n) = z(m,i,1,k)
                        n = n + 1
                    end do
                end do
            end do
            call MPI_IBSEND(zsent, Nc*localncols*localnlays, MPI_DOUBLE_PRECISION, myid-pncols, &!
                myid, MPI_COMM_WORLD, requestd, ierr)
            deallocate(zsent)
        end if

        if(prow /= pnrows) then
            allocate(zsent(Nc*localncols*localnlays))
            n = 1
            do k = 1, localnlays
                do i = 1, localncols
                    do m = 1, Nc
                        zsent(n) = z(m,i,localnrows,k)
                        n = n + 1
                    end do
                end do
            end do
            call MPI_IBSEND(zsent, Nc*localncols*localnlays, MPI_DOUBLE_PRECISION, myid+pncols, &!
                myid, MPI_COMM_WORLD, requestu, ierr)
            deallocate(zsent)
        end if

        if(play /= 1) then
            allocate(zsent(Nc*localncols*localnrows))
            n = 1
            do j = 1, localnrows
                do i = 1, localncols
                    do m = 1, Nc
                        zsent(n) = z(m,i,j,1)
                        n = n + 1
                    end do
                end do
            end do
            call MPI_IBSEND(zsent, Nc*localncols*localnrows, MPI_DOUBLE_PRECISION, myid-pncols*pnrows, &!
                myid, MPI_COMM_WORLD, requestf, ierr)
            deallocate(zsent)
        end if

        if(play /= pnlays) then
            allocate(zsent(Nc*localncols*localnrows))
            n = 1
            do j = 1, localnrows
                do i = 1, localncols
                    do m = 1, Nc
                        zsent(n) = z(m,i,j,localnlays)
                        n = n + 1
                    end do
                end do
            end do
            call MPI_IBSEND(zsent, Nc*localncols*localnrows, MPI_DOUBLE_PRECISION, myid+pncols*pnrows, &!
                myid, MPI_COMM_WORLD, requestb, ierr)
            deallocate(zsent)
        end if

        if(pcol /= 1) then
            allocate(recvbuffer(Nc*localnrows*localnlays))
            call MPI_RECV(recvbuffer, Nc*localnrows*localnlays, MPI_DOUBLE_PRECISION, myid-1, myid-1, &!
                MPI_COMM_WORLD, status, ierr)
            n = 1
            do k = 1, localnlays
                do j = 1, localnrows
                    do m = 1, Nc
                        z(m,0,j,k) = recvbuffer(n)
                        n = n + 1
                    end do
                end do
            end do
            deallocate(recvbuffer)
        end if

        if(pcol /= pncols) then
            allocate(recvbuffer(Nc*localnrows*localnlays))
            call MPI_RECV(recvbuffer, Nc*localnrows*localnlays, MPI_DOUBLE_PRECISION, myid+1, myid+1, &!
                MPI_COMM_WORLD, status, ierr)
            n = 1
            do k = 1, localnlays
                do j = 1, localnrows
                    do m = 1, Nc
                        z(m,localncols+1,j,k) = recvbuffer(n)
                        n = n + 1
                    end do
                end do
            end do
            deallocate(recvbuffer)
        end if

        if(prow /= 1) then
            allocate(recvbuffer(Nc*localncols*localnlays))
            call MPI_RECV(recvbuffer, Nc*localncols*localnlays, MPI_DOUBLE_PRECISION, myid-pncols, myid-pncols, &!
                MPI_COMM_WORLD, status, ierr)
            n = 1
            do k = 1, localnlays
                do i = 1, localncols
                    do m = 1, Nc
                        z(m,i,0,k) = recvbuffer(n)
                        n = n + 1
                    end do
                end do
            end do
            deallocate(recvbuffer)
        end if

        if(prow /= pnrows) then
            allocate(recvbuffer(Nc*localncols*localnlays))
            call MPI_RECV(recvbuffer, Nc*localncols*localnlays, MPI_DOUBLE_PRECISION, myid+pncols, myid+pncols, &!
                MPI_COMM_WORLD, status, ierr)
            n = 1
            do k = 1, localnlays
                do i = 1, localncols
                    do m = 1, Nc
                        z(m,i,localnrows+1,k) = recvbuffer(n)
                        n = n + 1
                    end do
                end do
            end do
            deallocate(recvbuffer)
        end if

        if(play /= 1) then
            allocate(recvbuffer(Nc*localncols*localnrows))
            call MPI_RECV(recvbuffer, Nc*localncols*localnrows, MPI_DOUBLE_PRECISION, myid-pncols*pnrows, &!
                myid-pncols*pnrows, MPI_COMM_WORLD, status, ierr)
            n = 1
            do j = 1, localnrows
                do i = 1, localncols
                    do m = 1, Nc
                        z(m,i,j,0) = recvbuffer(n)
                        n = n + 1
                    end do
                end do
            end do
            deallocate(recvbuffer)
        end if

        if(play /= pnlays) then
            allocate(recvbuffer(Nc*localncols*localnrows))
            call MPI_RECV(recvbuffer, Nc*localncols*localnrows, MPI_DOUBLE_PRECISION, myid+pncols*pnrows, &!
                myid+pncols*pnrows, MPI_COMM_WORLD, status, ierr)
            n = 1
            do j = 1, localnrows
                do i = 1, localncols
                    do m = 1, Nc
                        z(m,i,j,localnlays+1) = recvbuffer(n)
                        n = n + 1
                    end do
                end do
            end do
            deallocate(recvbuffer)
        end if

        if(pcol /= 1) then
            call MPI_WAIT(requestl, status, ierr)
        end if
        if(pcol /= pncols) then
            call MPI_WAIT(requestr, status, ierr)
        end if
        if(prow /= 1) then
            call MPI_WAIT(requestd, status, ierr)
        end if
        if(prow /= pnrows) then
            call MPI_WAIT(requestu, status, ierr)
        end if
        if(play /= 1) then
            call MPI_WAIT(requestf, status, ierr)
        end if
        if(play /= pnlays) then
            call MPI_WAIT(requestb, status, ierr)
        end if

        call MPI_BUFFER_DETACH(buffer,buffer_size,ierr)

        commtimefinish = MPI_Wtime()
        commtime = commtime + commtimefinish - commtimestart

        deallocate(xic)

    end subroutine computeMoleFrac

    subroutine computeVel()

        implicit none
      
        integer :: i, j, k, ierr

        ! compute the total new velocities in x direction
        if(pncols==1) then
            do k = 1, localnlays
                do j = 1, localnrows
                    if(isDiriX(1,ylower+j-1,zlower+k-1) == 1) then
                        Uwx(1,j,k) = -lambdawx(1,j,k)*((Pw(1,j,k)-Pw(0,j,k))*2/(xs(2)-xs(1)) - densiWbarx(1,j,k)*gravX)
                        Unx(1,j,k) = -lambdanx(1,j,k)*((Pw(1,j,k)-Pw(0,j,k))*2/(xs(2)-xs(1)) - densiNbarx(1,j,k)*gravX)
                    else
                        Uwx(1,j,k) = UwBdryX(1,ylower+j-1,zlower+k-1)
                        Unx(1,j,k) = UnBdryX(1,ylower+j-1,zlower+k-1)
                    end if
                end do
            end do
            do k = 1, localnlays
                do j = 1, localnrows
                    if(isDiriX(2,ylower+j-1,zlower+k-1) == 1) then
                        Uwx(localncols+1,j,k) = -lambdawx(localncols+1,j,k)*((Pw(localncols+1,j,k)-Pw(localncols,j,k)) &!
                            *2/(xs(nx+1)-xs(nx)) - densiWbarx(localncols+1,j,k)*gravX)
                        Unx(localncols+1,j,k) = -lambdanx(localncols+1,j,k)*((Pw(localncols+1,j,k)-Pw(localncols,j,k)) &!
                            *2/(xs(nx+1)-xs(nx)) - densiNbarx(localncols+1,j,k)*gravX)
                        if(Uwx(localncols+1,j,k) < 0) then
                            Uwx(localncols+1,j,k) = 0
                        end if
                        if(Unx(localncols+1,j,k) < 0) then
                            Unx(localncols+1,j,k) = 0
                        end if
                    else
                        Uwx(localncols+1,j,k) = UwBdryX(2,ylower+j-1,zlower+k-1)
                        Unx(localncols+1,j,k) = UnBdryX(2,ylower+j-1,zlower+k-1)
                    end if
                end do
            end do
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 2, localncols
                        Uwx(i,j,k) = -lambdawx(i,j,k)*((Pw(i,j,k)-Pw(i-1,j,k))*2/(xs(xlower+i)-xs(xlower+i-2)) &!
                            - densiWbarx(i,j,k)*gravX)
                        Unx(i,j,k) = -lambdanx(i,j,k)*((Pw(i,j,k)-Pw(i-1,j,k))*2/(xs(xlower+i)-xs(xlower+i-2)) &!
                            - densiNbarx(i,j,k)*gravX)
                    end do
                end do
            end do
        elseif(pcol==1) then
            do k = 1, localnlays
                do j = 1, localnrows
                    if(isDiriX(1,ylower+j-1,zlower+k-1) == 1) then
                        Uwx(1,j,k) = -lambdawx(1,j,k)*((Pw(1,j,k)-Pw(0,j,k))*2/(xs(2)-xs(1)) - densiWbarx(1,j,k)*gravX)
                        Unx(1,j,k) = -lambdanx(1,j,k)*((Pw(1,j,k)-Pw(0,j,k))*2/(xs(2)-xs(1)) - densiNbarx(1,j,k)*gravX)
                    else
                        Uwx(1,j,k) = UwBdryX(1,ylower+j-1,zlower+k-1)
                        Unx(1,j,k) = UnBdryX(1,ylower+j-1,zlower+k-1)
                    end if
                end do
            end do
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 2, localncols+1
                        Uwx(i,j,k) = -lambdawx(i,j,k)*((Pw(i,j,k)-Pw(i-1,j,k))*2/(xs(xlower+i)-xs(xlower+i-2)) &!
                            - densiWbarx(i,j,k)*gravX)
                        Unx(i,j,k) = -lambdanx(i,j,k)*((Pw(i,j,k)-Pw(i-1,j,k))*2/(xs(xlower+i)-xs(xlower+i-2)) &!
                            - densiNbarx(i,j,k)*gravX)
                    end do
                end do
            end do
        elseif(pcol==pncols) then
            do k = 1, localnlays
                do j = 1, localnrows
                    if(isDiriX(2,ylower+j-1,zlower+k-1) == 1) then
                        Uwx(localncols+1,j,k) = -lambdawx(localncols+1,j,k)*((Pw(localncols+1,j,k)-Pw(localncols,j,k)) &!
                            *2/(xs(nx+1)-xs(nx)) - densiWbarx(localncols+1,j,k)*gravX)
                        Unx(localncols+1,j,k) = -lambdanx(localncols+1,j,k)*((Pw(localncols+1,j,k)-Pw(localncols,j,k)) &!
                            *2/(xs(nx+1)-xs(nx)) - densiNbarx(localncols+1,j,k)*gravX)
                        if(Uwx(localncols+1,j,k) < 0) then
                            Uwx(localncols+1,j,k) = 0
                        end if
                        if(Unx(localncols+1,j,k) < 0) then
                            Unx(localncols+1,j,k) = 0
                        end if
                    else
                        Uwx(localncols+1,j,k) = UwBdryX(2,ylower+j-1,zlower+k-1)
                        Unx(localncols+1,j,k) = UnBdryX(2,ylower+j-1,zlower+k-1)
                    end if
                end do
            end do
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        Uwx(i,j,k) = -lambdawx(i,j,k)*((Pw(i,j,k)-Pw(i-1,j,k))*2/(xs(xlower+i)-xs(xlower+i-2)) &!
                            - densiWbarx(i,j,k)*gravX)
                        Unx(i,j,k) = -lambdanx(i,j,k)*((Pw(i,j,k)-Pw(i-1,j,k))*2/(xs(xlower+i)-xs(xlower+i-2)) &!
                            - densiNbarx(i,j,k)*gravX)
                    end do
                end do
            end do
        else
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols+1
                        Uwx(i,j,k) = -lambdawx(i,j,k)*((Pw(i,j,k)-Pw(i-1,j,k))*2/(xs(xlower+i)-xs(xlower+i-2)) &!
                            - densiWbarx(i,j,k)*gravX)
                        Unx(i,j,k) = -lambdanx(i,j,k)*((Pw(i,j,k)-Pw(i-1,j,k))*2/(xs(xlower+i)-xs(xlower+i-2)) &!
                            - densiNbarx(i,j,k)*gravX)
                    end do
                end do
            end do
        end if

        ! compute the total new velocities in y direction
        if(pnrows==1) then
            do k = 1, localnlays
                do i = 1, localncols
                    if(isDiriY(xlower+i-1,1,zlower+k-1) == 1) then
                        Uwy(i,1,k) = -lambdawy(i,1,k)*((Pw(i,1,k)-Pw(i,0,k))*2/(ys(2)-ys(1)) - densiWbary(i,1,k)*gravY)
                        Uny(i,1,k) = -lambdany(i,1,k)*((Pw(i,1,k)-Pw(i,0,k))*2/(ys(2)-ys(1)) - densiNbary(i,1,k)*gravY)
                    else
                        Uwy(i,1,k) = UwBdryY(xlower+i-1,1,zlower+k-1)
                        Uny(i,1,k) = UnBdryY(xlower+i-1,1,zlower+k-1)
                    end if
                end do
            end do
            do k = 1, localnlays
                do i = 1, localncols
                    if(isDiriY(xlower+i-1,2,zlower+k-1) == 1) then
                        Uwy(i,localnrows+1,k) = -lambdawy(i,localnrows+1,k)*((Pw(i,localnrows+1,k)-Pw(i,localnrows,k)) &!
                            *2/(ys(ny+1)-ys(ny)) - densiWbary(i,localnrows+1,k)*gravY)
                        Uny(i,localnrows+1,k) = -lambdany(i,localnrows+1,k)*((Pw(i,localnrows+1,k)-Pw(i,localnrows,k)) &!
                            *2/(ys(ny+1)-ys(ny)) - densiNbary(i,localnrows+1,k)*gravY)
                    else
                        Uwy(i,localnrows+1,k) = UwBdryY(xlower+i-1,2,zlower+k-1)
                        Uny(i,localnrows+1,k) = UnBdryY(xlower+i-1,2,zlower+k-1)
                    end if
                end do
            end do
            do k = 1, localnlays
                do j = 2, localnrows
                    do i = 1, localncols
                        Uwy(i,j,k) = -lambdawy(i,j,k)*((Pw(i,j,k)-Pw(i,j-1,k))*2/(ys(ylower+j)-ys(ylower+j-2)) &!
                            - densiWbary(i,j,k)*gravY)
                        Uny(i,j,k) = -lambdany(i,j,k)*((Pw(i,j,k)-Pw(i,j-1,k))*2/(ys(ylower+j)-ys(ylower+j-2)) &!
                            - densiNbary(i,j,k)*gravY)
                    end do
                end do
            end do
        elseif(prow==1) then
            do k = 1, localnlays
                do i = 1, localncols
                    if(isDiriY(xlower+i-1,1,zlower+k-1) == 1) then
                        Uwy(i,1,k) = -lambdawy(i,1,k)*((Pw(i,1,k)-Pw(i,0,k))*2/(ys(2)-ys(1)) - densiWbary(i,1,k)*gravY)
                        Uny(i,1,k) = -lambdany(i,1,k)*((Pw(i,1,k)-Pw(i,0,k))*2/(ys(2)-ys(1)) - densiNbary(i,1,k)*gravY)
                    else
                        Uwy(i,1,k) = UwBdryY(xlower+i-1,1,zlower+k-1)
                        Uny(i,1,k) = UnBdryY(xlower+i-1,1,zlower+k-1)
                    end if
                end do
            end do
            do k = 1, localnlays
                do j = 2, localnrows+1
                    do i = 1, localncols
                        Uwy(i,j,k) = -lambdawy(i,j,k)*((Pw(i,j,k)-Pw(i,j-1,k))*2/(ys(ylower+j)-ys(ylower+j-2)) &!
                            - densiWbary(i,j,k)*gravY)
                        Uny(i,j,k) = -lambdany(i,j,k)*((Pw(i,j,k)-Pw(i,j-1,k))*2/(ys(ylower+j)-ys(ylower+j-2)) &!
                            - densiNbary(i,j,k)*gravY)
                    end do
                end do
            end do
        elseif(prow==pnrows) then
            do k = 1, localnlays
                do i = 1, localncols
                    if(isDiriY(xlower+i-1,2,zlower+k-1) == 1) then
                        Uwy(i,localnrows+1,k) = -lambdawy(i,localnrows+1,k)*((Pw(i,localnrows+1,k)-Pw(i,localnrows,k)) &!
                            *2/(ys(ny+1)-ys(ny)) - densiWbary(i,localnrows+1,k)*gravY)
                        Uny(i,localnrows+1,k) = -lambdany(i,localnrows+1,k)*((Pw(i,localnrows+1,k)-Pw(i,localnrows,k)) &!
                            *2/(ys(ny+1)-ys(ny)) - densiNbary(i,localnrows+1,k)*gravY)
                    else
                        Uwy(i,localnrows+1,k) = UwBdryY(xlower+i-1,2,zlower+k-1)
                        Uny(i,localnrows+1,k) = UnBdryY(xlower+i-1,2,zlower+k-1)
                    end if
                end do
            end do
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        Uwy(i,j,k) = -lambdawy(i,j,k)*((Pw(i,j,k)-Pw(i,j-1,k))*2/(ys(ylower+j)-ys(ylower+j-2)) &!
                            - densiWbary(i,j,k)*gravY)
                        Uny(i,j,k) = -lambdany(i,j,k)*((Pw(i,j,k)-Pw(i,j-1,k))*2/(ys(ylower+j)-ys(ylower+j-2)) &!
                            - densiNbary(i,j,k)*gravY)
                    end do
                end do
            end do
        else
            do k = 1, localnlays
                do j = 1, localnrows+1
                    do i = 1, localncols
                        Uwy(i,j,k) = -lambdawy(i,j,k)*((Pw(i,j,k)-Pw(i,j-1,k))*2/(ys(ylower+j)-ys(ylower+j-2)) &!
                            - densiWbary(i,j,k)*gravY)
                        Uny(i,j,k) = -lambdany(i,j,k)*((Pw(i,j,k)-Pw(i,j-1,k))*2/(ys(ylower+j)-ys(ylower+j-2)) &!
                            - densiNbary(i,j,k)*gravY)
                    end do
                end do
            end do
        end if

        ! compute the total new velocities in z direction
        if(pnlays==1) then
            do j = 1, localnrows
                do i = 1, localncols
                    if(isDiriZ(xlower+i-1,ylower+j-1,1) == 1) then
                        Uwz(i,j,1) = -lambdawz(i,j,1)*((Pw(i,j,1)-Pw(i,j,0))*2/(zs(2)-zs(1)) - densiWbarz(i,j,1)*gravZ)
                        Unz(i,j,1) = -lambdanz(i,j,1)*((Pw(i,j,1)-Pw(i,j,0))*2/(zs(2)-zs(1)) - densiNbarz(i,j,1)*gravZ)
                    else
                        Uwz(i,j,1) = UwBdryZ(xlower+i-1,ylower+j-1,1)
                        Unz(i,j,1) = UnBdryZ(xlower+i-1,ylower+j-1,1)
                    end if
                end do
            end do
            do j = 1, localnrows
                do i = 1, localncols
                    if(isDiriZ(xlower+i-1,ylower+j-1,2) == 1) then
                        Uwz(i,j,localnlays+1) = -lambdawz(i,j,localnlays+1)*((Pw(i,j,localnlays+1)-Pw(i,j,localnlays)) &!
                            *2/(zs(nz+1)-zs(nz)) - densiWbarz(i,j,localnlays+1)*gravZ)
                        Unz(i,j,localnlays+1) = -lambdanz(i,j,localnlays+1)*((Pw(i,j,localnlays+1)-Pw(i,j,localnlays)) &!
                            *2/(zs(nz+1)-zs(nz)) - densiNbarz(i,j,localnlays+1)*gravZ)
                    else
                        Uwz(i,j,localnlays+1) = UwBdryZ(xlower+i-1,ylower+j-1,2)
                        Unz(i,j,localnlays+1) = UnBdryZ(xlower+i-1,ylower+j-1,2)
                    end if
                end do
            end do
            do k = 2, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        Uwz(i,j,k) = -lambdawz(i,j,k)*((Pw(i,j,k)-Pw(i,j,k-1))*2/(zs(zlower+k)-zs(zlower+k-2)) &!
                            - densiWbarz(i,j,k)*gravZ)
                        Unz(i,j,k) = -lambdanz(i,j,k)*((Pw(i,j,k)-Pw(i,j,k-1))*2/(zs(zlower+k)-zs(zlower+k-2)) &!
                            - densiNbarz(i,j,k)*gravZ)
                    end do
                end do
            end do
        elseif(play==1) then
            do j = 1, localnrows
                do i = 1, localncols
                    if(isDiriZ(xlower+i-1,ylower+j-1,1) == 1) then
                        Uwz(i,j,1) = -lambdawz(i,j,1)*((Pw(i,j,1)-Pw(i,j,0))*2/(zs(2)-zs(1)) - densiWbarz(i,j,1)*gravZ)
                        Unz(i,j,1) = -lambdanz(i,j,1)*((Pw(i,j,1)-Pw(i,j,0))*2/(zs(2)-zs(1)) - densiNbarz(i,j,1)*gravZ)
                    else
                        Uwz(i,j,1) = UwBdryZ(xlower+i-1,ylower+j-1,1)
                        Unz(i,j,1) = UnBdryZ(xlower+i-1,ylower+j-1,1)
                    end if
                end do
            end do
            do k = 2, localnlays+1
                do j = 1, localnrows
                    do i = 1, localncols
                        Uwz(i,j,k) = -lambdawz(i,j,k)*((Pw(i,j,k)-Pw(i,j,k-1))*2/(zs(zlower+k)-zs(zlower+k-2)) &!
                            - densiWbarz(i,j,k)*gravZ)
                        Unz(i,j,k) = -lambdanz(i,j,k)*((Pw(i,j,k)-Pw(i,j,k-1))*2/(zs(zlower+k)-zs(zlower+k-2)) &!
                            - densiNbarz(i,j,k)*gravZ)
                    end do
                end do
            end do
        elseif(play==pnlays) then
            do j = 1, localnrows
                do i = 1, localncols
                    if(isDiriZ(xlower+i-1,ylower+j-1,2) == 1) then
                        Uwz(i,j,localnlays+1) = -lambdawz(i,j,localnlays+1)*((Pw(i,j,localnlays+1)-Pw(i,j,localnlays)) &!
                            *2/(zs(nz+1)-zs(nz)) - densiWbarz(i,j,localnlays+1)*gravZ)
                        Unz(i,j,localnlays+1) = -lambdanz(i,j,localnlays+1)*((Pw(i,j,localnlays+1)-Pw(i,j,localnlays)) &!
                            *2/(zs(nz+1)-zs(nz)) - densiNbarz(i,j,localnlays+1)*gravZ)
                    else
                        Uwz(i,j,localnlays+1) = UwBdryZ(xlower+i-1,ylower+j-1,2)
                        Unz(i,j,localnlays+1) = UnBdryZ(xlower+i-1,ylower+j-1,2)
                    end if
                end do
            end do
            do k = 1, localnlays
                do j = 1, localnrows
                    do i = 1, localncols
                        Uwz(i,j,k) = -lambdawz(i,j,k)*((Pw(i,j,k)-Pw(i,j,k-1))*2/(zs(zlower+k)-zs(zlower+k-2)) &!
                            - densiWbarz(i,j,k)*gravZ)
                        Unz(i,j,k) = -lambdanz(i,j,k)*((Pw(i,j,k)-Pw(i,j,k-1))*2/(zs(zlower+k)-zs(zlower+k-2)) &!
                            - densiNbarz(i,j,k)*gravZ)
                    end do
                end do
            end do
        else
            do k = 1, localnlays+1
                do j = 1, localnrows
                    do i = 1, localncols
                        Uwz(i,j,k) = -lambdawz(i,j,k)*((Pw(i,j,k)-Pw(i,j,k-1))*2/(zs(zlower+k)-zs(zlower+k-2)) &!
                            - densiWbarz(i,j,k)*gravZ)
                        Unz(i,j,k) = -lambdanz(i,j,k)*((Pw(i,j,k)-Pw(i,j,k-1))*2/(zs(zlower+k)-zs(zlower+k-2)) &!
                            - densiNbarz(i,j,k)*gravZ)
                    end do
                end do
            end do
        end if

    end subroutine computeVel

    subroutine finalize()

        implicit none
       
        integer :: ierr
        real(kind=8) :: timefinish

        deallocate(xs)
        deallocate(ys)
        deallocate(zs)
        deallocate(ts)
        deallocate(Kxx)
        deallocate(Kyy)
        deallocate(Kzz)
        deallocate(poro)
        deallocate(src)
        deallocate(isDiriX)
        deallocate(isDiriY)
        deallocate(isDiriZ)
        deallocate(PwBdryX)
        deallocate(PwBdryY)
        deallocate(PwBdryZ)
        deallocate(PwInit)
        deallocate(zBdryX)
        deallocate(zBdryY)
        deallocate(zBdryZ)
        deallocate(zInit)
        deallocate(UwBdryX)
        deallocate(UwBdryY)
        deallocate(UwBdryZ)
        deallocate(UnBdryX)
        deallocate(UnBdryY)
        deallocate(UnBdryZ)
        deallocate(ct)
        deallocate(cp)
        deallocate(af)
        deallocate(mw)
        deallocate(cv)
        deallocate(psatA)
        deallocate(psatB)
        deallocate(psatC)
        deallocate(delta)

        deallocate(Pw)
        deallocate(Uwx)
        deallocate(Uwy)
        deallocate(Uwz)
        deallocate(Unx)
        deallocate(Uny)
        deallocate(Unz)
        deallocate(Sw)
        deallocate(lambdawx)
        deallocate(lambdawy)
        deallocate(lambdawz)
        deallocate(lambdanx)
        deallocate(lambdany)
        deallocate(lambdanz)
        deallocate(Kxxbar)
        deallocate(Kyybar)
        deallocate(Kzzbar)
        deallocate(z)
        deallocate(densiW)
        deallocate(densiN)
        deallocate(densiWbarx)
        deallocate(densiWbary)
        deallocate(densiWbarz)
        deallocate(densiNbarx)
        deallocate(densiNbary)
        deallocate(densiNbarz)
        deallocate(xW)
        deallocate(xN)
        deallocate(xWbarx)
        deallocate(xWbary)
        deallocate(xWbarz)
        deallocate(xNbarx)
        deallocate(xNbary)
        deallocate(xNbarz)
        deallocate(xiW)
        deallocate(xiN)
        deallocate(xiWbarx)
        deallocate(xiWbary)
        deallocate(xiWbarz)
        deallocate(xiNbarx)
        deallocate(xiNbary)
        deallocate(xiNbarz)
        deallocate(v)
        deallocate(Cf)
        deallocate(viscW)
        deallocate(viscN)

#ifdef SPARSE
        deallocate(xtable)
        deallocate(ytable)
        deallocate(xiLtable)
        deallocate(xiGtable)
        deallocate(rhoLtable)
        deallocate(rhoGtable)
        deallocate(sLtable)
        deallocate(vtable)
        deallocate(Cftable)
        deallocate(isWtable)
        deallocate(isNtable)
#endif

        deallocate(initial_x_guess)

        call HYPRE_StructGridDestroy(grid, ierr)
        call HYPRE_StructStencilDestroy(stencil, ierr)
        call HYPRE_StructMatrixDestroy(global_A, ierr)
        call HYPRE_StructVectorDestroy(global_b, ierr)
        call HYPRE_StructVectorDestroy(global_x, ierr)
        call HYPRE_StructSMGDestroy(solver, ierr)

        if(myid == 0) then
            close(20)
            close(30)
        end if

#ifdef SPARSE
        close(120)
        close(130)
        close(150)
        close(160)
        close(180)
        close(190)
        close(200)
        close(210)
        close(220)
        close(230)
        close(240)
        close(260)
        close(270)
        close(280)

        if(Nc == 3) then
            close(140)
            close(170)
            close(250)
        end if
#endif

        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        timefinish = MPI_Wtime()
        print *, 'The elapsed time = ', timefinish-timestart, ' seconds.'
        print *, 'The flash time = ', flashtime
        print *, 'The solver time = ', solvertime
        print *, 'The communication time = ', commtime

        call MPI_Finalize(ierr)

    end subroutine finalize

end module RST_compositionalTwoPhaseFlow
