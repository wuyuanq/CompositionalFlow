
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
        character(len=50), dimension(:), pointer :: fadtxt
        character :: charm
        integer :: indexl, indexr, indexu, indexd
        integer :: ierr, errorcode
        integer :: i, j, k, m

        call MPI_INIT(ierr)

        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        timestart = MPI_Wtime()
        flashtime = 0.0
        solvertime = 0.0
        commtime = 0.0

        pncols = modelCase%pncols
        pnrows = modelCase%pnrows
        Nc = modelCase%Nc
        Temp = modelCase%Temp
        Lx = modelCase%Lx
        Ly = modelCase%Ly
        timeEnd = modelCase%timeEnd
        nx = modelCase%nx
        ny = modelCase%ny
        nt = modelCase%nt
        gravX = modelCase%gravX
        gravY = modelCase%gravY

        allocate(xs(nx+1))
        allocate(ys(ny+1))
        allocate(ts(nt+1))
        allocate(Kxx(nx,ny))
        allocate(Kyy(nx,ny))
        allocate(poro(nx,ny))
        allocate(src(Nc,nx,ny))
        allocate(isDiriX(2,ny))
        allocate(isDiriY(nx,2))
        allocate(PwBdryX(2,ny))
        allocate(PwBdryY(nx,2))
        allocate(PwInit(nx,ny))
        allocate(zBdryX(Nc,2,ny))
        allocate(zBdryY(Nc,nx,2))
        allocate(zInit(Nc,nx,ny))
        allocate(UwBdryX(2,ny))
        allocate(UwBdryY(nx,2))
        allocate(UnBdryX(2,ny))
        allocate(UnBdryY(nx,2))
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
        ts = modelCase%ts
        Kxx = modelCase%Kxx
        Kyy = modelCase%Kyy
        poro = modelCase%poro
        src = modelCase%src
        isDiriX = modelCase%isDiriX
        isDiriY = modelCase%isDiriY
        PwBdryX = modelCase%PwBdryX
        PwBdryY = modelCase%PwBdryY
        PwInit = modelCase%PwInit
        zBdryX = modelCase%zBdryX
        zBdryY = modelCase%zBdryY
        zInit = modelCase%zInit
        UwBdryX = modelCase%UwBdryX
        UwBdryY = modelCase%UwBdryY
        UnBdryX = modelCase%UnBdryX
        UnBdryY = modelCase%UnBdryY
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

        UwBdryX(1, 1:ny) = -UwBdryX(1, 1:ny)
        UwBdryY(1:nx, 1) = -UwBdryY(1:nx, 1)
        UnBdryX(1, 1:ny) = -UnBdryX(1, 1:ny)
        UnBdryY(1:nx, 1) = -UnBdryY(1:nx, 1)

        call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

        localncols = nx/pncols
        localnrows = ny/pnrows

        pcol = myid/pnrows+1
        prow = mod(myid,pnrows)+1

        xlower = (pcol-1)*localncols+1
        xupper = pcol*localncols
        ylower = (prow-1)*localnrows+1
        yupper = prow*localnrows

        allocate(Pw(0:localnrows+1, 0:localncols+1))

        allocate(Uwx(localnrows, localncols+1))
        Uwx = 0
        allocate(Uwy(localnrows+1, localncols))
        Uwy = 0
        allocate(Unx(localnrows, localncols+1))
        Unx = 0
        allocate(Uny(localnrows+1, localncols))
        Uny = 0

        allocate(Sw(0:localnrows+1, 0:localncols+1))

        allocate(lambdawx(localnrows, localncols+1))
        allocate(lambdawy(localnrows+1, localncols))
        allocate(lambdanx(localnrows, localncols+1))
        allocate(lambdany(localnrows+1, localncols))

        allocate(Kxxbar(localnrows, localncols+1))
        allocate(Kyybar(localnrows+1, localncols))

        allocate(z(Nc, 0:localnrows+1, 0:localncols+1))

        allocate(densiW(0:localnrows+1, 0:localncols+1))
        allocate(densiN(0:localnrows+1, 0:localncols+1))
        allocate(densiWbarx(localnrows, localncols+1))
        allocate(densiWbary(localnrows+1, localncols))
        allocate(densiNbarx(localnrows, localncols+1))
        allocate(densiNbary(localnrows+1, localncols))

        allocate(xW(Nc, 0:localnrows+1, 0:localncols+1))
        allocate(xN(Nc, 0:localnrows+1, 0:localncols+1))
        allocate(xWbarx(Nc, localnrows, localncols+1))
        allocate(xWbary(Nc, localnrows+1, localncols))
        allocate(xNbarx(Nc, localnrows, localncols+1))
        allocate(xNbary(Nc, localnrows+1, localncols))

        allocate(xiW(0:localnrows+1, 0:localncols+1))
        allocate(xiN(0:localnrows+1, 0:localncols+1))
        allocate(xiWbarx(localnrows, localncols+1))
        allocate(xiWbary(localnrows+1, localncols))
        allocate(xiNbarx(localnrows, localncols+1))
        allocate(xiNbary(localnrows+1, localncols))

        allocate(v(Nc, 0:localnrows+1, 0:localncols+1))

        allocate(Cf(0:localnrows+1, 0:localncols+1))

        allocate(viscW(0:localnrows+1, 0:localncols+1))
        allocate(viscN(0:localnrows+1, 0:localncols+1))

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
#endif

        if(pcol == 1) then
            Uwx(1:localnrows, 1) = UwBdryX(1, ylower:ylower+localnrows-1)
            Unx(1:localnrows, 1) = UnBdryX(1, ylower:ylower+localnrows-1)
        end if
        if(pcol == pncols) then
            Uwx(1:localnrows, localncols+1) = UwBdryX(2, ylower:ylower+localnrows-1)
            Unx(1:localnrows, localncols+1) = UnBdryX(2, ylower:ylower+localnrows-1)
        end if
        if(prow == 1) then
            Uwy(1, 1:localncols) = UwBdryY(xlower:xlower+localncols-1, 1)
            Uny(1, 1:localncols) = UnBdryY(xlower:xlower+localncols-1, 1)
        end if
        if(prow == pnrows) then
            Uwy(localnrows+1, 1:localncols) = UwBdryY(xlower:xlower+localncols-1, 2)
            Uny(localnrows+1, 1:localncols) = UnBdryY(xlower:xlower+localncols-1, 2)
        end if

        if(pncols == 1) then ! only one process column
            do i = 1, localnrows
                Kxxbar(i,1) = Kxx(1,ylower+i-1)
            end do
            do i = 1, localnrows
                Kxxbar(i,localncols+1) = Kxx(nx,ylower+i-1)
            end do
            do i = 1, localnrows
                do j = 2, localncols
                    Kxxbar(i,j) = (xs(xlower+j)-xs(xlower+j-2)) / ((xs(xlower+j-1)-xs(xlower+j-2))/&!
                        Kxx(xlower+j-2,ylower+i-1)+(xs(xlower+j)-xs(xlower+j-1))/Kxx(xlower+j-1,ylower+i-1))
                end do
            end do
        elseif(pcol == 1) then
            do i = 1, localnrows
                Kxxbar(i,1) = Kxx(1,ylower+i-1)
            end do
            do i = 1, localnrows
                do j = 2, localncols+1
                    Kxxbar(i,j) = (xs(xlower+j)-xs(xlower+j-2)) / ((xs(xlower+j-1)-xs(xlower+j-2))/&!
                        Kxx(xlower+j-2,ylower+i-1)+(xs(xlower+j)-xs(xlower+j-1))/Kxx(xlower+j-1,ylower+i-1))
                end do
            end do
        elseif(pcol == pncols) then
            do i = 1, localnrows
                Kxxbar(i,localncols+1) = Kxx(nx,ylower+i-1)
            end do
            do i = 1, localnrows
                do j = 1, localncols
                    Kxxbar(i,j) = (xs(xlower+j)-xs(xlower+j-2)) / ((xs(xlower+j-1)-xs(xlower+j-2))/&!
                        Kxx(xlower+j-2,ylower+i-1)+(xs(xlower+j)-xs(xlower+j-1))/Kxx(xlower+j-1,ylower+i-1))
                end do
            end do
        else
            do i = 1, localnrows
                do j = 1, localncols+1
                    Kxxbar(i,j) = (xs(xlower+j)-xs(xlower+j-2)) / ((xs(xlower+j-1)-xs(xlower+j-2))/&!
                        Kxx(xlower+j-2,ylower+i-1)+(xs(xlower+j)-xs(xlower+j-1))/Kxx(xlower+j-1,ylower+i-1))
                end do
            end do
        end if

        if(pnrows == 1) then ! only one process row
            do i = 1, localncols
                Kyybar(1,i) = Kyy(xlower+i-1,1)
            end do
            do i = 1, localncols
                Kyybar(localnrows+1,i) = Kyy(xlower+i-1, ny)
            end do
            do i = 2, localnrows
                do j = 1, localncols
                    Kyybar(i,j) = (ys(ylower+i)-ys(ylower+i-2)) / ((ys(ylower+i-1)-ys(ylower+i-2))/&!
                    Kyy(xlower+j-1,ylower+i-2)+(ys(ylower+i)-ys(ylower+i-1))/Kyy(xlower+j-1,ylower+i-1))
                end do
            end do
        elseif(prow == 1) then
            do i = 1, localncols
                Kyybar(1,i) = Kyy(xlower+i-1,1)
            end do
            do i = 2, localnrows+1
                do j = 1, localncols
                    Kyybar(i,j) = (ys(ylower+i)-ys(ylower+i-2)) / ((ys(ylower+i-1)-ys(ylower+i-2))/&!
                        Kyy(xlower+j-1,ylower+i-2)+(ys(ylower+i)-ys(ylower+i-1))/Kyy(xlower+j-1,ylower+i-1))
                end do
            end do
        elseif(prow == pnrows) then
            do i = 1, localncols
                Kyybar(localnrows+1,i) = Kyy(xlower+i-1, ny)
            end do
            do i = 1, localnrows
                do j = 1, localncols
                    Kyybar(i,j) = (ys(ylower+i)-ys(ylower+i-2)) / ((ys(ylower+i-1)-ys(ylower+i-2))/&!
                        Kyy(xlower+j-1,ylower+i-2)+(ys(ylower+i)-ys(ylower+i-1))/Kyy(xlower+j-1,ylower+i-1))
                end do
            end do
        else
            do i = 1, localnrows+1
                do j = 1, localncols
                    Kyybar(i,j) = (ys(ylower+i)-ys(ylower+i-2)) / ((ys(ylower+i-1)-ys(ylower+i-2))/&!
                        Kyy(xlower+j-1,ylower+i-2)+(ys(ylower+i)-ys(ylower+i-1))/Kyy(xlower+j-1,ylower+i-1))
                end do
            end do
        end if

        indexd = 0
        indexu = localnrows+1
        indexl = 0
        indexr = localncols+1
        if(prow == 1) then
            indexd = 1
        end if
        if(prow == pnrows) then
            indexu = localnrows
        end if
        if(pcol == 1) then
            indexl = 1
        end if
        if(pcol == pncols) then
            indexr = localncols
        end if

        do i = indexd, indexu
            do j = indexl, indexr
                Pw(i,j) = PwInit(xlower+j-1,ylower+i-1)
            end do
        end do
        if(prow==1) then
            do i = 1, localncols
                Pw(0, i) = PwBdryY(xlower+i-1, 1)
            end do
        end if
        if(prow==pnrows) then
            do i = 1, localncols
                Pw(localnrows+1, i) = PwBdryY(xlower+i-1, 2)
            end do
        end if
        if(pcol==1) then
            do i = 1, localnrows
                Pw(i, 0) = PwBdryX(1, ylower+i-1)
            end do
        end if
        if(pcol==pncols) then
            do i = 1, localnrows
                Pw(i, localncols+1) = PwBdryX(2, ylower+i-1)
            end do
        end if

        do m = 1, Nc
            do i = indexd, indexu
                do j = indexl, indexr
                    z(m,i,j) = zInit(m,xlower+j-1,ylower+i-1)
                end do
            end do
        end do

#ifdef SPARSE
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
            open(unit=140, file=FULLGRIDPREFIX//"xW3.txt", status='old')
            open(unit=170, file=FULLGRIDPREFIX//"xN3.txt", status='old')
            open(unit=250, file=FULLGRIDPREFIX//"v3.txt", status='old')
            read(140,*) xtable(3,1:TABLESIZE)
            read(170,*) ytable(3,1:TABLESIZE)
            read(250,*) vtable(3,1:TABLESIZE)
        end if

#endif

        totalmole = 0.0
        t = 2

        local_size = nx*ny/num_procs
        allocate(initial_x_guess(local_size))
        k = 1
        do i = 1, localnrows
            do j = 1, localncols
                initial_x_guess(k) = Pw(i,j)
                k = k + 1
            end do
        end do

        stencil_indices(1:5) = (/0, 1, 2, 3, 4/)
        offsets(1:5,1:2) = reshape([0,-1,1,0,0,0,0,0,-1,1], [5,2])

        ilower(1) = xlower
        iupper(1) = xupper
        ilower(2) = ylower
        iupper(2) = yupper

        call HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, grid, ierr)
        call HYPRE_StructGridSetExtents(grid, ilower, iupper, ierr)
        call HYPRE_StructGridAssemble(grid, ierr)

        call HYPRE_StructStencilCreate(2, 5, stencil, ierr)
        do k = 1, 5
            call HYPRE_StructStencilSetElement(stencil, k-1, offsets(k,1:2), ierr)
        end do

        call HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, global_A, ierr)
        call HYPRE_StructMatrixInitialize(global_A, ierr)

        call HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, global_b, ierr)
        call HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, global_x, ierr)
        call HYPRE_StructVectorInitialize(global_b, ierr)
        call HYPRE_StructVectorInitialize(global_x, ierr)

        call HYPRE_StructSMGCreate(MPI_COMM_WORLD, solver, ierr)
        call HYPRE_StructSMGSetMemoryUse(solver, 0, ierr)
        call HYPRE_StructSMGSetMaxIter(solver, 50, ierr)
        call HYPRE_StructSMGSetTol(solver, 1.0e-07, ierr)
        call HYPRE_StructSMGSetRelChange(solver, 0, ierr)
        call HYPRE_StructSMGSetNumPreRelax(solver, 1, ierr)
        call HYPRE_StructSMGSetNumPostRelax(solver, 1, ierr)
        call HYPRE_StructSMGSetLogging(solver, 1, ierr)

        if(myid == 0) then

            fmhtxt = trim(adjustl(soludoc))//"/soln_cmp2PhFlw_moleHistory.txt"
            fmrtxt = trim(adjustl(soludoc))//"/soln_cmp2PhFlw_moleRatioHistory.txt"
            allocate(fadtxt(Nc))
            do m = 1, Nc
                write(charm,'(i1)') m
                fadtxt(m) = trim(adjustl(soludoc))//"/soln_cmp2PhFlw_aveDiff"//charm//".txt"
            end do

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
            do m = 1, Nc
                open(unit=40+m, file=trim(adjustl(fadtxt(m))), status='replace', iostat=ierr)
                if(ierr /= 0) then
                    print *, 'open file error. ', ierr
                    call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
                end if
            end do
            deallocate(fadtxt)
        end if

        deallocate(modelCase%xs)
        deallocate(modelCase%ys)
        deallocate(modelCase%ts)
        deallocate(modelCase%Kxx)
        deallocate(modelCase%Kyy)
        deallocate(modelCase%poro)
        deallocate(modelCase%src)
        deallocate(modelCase%isDiriX)
        deallocate(modelCase%isDiriY)
        deallocate(modelCase%PwBdryX)
        deallocate(modelCase%PwBdryY)
        deallocate(modelCase%PwInit)
        deallocate(modelCase%zBdryX)
        deallocate(modelCase%zBdryY)
        deallocate(modelCase%zInit)
        deallocate(modelCase%UwBdryX)
        deallocate(modelCase%UwBdryY)
        deallocate(modelCase%UnBdryX)
        deallocate(modelCase%UnBdryY)
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
        real(kind=8), dimension(:,:,:), pointer :: moleincell
        real(kind=8), dimension(:), pointer :: leftmole, ad, adsum, sentbuffer, recvbuffer
        real(kind=8) :: totaldesiredleftmole, sum, leftmole1sum
        integer :: indexd, indexu, indexl, indexr
        integer :: i, j, m, k, prociteration

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

        indexd = 0
        indexu = localnrows+1
        indexl = 0
        indexr = localncols+1
        if(prow == 1) then
            indexd = 1
        end if
        if(prow == pnrows) then
            indexu = localnrows
        end if
        if(pcol == 1) then
            indexl = 1
        end if
        if(pcol == pncols) then
            indexr = localncols
        end if

        do j = indexl, indexr
            do i = indexd, indexu

                if((.not.((i==0).and.(j==0))).and.(.not.((i==0).and.(j==localncols+1))).and.(.not.((i==localnrows+1).and.(j==0))) &!
                    .and.(.not.((i==localnrows+1).and.(j==localncols+1)))) then

                    ztemp(1:Nc) = z(1:Nc,i,j)

!if((Pw(i,j)>3.300D6).or.(Pw(i,j)<1.800D6)) then
!print *, 'Out of the range of sparse grids.', pw(i,j),i,j
!call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
!end if

                    flashtimestart = MPI_Wtime()

#ifdef SPARSE
                    call flashcalculation_fullgrid(Pw(i,j), ztemp, xWtemp, xNtemp, xiW(i,j), xiN(i,j), &!
                        densiW(i,j), densiN(i,j), Sw(i,j), vtemp, Cf(i,j), isW, isN)
#else
                    call flashcalculation( Pw(i,j), ztemp, xWtemp, xNtemp, xiW(i,j), xiN(i,j), &!
                        densiW(i,j), densiN(i,j), Sw(i,j), vtemp, Cf(i,j), isW, isN, isRea )
#endif

                    flashtimefinish = MPI_Wtime()
                    flashtime = flashtime + flashtimefinish - flashtimestart

                    do m = 1, Nc
                        xW(m,i,j) = xWtemp(m)
                        xN(m,i,j) = xNtemp(m)
                        v(m,i,j) = vtemp(m)
                    end do

                    if(isW) then
                        viscW(i,j) = viscosity( xWtemp, xiW(i,j), Pw(i,j), 'l' )
                    else
                        viscW(i,j) = 1.D12
                    end if
                
                    if(isN) then
                        viscN(i,j) = viscosity( xNtemp, xiN(i,j), Pw(i,j), 'g' )
                    else
                        viscN(i,j) = 1.D12
                    end if
                end if

            end do
        end do

        allocate(moleincell(Nc, localnrows, localncols))
        do j = 1, localncols
            do i = 1, localnrows
                do m = 1, Nc
                    moleincell(m,i,j) = z(m,i,j)*(xiW(i,j)*Sw(i,j) + xiN(i,j)*(1-Sw(i,j)))&!
                        *(xs(xlower+j)-xs(xlower+j-1))*(ys(ylower+i)-ys(ylower+i-1))*poro(xlower+j-1,ylower+i-1)
                end do
            end do
        end do

        allocate(leftmole(Nc))
        leftmole = 0
        do m = 1, Nc
            do j = 1, localncols
                do i = 1, localnrows
                    leftmole(m) = leftmole(m) + moleincell(m,i,j)
                end do
            end do
        end do

        totaldesiredleftmole = 0.0
        do m = 2, Nc
            totaldesiredleftmole = totaldesiredleftmole + leftmole(m)
        end do

        allocate(ad(Nc))
        ad = 0
        do m = 1, Nc
            do j = 1, localncols
                do i = 1, localnrows
                    ad(m) = ad(m) + abs((1/xiN(i,j)-v(m,i,j))/(1/xiN(i,j)))
                end do
            end do
        end do

        do m = 1, Nc
            ad(m) = ad(m) / (nx*ny/pncols/pnrows)
        end do

        call MPI_BUFFER_ATTACH(buffer,buffer_size,ierr)

        if(myid /= 0) then

            position = 0
            allocate(sentbuffer(2+Nc))
            call MPI_PACK(totaldesiredleftmole, 1, MPI_DOUBLE_PRECISION, sentbuffer, (2+Nc)*8, position, MPI_COMM_WORLD, ierr)
            call MPI_PACK(leftmole(1), 1, MPI_DOUBLE_PRECISION, sentbuffer, (2+Nc)*8, position, MPI_COMM_WORLD, ierr)
            do m = 1, Nc
                call MPI_PACK(ad(m), 1, MPI_DOUBLE_PRECISION, sentbuffer, (2+Nc)*8, position, MPI_COMM_WORLD, ierr)
            end do
            call MPI_IBSEND(sentbuffer,2+Nc,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,request,ierr)
            call MPI_WAIT(request, status, ierr)
            deallocate(sentbuffer)
            
        else

            sum = totaldesiredleftmole
            leftmole1sum = leftmole(1)
            allocate(adsum(Nc))
            do m = 1, Nc
                adsum(m) = ad(m)
            end do

            allocate(recvbuffer(2+Nc))
            do prociteration = 1, num_procs-1

                call MPI_RECV(recvbuffer, 2+Nc, MPI_DOUBLE_PRECISION, prociteration, prociteration, &!
                    MPI_COMM_WORLD, status, ierr)

                sum = sum + recvbuffer(1)
                leftmole1sum = leftmole1sum + recvbuffer(2)
                do m = 1, Nc
                    adsum(m) = adsum(m) + recvbuffer(2+m)
                end do

            end do
            deallocate(recvbuffer)

            if(t == 2) then
                totalmole = sum
            end if

            write(20, fmt="(es12.5)") (totalmole-sum)/totalmole
            write(30, fmt="(es12.5)") sum/leftmole1sum
            do m = 1, Nc
                write(40+m, fmt="(es12.5)") adsum(m)/num_procs
            end do

            deallocate(adsum)

        end if

        call MPI_BUFFER_DETACH(buffer,buffer_size,ierr)

        deallocate(moleincell)
        deallocate(leftmole)
        deallocate(ad)

        if(pncols == 1) then
            do i = 1, localnrows
                if((Uwx(i,1) > 0).or.(Unx(i,1) > 0)) then
                    ztemp(1:Nc) = zBdryX(1:Nc,1,ylower+i-1)
#ifdef SPARSE
                    call flashcalculation_fullgrid( Pw(i,1), ztemp, xWtemp, xNtemp, xiWbarx(i,1), xiNbarx(i,1), &!
                        densiWbarx(i,1), densiNbarx(i,1), Swtemp, vtemp, ignore1, isW, isN )
#else
                    call flashcalculation(Pw(i,1), ztemp, xWtemp, xNtemp, xiWbarx(i,1), xiNbarx(i,1), &!
                        densiWbarx(i,1), densiNbarx(i,1), Swtemp, vtemp, ignore1, isW, isN, isRea)
#endif
                    xWbarx(1:Nc,i,1) = xWtemp(1:Nc)
                    xNbarx(1:Nc,i,1) = xNtemp(1:Nc)
                    if(isW) then
                        viscWtemp = viscosity( xWtemp, xiWbarx(i,1), Pw(i,1), 'l' )
                    end if
                    if(isN) then
                        viscNtemp = viscosity( xNtemp, xiNbarx(i,1), Pw(i,1), 'g' )
                    end if
                end if
                if(Uwx(i,1) > 0) then
                    lambdawx(i,1) = Kxxbar(i,1)*computekr_W(Swtemp)/viscWtemp
                elseif(Uwx(i,1) < 0) then
                    lambdawx(i,1) = Kxxbar(i,1)*computekr_W(Sw(i,1))/viscW(i,1)
                    xiWbarx(i,1) = xiW(i,1)
                    densiWbarx(i,1) = densiW(i,1)
                    xWbarx(1:Nc,i,1) = xW(1:Nc,i,1)
                elseif(isDiriX(1,i) == 1) then
                    lambdawx(i,1) = Kxxbar(i,1)*computekr_W(Sw(i,1))/viscW(i,1)
                    xiWbarx(i,1) = xiW(i,1)
                    densiWbarx(i,1) = densiW(i,1)
                    xWbarx(1:Nc,i,1) = xW(1:Nc,i,1)
                end if
                if(Unx(i,1) > 0) then
                    lambdanx(i,1) = Kxxbar(i,1)*computekr_N(Swtemp)/viscNtemp
                elseif(Unx(i,1) < 0) then
                    lambdanx(i,1) = Kxxbar(i,1)*computekr_N(Sw(i,1))/viscN(i,1)
                    xiNbarx(i,1) = xiN(i,1)
                    densiNbarx(i,1) = densiN(i,1)
                    xNbarx(1:Nc,i,1) = xN(1:Nc,i,1)
                elseif(isDiriX(1,i) == 1) then
                    lambdanx(i,1) = Kxxbar(i,1)*computekr_N(Sw(i,1))/viscN(i,1)
                    xiNbarx(i,1) = xiN(i,1)
                    densiNbarx(i,1) = densiN(i,1)
                    xNbarx(1:Nc,i,1) = xN(1:Nc,i,1)
                end if
            end do
            do i = 1, localnrows
                if((Uwx(i,localncols+1) < 0).or.(Unx(i,localncols+1) < 0)) then
                    ztemp(1:Nc) = zBdryX(1:Nc,2,ylower+i-1)
#ifdef SPARSE
                    call flashcalculation_fullgrid( Pw(i,localncols), ztemp, xWtemp, xNtemp, xiWbarx(i,localncols+1), &!
                        xiNbarx(i,localncols+1), densiWbarx(i,localncols+1), densiNbarx(i,localncols+1), Swtemp, vtemp, &!
                        ignore1, isW, isN )
#else
                    call flashcalculation(Pw(i,localncols), ztemp, xWtemp, &!
                        xNtemp, xiWbarx(i,localncols+1), xiNbarx(i,localncols+1), densiWbarx(i,localncols+1), &!
                        densiNbarx(i,localncols+1), Swtemp, vtemp, ignore1, isW, isN, isRea)
#endif
                    xWbarx(1:Nc,i,localncols+1) = xWtemp(1:Nc)
                    xNbarx(1:Nc,i,localncols+1) = xNtemp(1:Nc)
                    if(isW) then
                        viscWtemp = viscosity( xWtemp, xiWbarx(i,localncols+1), Pw(i,localncols), 'l' )
                    end if
                    if(isN) then
                        viscNtemp = viscosity( xNtemp, xiNbarx(i,localncols+1), Pw(i,localncols), 'g' )
                    end if
                end if
                if(Uwx(i,localncols+1) < 0) then
                    lambdawx(i,localncols+1) = Kxxbar(i,localncols+1)*computekr_W(Swtemp)/viscWtemp
                elseif(Uwx(i,localncols+1) > 0) then
                    lambdawx(i,localncols+1) = Kxxbar(i,localncols+1)*computekr_W(Sw(i,localncols))/viscW(i,localncols)
                    xiWbarx(i,localncols+1) = xiW(i,localncols)
                    densiWbarx(i,localncols+1) = densiW(i,localncols)
                    xWbarx(1:Nc,i,localncols+1) = xW(1:Nc,i,localncols)
                elseif(isDiriX(2,ylower+i-1) == 1) then
                    lambdawx(i,localncols+1) = Kxxbar(i,localncols+1)*computekr_W(Sw(i,localncols))/viscW(i,localncols)
                    xiWbarx(i,localncols+1) = xiW(i,localncols)
                    densiWbarx(i,localncols+1) = densiW(i,localncols)
                    xWbarx(1:Nc,i,localncols+1) = xW(1:Nc,i,localncols)
                end if
                if(Unx(i,localncols+1) < 0) then
                    lambdanx(i,localncols+1) = Kxxbar(i,localncols+1)*computekr_N(Swtemp)/viscNtemp
                elseif(Unx(i,localncols+1) > 0) then
                    lambdanx(i,localncols+1) = Kxxbar(i,localncols+1)*computekr_N(Sw(i,localncols))/viscN(i,localncols)
                    xiNbarx(i,localncols+1) = xiN(i,localncols)
                    densiNbarx(i,localncols+1) = densiN(i,localncols)
                    xNbarx(1:Nc,i,localncols+1) = xN(1:Nc,i,localncols)
                elseif(isDiriX(2,ylower+i-1) == 1) then
                    lambdanx(i,localncols+1) = Kxxbar(i,localncols+1)*computekr_N(Sw(i,localncols))/viscN(i,localncols)
                    xiNbarx(i,localncols+1) = xiN(i,localncols)
                    densiNbarx(i,localncols+1) = densiN(i,localncols)
                    xNbarx(1:Nc,i,localncols+1) = xN(1:Nc,i,localncols)
                end if
            end do
            do j = 2, localncols
                do i = 1, localnrows
                    if(Uwx(i,j) > 0) then
                        lambdawx(i,j) = Kxxbar(i,j)*computekr_W(Sw(i,j-1))/viscW(i,j-1)
                        xiWbarx(i,j) = xiW(i,j-1)
                        densiWbarx(i,j) = densiW(i,j-1)
                        xWbarx(1:Nc,i,j) = xW(1:Nc,i,j-1)
                    else
                        lambdawx(i,j) = Kxxbar(i,j)*computekr_W(Sw(i,j))/viscW(i,j)
                        xiWbarx(i,j) = xiW(i,j)
                        densiWbarx(i,j) = densiW(i,j)
                        xWbarx(1:Nc,i,j) = xW(1:Nc,i,j)
                    end if
                    if(Unx(i,j) > 0) then
                        lambdanx(i,j) = Kxxbar(i,j)*computekr_N(Sw(i,j-1))/viscN(i,j-1)
                        xiNbarx(i,j) = xiN(i,j-1)
                        densiNbarx(i,j) = densiN(i,j-1)
                        xNbarx(1:Nc,i,j) = xN(1:Nc,i,j-1)
                    else
                        lambdanx(i,j) = Kxxbar(i,j)*computekr_N(Sw(i,j))/viscN(i,j)
                        xiNbarx(i,j) = xiN(i,j)
                        densiNbarx(i,j) = densiN(i,j)
                        xNbarx(1:Nc,i,j) = xN(1:Nc,i,j)
                    end if
                end do
            end do
        elseif(pcol==1) then
            do i = 1, localnrows
                if((Uwx(i,1) > 0).or.(Unx(i,1) > 0)) then
                    ztemp(1:Nc) = zBdryX(1:Nc,1,ylower+i-1)
#ifdef SPARSE
                    call flashcalculation_fullgrid( Pw(i,1), ztemp, xWtemp, xNtemp, xiWbarx(i,1), xiNbarx(i,1), &!
                        densiWbarx(i,1), densiNbarx(i,1), Swtemp, vtemp, ignore1, isW, isN )
#else
                    call flashcalculation(Pw(i,1), ztemp, xWtemp, xNtemp, xiWbarx(i,1), xiNbarx(i,1), &!
                        densiWbarx(i,1), densiNbarx(i,1), Swtemp, vtemp, ignore1, isW, isN, isRea)
#endif
                    xWbarx(1:Nc,i,1) = xWtemp(1:Nc)
                    xNbarx(1:Nc,i,1) = xNtemp(1:Nc)
                    if(isW) then
                        viscWtemp = viscosity( xWtemp, xiWbarx(i,1), Pw(i,1), 'l' )
                    end if
                    if(isN) then
                        viscNtemp = viscosity( xNtemp, xiNbarx(i,1), Pw(i,1), 'g' )
                    end if
                end if
                if(Uwx(i,1) > 0) then
                    lambdawx(i,1) = Kxxbar(i,1)*computekr_W(Swtemp)/viscWtemp
                elseif(Uwx(i,1) < 0) then
                    lambdawx(i,1) = Kxxbar(i,1)*computekr_W(Sw(i,1))/viscW(i,1)
                    xiWbarx(i,1) = xiW(i,1)
                    densiWbarx(i,1) = densiW(i,1)
                    xWbarx(1:Nc,i,1) = xW(1:Nc,i,1)
                elseif(isDiriX(1,ylower+i-1) == 1) then
                    lambdawx(i,1) = Kxxbar(i,1)*computekr_W(Sw(i,1))/viscW(i,1)
                    xiWbarx(i,1) = xiW(i,1)
                    densiWbarx(i,1) = densiW(i,1)
                    xWbarx(1:Nc,i,1) = xW(1:Nc,i,1)
                end if
                if(Unx(i,1) > 0) then
                    lambdanx(i,1) = Kxxbar(i,1)*computekr_N(Swtemp)/viscNtemp
                elseif(Unx(i,1) < 0) then
                    lambdanx(i,1) = Kxxbar(i,1)*computekr_N(Sw(i,1))/viscN(i,1)
                    xiNbarx(i,1) = xiN(i,1)
                    densiNbarx(i,1) = densiN(i,1)
                    xNbarx(1:Nc,i,1) = xN(1:Nc,i,1)
                elseif(isDiriX(1,ylower+i-1) == 1) then
                    lambdanx(i,1) = Kxxbar(i,1)*computekr_N(Sw(i,1))/viscN(i,1)
                    xiNbarx(i,1) = xiN(i,1)
                    densiNbarx(i,1) = densiN(i,1)
                    xNbarx(1:Nc,i,1) = xN(1:Nc,i,1)
                end if
            end do
            do j = 2, localncols+1
                do i = 1, localnrows
                    if(Uwx(i,j) > 0) then
                        lambdawx(i,j) = Kxxbar(i,j)*computekr_W(Sw(i,j-1))/viscW(i,j-1)
                        xiWbarx(i,j) = xiW(i,j-1)
                        densiWbarx(i,j) = densiW(i,j-1)
                        xWbarx(1:Nc,i,j) = xW(1:Nc,i,j-1)
                    else
                        lambdawx(i,j) = Kxxbar(i,j)*computekr_W(Sw(i,j))/viscW(i,j)
                        xiWbarx(i,j) = xiW(i,j)
                        densiWbarx(i,j) = densiW(i,j)
                        xWbarx(1:Nc,i,j) = xW(1:Nc,i,j)
                    end if
                    if(Unx(i,j) > 0) then
                        lambdanx(i,j) = Kxxbar(i,j)*computekr_N(Sw(i,j-1))/viscN(i,j-1)
                        xiNbarx(i,j) = xiN(i,j-1)
                        densiNbarx(i,j) = densiN(i,j-1)
                        xNbarx(1:Nc,i,j) = xN(1:Nc,i,j-1)
                    else
                        lambdanx(i,j) = Kxxbar(i,j)*computekr_N(Sw(i,j))/viscN(i,j)
                        xiNbarx(i,j) = xiN(i,j)
                        densiNbarx(i,j) = densiN(i,j)
                        xNbarx(1:Nc,i,j) = xN(1:Nc,i,j)
                    end if
                end do
            end do
        elseif(pcol==pncols) then
            do i = 1, localnrows
                if((Uwx(i,localncols+1) < 0).or.(Unx(i,localncols+1) < 0)) then
                    ztemp(1:Nc) = zBdryX(1:Nc,2,ylower+i-1)
#ifdef SPARSE
                    call flashcalculation_fullgrid( Pw(i,localncols), ztemp, xWtemp, &!
                        xNtemp, xiWbarx(i,localncols+1), xiNbarx(i,localncols+1), densiWbarx(i,localncols+1), &!
                        densiNbarx(i,localncols+1), Swtemp, vtemp, ignore1, isW, isN )
#else
                    call flashcalculation(Pw(i,localncols), ztemp, xWtemp, &!
                        xNtemp, xiWbarx(i,localncols+1), xiNbarx(i,localncols+1), densiWbarx(i,localncols+1), &!
                        densiNbarx(i,localncols+1), Swtemp, vtemp, ignore1, isW, isN, isRea)
#endif
                    xWbarx(1:Nc,i,localncols+1) = xWtemp(1:Nc)
                    xNbarx(1:Nc,i,localncols+1) = xNtemp(1:Nc)
                    if(isW) then
                        viscWtemp = viscosity( xWtemp, xiWbarx(i,localncols+1), Pw(i,localncols), 'l' )
                    end if
                    if(isN) then
                        viscNtemp = viscosity( xNtemp, xiNbarx(i,localncols+1), Pw(i,localncols), 'g' )
                    end if
                end if
                if(Uwx(i,localncols+1) < 0) then
                    lambdawx(i,localncols+1) = Kxxbar(i,localncols+1)*computekr_W(Swtemp)/viscWtemp
                elseif(Uwx(i,localncols+1) > 0) then
                    lambdawx(i,localncols+1) = Kxxbar(i,localncols+1)*computekr_W(Sw(i,localncols))/viscW(i,localncols)
                    xiWbarx(i,localncols+1) = xiW(i,localncols)
                    densiWbarx(i,localncols+1) = densiW(i,localncols)
                    xWbarx(1:Nc,i,localncols+1) = xW(1:Nc,i,localncols)
                elseif(isDiriX(2,ylower+i-1) == 1) then
                    lambdawx(i,localncols+1) = Kxxbar(i,localncols+1)*computekr_W(Sw(i,localncols))/viscW(i,localncols)
                    xiWbarx(i,localncols+1) = xiW(i,localncols)
                    densiWbarx(i,localncols+1) = densiW(i,localncols)
                    xWbarx(1:Nc,i,localncols+1) = xW(1:Nc,i,localncols)
                end if
                if(Unx(i,localncols+1) < 0) then
                    lambdanx(i,localncols+1) = Kxxbar(i,localncols+1)*computekr_N(Swtemp)/viscNtemp
                elseif(Unx(i,localncols+1) > 0) then
                    lambdanx(i,localncols+1) = Kxxbar(i,localncols+1)*computekr_N(Sw(i,localncols))/viscN(i,localncols)
                    xiNbarx(i,localncols+1) = xiN(i,localncols)
                    densiNbarx(i,localncols+1) = densiN(i,localncols)
                    xNbarx(1:Nc,i,localncols+1) = xN(1:Nc,i,localncols)
                elseif(isDiriX(2,ylower+i-1) == 1) then
                    lambdanx(i,localncols+1) = Kxxbar(i,localncols+1)*computekr_N(Sw(i,localncols))/viscN(i,localncols)
                    xiNbarx(i,localncols+1) = xiN(i,localncols)
                    densiNbarx(i,localncols+1) = densiN(i,localncols)
                    xNbarx(1:Nc,i,localncols+1) = xN(1:Nc,i,localncols)
                end if
            end do
            do j = 1, localncols
                do i = 1, localnrows
                    if(Uwx(i,j) > 0) then
                        lambdawx(i,j) = Kxxbar(i,j)*computekr_W(Sw(i,j-1))/viscW(i,j-1)
                        xiWbarx(i,j) = xiW(i,j-1)
                        densiWbarx(i,j) = densiW(i,j-1)
                        xWbarx(1:Nc,i,j) = xW(1:Nc,i,j-1)
                    else
                        lambdawx(i,j) = Kxxbar(i,j)*computekr_W(Sw(i,j))/viscW(i,j)
                        xiWbarx(i,j) = xiW(i,j)
                        densiWbarx(i,j) = densiW(i,j)
                        xWbarx(1:Nc,i,j) = xW(1:Nc,i,j)
                    end if
                    if(Unx(i,j) > 0) then
                        lambdanx(i,j) = Kxxbar(i,j)*computekr_N(Sw(i,j-1))/viscN(i,j-1)
                        xiNbarx(i,j) = xiN(i,j-1)
                        densiNbarx(i,j) = densiN(i,j-1)
                        xNbarx(1:Nc,i,j) = xN(1:Nc,i,j-1)
                    else
                        lambdanx(i,j) = Kxxbar(i,j)*computekr_N(Sw(i,j))/viscN(i,j)
                        xiNbarx(i,j) = xiN(i,j)
                        densiNbarx(i,j) = densiN(i,j)
                        xNbarx(1:Nc,i,j) = xN(1:Nc,i,j)
                    end if
                end do
            end do
        else
            do j = 1, localncols+1
                do i = 1, localnrows
                    if(Uwx(i,j) > 0) then
                        lambdawx(i,j) = Kxxbar(i,j)*computekr_W(Sw(i,j-1))/viscW(i,j-1)
                        xiWbarx(i,j) = xiW(i,j-1)
                        densiWbarx(i,j) = densiW(i,j-1)
                        xWbarx(1:Nc,i,j) = xW(1:Nc,i,j-1)
                    else
                        lambdawx(i,j) = Kxxbar(i,j)*computekr_W(Sw(i,j))/viscW(i,j)
                        xiWbarx(i,j) = xiW(i,j)
                        densiWbarx(i,j) = densiW(i,j)
                        xWbarx(1:Nc,i,j) = xW(1:Nc,i,j)
                    end if
                    if(Unx(i,j) > 0) then
                        lambdanx(i,j) = Kxxbar(i,j)*computekr_N(Sw(i,j-1))/viscN(i,j-1)
                        xiNbarx(i,j) = xiN(i,j-1)
                        densiNbarx(i,j) = densiN(i,j-1)
                        xNbarx(1:Nc,i,j) = xN(1:Nc,i,j-1)
                    else
                        lambdanx(i,j) = Kxxbar(i,j)*computekr_N(Sw(i,j))/viscN(i,j)
                        xiNbarx(i,j) = xiN(i,j)
                        densiNbarx(i,j) = densiN(i,j)
                        xNbarx(1:Nc,i,j) = xN(1:Nc,i,j)
                    end if
                end do
            end do
        end if

        if(pnrows == 1) then
            do i = 1, localncols
                if((Uwy(1,i) > 0).or.(Uny(1,i) > 0)) then
                    ztemp(1:Nc) = zBdryY(1:Nc,xlower+i-1,1)
#ifdef SPARSE
                    call flashcalculation_fullgrid( Pw(1,i), ztemp, xWtemp, xNtemp, xiWbary(1,i), xiNbary(1,i), &!
                        densiWbary(1,i), densiNbary(1,i), Swtemp, vtemp, ignore1, isW, isN )
#else
                    call flashcalculation(Pw(1,i), ztemp, xWtemp, xNtemp, xiWbary(1,i), xiNbary(1,i), &!
                        densiWbary(1,i), densiNbary(1,i), Swtemp, vtemp, ignore1, isW, isN, isRea)
#endif
                    xWbary(1:Nc,1,i) = xWtemp(1:Nc)
                    xNbary(1:Nc,1,i) = xNtemp(1:Nc)
                    if(isW) then
                        viscWtemp = viscosity( xWtemp, xiWbary(1,i), Pw(1,i), 'l' )
                    end if
                    if(isN) then
                        viscNtemp = viscosity( xNtemp, xiNbary(1,i), Pw(1,i), 'g' )
                    end if
                end if
                if(Uwy(1,i) > 0) then
                    lambdawy(1,i) = Kyybar(1,i)*computekr_W(Swtemp)/viscWtemp
                elseif(Uwy(1,i) < 0) then
                    lambdawy(1,i) = Kyybar(1,i)*computekr_W(Sw(1,i))/viscW(1,i)
                    xiWbary(1,i) = xiW(1,i)
                    densiWbary(1,i) = densiW(1,i)
                    xWbary(1:Nc,1,i) = xW(1:Nc,1,i)
                elseif(isDiriY(xlower+i-1,1) == 1) then
                    lambdawy(1,i) = Kyybar(1,i)*computekr_W(Sw(1,i))/viscW(1,i)
                    xiWbary(1,i) = xiW(1,i)
                    densiWbary(1,i) = densiW(1,i)
                    xWbary(1:Nc,1,i) = xW(1:Nc,1,i)
                end if
                if(Uny(1,i) > 0) then
                    lambdany(1,i) = Kyybar(1,i)*computekr_N(Swtemp)/viscNtemp
                elseif(Uny(1,i) < 0) then
                    lambdany(1,i) = Kyybar(1,i)*computekr_N(Sw(1,i))/viscN(1,i)
                    xiNbary(1,i) = xiN(1,i)
                    densiNbary(1,i) = densiN(1,i)
                    xNbary(1:Nc,1,i) = xN(1:Nc,1,i)
                elseif(isDiriY(xlower+i-1,1) == 1) then
                    lambdany(1,i) = Kyybar(1,i)*computekr_N(Sw(1,i))/viscN(1,i)
                    xiNbary(1,i) = xiN(1,i)
                    densiNbary(1,i) = densiN(1,i)
                    xNbary(1:Nc,1,i) = xN(1:Nc,1,i)
                end if
            end do
            do i = 1, localncols
                if((Uwy(localnrows+1,i) < 0).or.(Uny(localnrows+1,i) < 0)) then
                    ztemp(1:Nc) = zBdryY(1:Nc,xlower+i-1,2)
#ifdef SPARSE
                    call flashcalculation_fullgrid( Pw(localnrows,i), ztemp, xWtemp, xNtemp, xiWbary(localnrows+1,i), &!
                        xiNbary(localnrows+1,i), densiWbary(localnrows+1,i), densiNbary(localnrows+1,i), Swtemp, vtemp, &!
                        ignore1, isW, isN )
#else
                    call flashcalculation(Pw(localnrows,i), ztemp, xWtemp, xNtemp, xiWbary(localnrows+1,i), &!
                        xiNbary(localnrows+1,i), densiWbary(localnrows+1,i), densiNbary(localnrows+1,i), Swtemp, vtemp, &!
                        ignore1, isW, isN, isRea)
#endif
                    xWbary(1:Nc,localnrows+1,i) = xWtemp(1:Nc)
                    xNbary(1:Nc,localnrows+1,i) = xNtemp(1:Nc)
                    if(isW) then
                        viscWtemp = viscosity( xWtemp, xiWbary(localnrows+1,i), Pw(localnrows,i), 'l' )
                    end if
                    if(isN) then
                        viscNtemp = viscosity( xNtemp, xiNbary(localnrows+1,i), Pw(localnrows,i), 'g' )
                    end if
                end if
                if(Uwy(localnrows+1,i) < 0) then
                    lambdawy(localnrows+1,i) = Kyybar(localnrows+1,i)*computekr_W(Swtemp)/viscWtemp
                elseif(Uwy(localnrows+1,i) > 0) then
                    lambdawy(localnrows+1,i) = Kyybar(localnrows+1,i)*computekr_W(Sw(localnrows,i))/viscW(localnrows,i)
                    xiWbary(localnrows+1,i) = xiW(localnrows,i)
                    densiWbary(localnrows+1,i) = densiW(localnrows,i)
                    xWbary(1:Nc,localnrows+1,i) = xW(1:Nc,localnrows,i)
                elseif(isDiriY(xlower+i-1,2) == 1) then
                    lambdawy(localnrows+1,i) = Kyybar(localnrows+1,i)*computekr_W(Sw(localnrows,i))/viscW(localnrows,i)
                    xiWbary(localnrows+1,i) = xiW(localnrows,i)
                    densiWbary(localnrows+1,i) = densiW(localnrows,i)
                    xWbary(1:Nc,localnrows+1,i) = xW(1:Nc,localnrows,i)
                end if
                if(Uny(localnrows+1,i) < 0) then
                    lambdany(localnrows+1,i) = Kyybar(localnrows+1,i)*computekr_N(Swtemp)/viscNtemp
                elseif(Uny(localnrows+1,i) > 0) then
                    lambdany(localnrows+1,i) = Kyybar(localnrows+1,i)*computekr_N(Sw(localnrows,i))/viscN(localnrows,i)
                    xiNbary(localnrows+1,i) = xiN(localnrows,i)
                    densiNbary(localnrows+1,i) = densiN(localnrows,i)
                    xNbary(1:Nc,localnrows+1,i) = xN(1:Nc,localnrows,i)
                elseif(isDiriY(xlower+i-1,2) == 1) then
                    lambdany(localnrows+1,i) = Kyybar(localnrows+1,i)*computekr_N(Sw(localnrows,i))/viscN(localnrows,i)
                    xiNbary(localnrows+1,i) = xiN(localnrows,i)
                    densiNbary(localnrows+1,i) = densiN(localnrows,i)
                    xNbary(1:Nc,localnrows+1,i) = xN(1:Nc,localnrows,i)
                end if
            end do
            do j = 1, localncols
                do i = 2, localnrows
                    if(Uwy(i,j) > 0) then
                        lambdawy(i,j) = Kyybar(i,j)*computekr_W(Sw(i-1,j))/viscW(i-1,j)
                        xiWbary(i,j) = xiW(i-1,j)
                        densiWbary(i,j) = densiW(i-1,j)
                        xWbary(1:Nc,i,j) = xW(1:Nc,i-1,j)
                    else
                        lambdawy(i,j) = Kyybar(i,j)*computekr_W(Sw(i,j))/viscW(i,j)
                        xiWbary(i,j) = xiW(i,j)
                        densiWbary(i,j) = densiW(i,j)
                        xWbary(1:Nc,i,j) = xW(1:Nc,i,j)
                    end if
                    if(Uny(i,j) > 0) then
                        lambdany(i,j) = Kyybar(i,j)*computekr_N(Sw(i-1,j))/viscN(i-1,j)
                        xiNbary(i,j) = xiN(i-1,j)
                        densiNbary(i,j) = densiN(i-1,j)
                        xNbary(1:Nc,i,j) = xN(1:Nc,i-1,j)
                    else
                        lambdany(i,j) = Kyybar(i,j)*computekr_N(Sw(i,j))/viscN(i,j)
                        xiNbary(i,j) = xiN(i,j)
                        densiNbary(i,j) = densiN(i,j)
                        xNbary(1:Nc,i,j) = xN(1:Nc,i,j)
                    end if
                end do
            end do
        elseif(prow==1) then
            do i = 1, localncols
                if((Uwy(1,i) > 0).or.(Uny(1,i) > 0)) then
                    ztemp(1:Nc) = zBdryY(1:Nc,xlower+i-1,1)
#ifdef SPARSE
                    call flashcalculation_fullgrid( Pw(1,i), ztemp, xWtemp, xNtemp, xiWbary(1,i), xiNbary(1,i), &!
                        densiWbary(1,i), densiNbary(1,i), Swtemp, vtemp, ignore1, isW, isN )
#else
                    call flashcalculation(Pw(1,i), ztemp, xWtemp, xNtemp, xiWbary(1,i), xiNbary(1,i), &!
                        densiWbary(1,i), densiNbary(1,i), Swtemp, vtemp, ignore1, isW, isN, isRea)
#endif
                    xWbary(1:Nc,1,i) = xWtemp(1:Nc)
                    xNbary(1:Nc,1,i) = xNtemp(1:Nc)
                    if(isW) then
                        viscWtemp = viscosity( xWtemp, xiWbary(1,i), Pw(1,i), 'l' )
                    end if
                    if(isN) then
                        viscNtemp = viscosity( xNtemp, xiNbary(1,i), Pw(1,i), 'g' )
                    end if
                end if
                if(Uwy(1,i) > 0) then
                    lambdawy(1,i) = Kyybar(1,i)*computekr_W(Swtemp)/viscWtemp
                elseif(Uwy(1,i) < 0) then
                    lambdawy(1,i) = Kyybar(1,i)*computekr_W(Sw(1,i))/viscW(1,i)
                    xiWbary(1,i) = xiW(1,i)
                    densiWbary(1,i) = densiW(1,i)
                    xWbary(1:Nc,1,i) = xW(1:Nc,1,i)
                elseif(isDiriY(xlower+i-1,1) == 1) then
                    lambdawy(1,i) = Kyybar(1,i)*computekr_W(Sw(1,i))/viscW(1,i)
                    xiWbary(1,i) = xiW(1,i)
                    densiWbary(1,i) = densiW(1,i)
                    xWbary(1:Nc,1,i) = xW(1:Nc,1,i)
                end if
                if(Uny(1,i) > 0) then
                    lambdany(1,i) = Kyybar(1,i)*computekr_N(Swtemp)/viscNtemp
                elseif(Uny(1,i) < 0) then
                    lambdany(1,i) = Kyybar(1,i)*computekr_N(Sw(1,i))/viscN(1,i)
                    xiNbary(1,i) = xiN(1,i)
                    densiNbary(1,i) = densiN(1,i)
                    xNbary(1:Nc,1,i) = xN(1:Nc,1,i)
                elseif(isDiriY(xlower+i-1,1) == 1) then
                    lambdany(1,i) = Kyybar(1,i)*computekr_N(Sw(1,i))/viscN(1,i)
                    xiNbary(1,i) = xiN(1,i)
                    densiNbary(1,i) = densiN(1,i)
                    xNbary(1:Nc,1,i) = xN(1:Nc,1,i)
                end if
            end do
            do j = 1, localncols
                do i = 2, localnrows+1
                    if(Uwy(i,j) > 0) then
                        lambdawy(i,j) = Kyybar(i,j)*computekr_W(Sw(i-1,j))/viscW(i-1,j)
                        xiWbary(i,j) = xiW(i-1,j)
                        densiWbary(i,j) = densiW(i-1,j)
                        xWbary(1:Nc,i,j) = xW(1:Nc,i-1,j)
                    else
                        lambdawy(i,j) = Kyybar(i,j)*computekr_W(Sw(i,j))/viscW(i,j)
                        xiWbary(i,j) = xiW(i,j)
                        densiWbary(i,j) = densiW(i,j)
                        xWbary(1:Nc,i,j) = xW(1:Nc,i,j)
                    end if
                    if(Uny(i,j) > 0) then
                        lambdany(i,j) = Kyybar(i,j)*computekr_N(Sw(i-1,j))/viscN(i-1,j)
                        xiNbary(i,j) = xiN(i-1,j)
                        densiNbary(i,j) = densiN(i-1,j)
                        xNbary(1:Nc,i,j) = xN(1:Nc,i-1,j)
                    else
                        lambdany(i,j) = Kyybar(i,j)*computekr_N(Sw(i,j))/viscN(i,j)
                        xiNbary(i,j) = xiN(i,j)
                        densiNbary(i,j) = densiN(i,j)
                        xNbary(1:Nc,i,j) = xN(1:Nc,i,j)
                    end if
                end do
            end do
        elseif(prow==pnrows) then
            do i = 1, localncols
                if((Uwy(localnrows+1,i) < 0).or.(Uny(localnrows+1,i) < 0)) then
                    ztemp(1:Nc) = zBdryY(1:Nc,xlower+i-1,2)
#ifdef SPARSE
                    call flashcalculation_fullgrid( Pw(localnrows,i), ztemp, xWtemp, &!
                        xNtemp, xiWbary(localnrows+1,i), xiNbary(localnrows+1,i), densiWbary(localnrows+1,i), &!
                        densiNbary(localnrows+1,i), Swtemp, vtemp, ignore1, isW, isN )
#else
                    call flashcalculation(Pw(localnrows,i), ztemp, xWtemp, &!
                        xNtemp, xiWbary(localnrows+1,i), xiNbary(localnrows+1,i), densiWbary(localnrows+1,i), &!
                        densiNbary(localnrows+1,i), Swtemp, vtemp, ignore1, isW, isN, isRea)
#endif
                    xWbary(1:Nc,localnrows+1,i) = xWtemp(1:Nc)
                    xNbary(1:Nc,localnrows+1,i) = xNtemp(1:Nc)
                    if(isW) then
                        viscWtemp = viscosity( xWtemp, xiWbary(localnrows+1,i), Pw(localnrows,i), 'l' )
                    end if
                    if(isN) then
                        viscNtemp = viscosity( xNtemp, xiNbary(localnrows+1,i), Pw(localnrows,i), 'g' )
                    end if
                end if
                if(Uwy(localnrows+1,i) < 0) then
                    lambdawy(localnrows+1,i) = Kyybar(localnrows+1,i)*computekr_W(Swtemp)/viscWtemp
                elseif(Uwy(localnrows+1,i) > 0) then
                    lambdawy(localnrows+1,i) = Kyybar(localnrows+1,i)*computekr_W(Sw(localnrows,i))/viscW(localnrows,i)
                    xiWbary(localnrows+1,i) = xiW(localnrows,i)
                    densiWbary(localnrows+1,i) = densiW(localnrows,i)
                    xWbary(1:Nc,localnrows+1,i) = xW(1:Nc,localnrows,i)
                elseif(isDiriY(xlower+i-1,2) == 1) then
                    lambdawy(localnrows+1,i) = Kyybar(localnrows+1,i)*computekr_W(Sw(localnrows,i))/viscW(localnrows,i)
                    xiWbary(localnrows+1,i) = xiW(localnrows,i)
                    densiWbary(localnrows+1,i) = densiW(localnrows,i)
                    xWbary(1:Nc,localnrows+1,i) = xW(1:Nc,localnrows,i)
                end if
                if(Uny(localnrows+1,i) < 0) then
                    lambdany(localnrows+1,i) = Kyybar(localnrows+1,i)*computekr_N(Swtemp)/viscNtemp
                elseif(Uny(localnrows+1,i) > 0) then
                    lambdany(localnrows+1,i) = Kyybar(localnrows+1,i)*computekr_N(Sw(localnrows,i))/viscN(localnrows,i)
                    xiNbary(localnrows+1,i) = xiN(localnrows,i)
                    densiNbary(localnrows+1,i) = densiN(localnrows,i)
                    xNbary(1:Nc,localnrows+1,i) = xN(1:Nc,localnrows,i)
                elseif(isDiriY(xlower+i-1,2) == 1) then
                    lambdany(localnrows+1,i) = Kyybar(localnrows+1,i)*computekr_N(Sw(localnrows,i))/viscN(localnrows,i)
                    xiNbary(localnrows+1,i) = xiN(localnrows,i)
                    densiNbary(localnrows+1,i) = densiN(localnrows,i)
                    xNbary(1:Nc,localnrows+1,i) = xN(1:Nc,localnrows,i)
                end if
            end do
            do j = 1, localncols
                do i = 1, localnrows
                    if(Uwy(i,j) > 0) then
                        lambdawy(i,j) = Kyybar(i,j)*computekr_W(Sw(i-1,j))/viscW(i-1,j)
                        xiWbary(i,j) = xiW(i-1,j)
                        densiWbary(i,j) = densiW(i-1,j)
                        xWbary(1:Nc,i,j) = xW(1:Nc,i-1,j)
                    else
                        lambdawy(i,j) = Kyybar(i,j)*computekr_W(Sw(i,j))/viscW(i,j)
                        xiWbary(i,j) = xiW(i,j)
                        densiWbary(i,j) = densiW(i,j)
                        xWbary(1:Nc,i,j) = xW(1:Nc,i,j)
                    end if
                    if(Uny(i,j) > 0) then
                        lambdany(i,j) = Kyybar(i,j)*computekr_N(Sw(i-1,j))/viscN(i-1,j)
                        xiNbary(i,j) = xiN(i-1,j)
                        densiNbary(i,j) = densiN(i-1,j)
                        xNbary(1:Nc,i,j) = xN(1:Nc,i-1,j)
                    else
                        lambdany(i,j) = Kyybar(i,j)*computekr_N(Sw(i,j))/viscN(i,j)
                        xiNbary(i,j) = xiN(i,j)
                        densiNbary(i,j) = densiN(i,j)
                        xNbary(1:Nc,i,j) = xN(1:Nc,i,j)
                    end if
                end do
            end do
        else
            do j = 1, localncols
                do i = 1, localnrows+1
                    if(Uwy(i,j) > 0) then
                        lambdawy(i,j) = Kyybar(i,j)*computekr_W(Sw(i-1,j))/viscW(i-1,j)
                        xiWbary(i,j) = xiW(i-1,j)
                        densiWbary(i,j) = densiW(i-1,j)
                        xWbary(1:Nc,i,j) = xW(1:Nc,i-1,j)
                    else
                        lambdawy(i,j) = Kyybar(i,j)*computekr_W(Sw(i,j))/viscW(i,j)
                        xiWbary(i,j) = xiW(i,j)
                        densiWbary(i,j) = densiW(i,j)
                        xWbary(1:Nc,i,j) = xW(1:Nc,i,j)
                    end if
                    if(Uny(i,j) > 0) then
                        lambdany(i,j) = Kyybar(i,j)*computekr_N(Sw(i-1,j))/viscN(i-1,j)
                        xiNbary(i,j) = xiN(i-1,j)
                        densiNbary(i,j) = densiN(i-1,j)
                        xNbary(1:Nc,i,j) = xN(1:Nc,i-1,j)
                    else
                        lambdany(i,j) = Kyybar(i,j)*computekr_N(Sw(i,j))/viscN(i,j)
                        xiNbary(i,j) = xiN(i,j)
                        densiNbary(i,j) = densiN(i,j)
                        xNbary(1:Nc,i,j) = xN(1:Nc,i,j)
                    end if
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

        real(kind=8) :: xedge, yedge, ledge, redge, uedge, dedge
        real(kind=8) :: cotwx2, cotnx2, cotwx1, cotnx1, cotwy1, cotny1, cotwy2, cotny2
        real(kind=8) :: up, down, left, right
        real(kind=8) :: sum1, sum2
        real(kind=8), dimension(:), pointer :: Pwsent
        real(kind=8), dimension(:), pointer :: sentbuffer, recvbuffer
        integer :: sentbuffersize, recvbuffersize
        integer :: i, j, m, r, k

        integer :: ierr, errorcode, num_iter
        integer :: status(MPI_STATUS_SIZE)
        integer :: request, requestl, requestr, requestu, requestd
        integer :: buffer_size = MAX_BUF
        real(kind=8) :: buffer(MAX_BUF)
        real(kind=8) :: solvertimestart, solvertimefinish
        real(kind=8) :: commtimestart, commtimefinish
        real(kind=8), dimension(:), pointer :: values
        real(kind=8), dimension(:), pointer :: rhs_values, x_values

        allocate(rhs_values(local_size))
        allocate(x_values(local_size))
        allocate(values(local_size*5))
        values = 0

        k = 1
        r = 1
        do i = 1, localnrows
            do j = 1, localncols
            
                xedge = xs(xlower+j) - xs(xlower+j-1)
                yedge = ys(ylower+i) - ys(ylower+i-1)

                if((pcol==1).and.(j==1)) then
                    ledge = 0
                else
                    ledge = xs(xlower+j-1) - xs(xlower+j-2)
                end if

                if((pcol==pncols).and.(j==localncols)) then
                    redge = 0
                else
                    redge = xs(xlower+j+1) - xs(xlower+j)
                end if

                if((prow==1).and.(i==1)) then
                    dedge = 0
                else
                    dedge = ys(ylower+i-1) - ys(ylower+i-2)
                end if

                if((prow==pnrows).and.(i==localnrows)) then
                    uedge = 0
                else
                    uedge = ys(ylower+i+1) - ys(ylower+i)
                end if

                values(r) = poro(xlower+j-1,ylower+i-1)*Cf(i,j)/(timeEnd/nt)

                rhs_values(k) = poro(xlower+j-1,ylower+i-1)*Cf(i,j)/(timeEnd/nt)*Pw(i,j)
                do m = 1, Nc
                    rhs_values(k) = rhs_values(k) + v(m,i,j)*src(m,xlower+j-1,ylower+i-1)
                end do

                cotwx2 = 0
                cotnx2 = 0
                cotwx1 = 0
                cotnx1 = 0
                cotwy1 = 0
                cotny1 = 0
                cotwy2 = 0
                cotny2 = 0
                do m = 1, Nc
                    cotwx2 = cotwx2 + lambdawx(i,j)*xiWbarx(i,j)*xWbarx(m,i,j)*v(m,i,j)
                    cotnx2 = cotnx2 + lambdanx(i,j)*xiNbarx(i,j)*xNbarx(m,i,j)*v(m,i,j)
                    cotwx1 = cotwx1 + lambdawx(i,j+1)*xiWbarx(i,j+1)*xWbarx(m,i,j+1)*v(m,i,j)
                    cotnx1 = cotnx1 + lambdanx(i,j+1)*xiNbarx(i,j+1)*xNbarx(m,i,j+1)*v(m,i,j)
                    cotwy1 = cotwy1 + lambdawy(i+1,j)*xiWbary(i+1,j)*xWbary(m,i+1,j)*v(m,i,j)
                    cotny1 = cotny1 + lambdany(i+1,j)*xiNbary(i+1,j)*xNbary(m,i+1,j)*v(m,i,j)
                    cotwy2 = cotwy2 + lambdawy(i,j)*xiWbary(i,j)*xWbary(m,i,j)*v(m,i,j)
                    cotny2 = cotny2 + lambdany(i,j)*xiNbary(i,j)*xNbary(m,i,j)*v(m,i,j)
                end do

                up = -2*(cotwy1+cotny1)/yedge/(yedge+uedge)
                down = -2*(cotwy2+cotny2)/yedge/(yedge+dedge)
                left = -2*(cotwx2+cotnx2)/xedge/(xedge+ledge)
                right = -2*(cotwx1+cotnx1)/xedge/(xedge+redge)

                if((prow == 1).and.(i == 1).and.(isDiriY(xlower+j-1,1) == 0)) then
                    sum1 = 0
                    sum2 = 0
                    do m = 1, Nc
                        sum1 = sum1 + xWbary(m,i,j)*v(m,i,j)
                        sum2 = sum2 + xNbary(m,i,j)*v(m,i,j)
                    end do
                    rhs_values(k) = rhs_values(k) + UwBdryY(xlower+j-1,1)*xiWbary(i,j)*sum1/yedge + &!
                        UnBdryY(xlower+j-1,1)*xiNbary(i,j)*sum2/yedge
                elseif((prow == 1).and.(i == 1).and.(isDiriY(xlower+j-1,1) == 1)) then
                    rhs_values(k) = rhs_values(k) - down*Pw(i-1,j) + cotwy2*densiWbary(i,j)* &!
                        gravY/yedge + cotny2*densiNbary(i,j)*gravY/yedge
                    values(r) = values(r) - down
                else
                    values(r+3) = down
                    values(r) = values(r) - down
                    rhs_values(k) = rhs_values(k) + cotwy2*densiWbary(i,j)*gravY/yedge + &!
                        cotny2*densiNbary(i,j)*gravY/yedge
                end if

                if((pcol == 1).and.(j == 1).and.(isDiriX(1,ylower+i-1) == 0)) then
                    sum1 = 0
                    sum2 = 0
                    do m = 1, Nc
                        sum1 = sum1 + xWbarx(m,i,j)*v(m,i,j)
                        sum2 = sum2 + xNbarx(m,i,j)*v(m,i,j)
                    end do
                    rhs_values(k) = rhs_values(k) + UwBdryX(1,ylower+i-1)*xiWbarx(i,j)*sum1/xedge + &!
                        UnBdryX(1,ylower+i-1)*xiNbarx(i,j)*sum2/xedge
                else if((pcol == 1).and.(j == 1).and.(isDiriX(1,ylower+i-1) == 1)) then
                    rhs_values(k) = rhs_values(k) - left*Pw(i,j-1) + cotwx2*densiWbarx(i,j)* &!
                        gravX/xedge + cotnx2*densiNbarx(i,j)*gravX/xedge
                    values(r) = values(r) - left
                else
                    values(r+1) = left
                    values(r) = values(r) - left
                    rhs_values(k) = rhs_values(k) + cotwx2*densiWbarx(i,j)*gravX/xedge + &!
                        cotnx2*densiNbarx(i,j)*gravX/xedge
                end if

                if((pcol == pncols).and.(j == localncols).and.(isDiriX(2,ylower+i-1) == 0)) then
                    sum1 = 0
                    sum2 = 0
                    do m = 1, Nc
                        sum1 = sum1 + xWbarx(m,i,j+1)*v(m,i,j)
                        sum2 = sum2 + xNbarx(m,i,j+1)*v(m,i,j)
                    end do
                    rhs_values(k) = rhs_values(k) - UwBdryX(2,ylower+i-1)*xiWbarx(i,j+1)*sum1/xedge - &!
                        UnBdryX(2,ylower+i-1)*xiNbarx(i,j+1)*sum2/xedge
                else if((pcol == pncols).and.(j == localncols).and.(isDiriX(2,ylower+i-1) == 1)) then
                    rhs_values(k) = rhs_values(k) - right*Pw(i,j+1)-cotwx1*densiWbarx(i,j+1)* &!
                        gravX/xedge - cotnx1*densiNbarx(i,j+1)*gravX/xedge
                    values(r) = values(r) - right
                else
                    values(r+2) = right
                    values(r) = values(r) - right
                    rhs_values(k) = rhs_values(k) - cotwx1*densiWbarx(i,j+1)*gravX/xedge - &!
                        cotnx1*densiNbarx(i,j+1)*gravX/xedge
                end if

                if((prow == pnrows).and.(i == localnrows).and.(isDiriY(xlower+j-1,2) == 0)) then
                    sum1 = 0
                    sum2 = 0
                    do m = 1, Nc
                        sum1 = sum1 + xWbary(m,i+1,j)*v(m,i,j)
                        sum2 = sum2 + xNbary(m,i+1,j)*v(m,i,j)
                    end do
                    rhs_values(k) = rhs_values(k) - UwBdryY(xlower+j-1,2)*xiWbary(i+1,j)*sum1/yedge - &!
                        UnBdryY(xlower+j-1,2)*xiNbary(i+1,j)*sum2/yedge
                else if((prow == pnrows).and.(i == localnrows).and.(isDiriY(xlower+j-1,2) == 1)) then
                    rhs_values(k) = rhs_values(k) - up*Pw(i+1,j) - cotwy1*densiWbary(i+1,j)*gravY&!
                        /yedge - cotny1*densiNbary(i+1,j)*gravY/yedge
                    values(r) = values(r) - up
                else
                    values(r+4) = up
                    values(r) = values(r) - up
                    rhs_values(k) = rhs_values(k) - cotwy1*densiWbary(i+1,j)*gravY/yedge - &!
                        cotny1*densiNbary(i+1,j)*gravY/yedge
                end if

                r = r + 5
                k = k + 1

            end do
        end do

        call HYPRE_StructMatrixSetBoxValues(global_A, ilower, iupper, 5, stencil_indices, values, ierr)
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

        k = 1
        do i = 1, localnrows
            do j = 1, localncols
                Pw(i,j) = x_values(k)
                k = k + 1
            end do
        end do

        initial_x_guess(:) = x_values(:)

        deallocate(rhs_values)
        deallocate(x_values)
        deallocate(values)

        commtimestart = MPI_Wtime()

        call MPI_BUFFER_ATTACH(buffer,buffer_size,ierr)

        if(pcol /= 1) then
            allocate(Pwsent(localnrows))
            Pwsent = Pw(1:localnrows,1)
            call MPI_IBSEND(Pwsent, localnrows, MPI_DOUBLE_PRECISION, myid-pnrows, myid, MPI_COMM_WORLD, requestl, ierr)
            deallocate(Pwsent)
        end if

        if(pcol /= pncols) then
            allocate(Pwsent(localnrows))
            Pwsent = Pw(1:localnrows,localncols)
            call MPI_IBSEND(Pwsent, localnrows, MPI_DOUBLE_PRECISION, myid+pnrows, myid, MPI_COMM_WORLD, requestr, ierr)
            deallocate(Pwsent)
        end if

        if(prow /= 1) then
            allocate(Pwsent(localncols))
            Pwsent = Pw(1,1:localncols)
            call MPI_IBSEND(Pwsent, localncols, MPI_DOUBLE_PRECISION, myid-1, myid, MPI_COMM_WORLD, requestd, ierr)
            deallocate(Pwsent)
        end if

        if(prow /= pnrows) then
            allocate(Pwsent(localncols))
            Pwsent = Pw(localnrows,1:localncols)
            call MPI_IBSEND(Pwsent, localncols, MPI_DOUBLE_PRECISION, myid+1, myid, MPI_COMM_WORLD, requestu, ierr)
            deallocate(Pwsent)
        end if

        if(pcol /= 1) then
            recvbuffersize = localnrows
            allocate(recvbuffer(recvbuffersize))
            call MPI_RECV(recvbuffer, recvbuffersize, MPI_DOUBLE_PRECISION, myid-pnrows, myid-pnrows, &!
                MPI_COMM_WORLD, status, ierr)
            Pw(1:localnrows,0) = recvbuffer(1:localnrows)
            deallocate(recvbuffer)
        end if

        if(pcol /= pncols) then
            recvbuffersize = localnrows
            allocate(recvbuffer(recvbuffersize))
            call MPI_RECV(recvbuffer, recvbuffersize, MPI_DOUBLE_PRECISION, myid+pnrows, myid+pnrows, &!
                MPI_COMM_WORLD, status, ierr)
            Pw(1:localnrows,localncols+1) = recvbuffer(1:localnrows)
            deallocate(recvbuffer)
        end if

        if(prow /= 1) then
            recvbuffersize = localncols
            allocate(recvbuffer(recvbuffersize))
            call MPI_RECV(recvbuffer, recvbuffersize, MPI_DOUBLE_PRECISION, myid-1, myid-1, &!
                MPI_COMM_WORLD, status, ierr)
            Pw(0,1:localncols) = recvbuffer(1:localncols)
            deallocate(recvbuffer)
        end if

        if(prow /= pnrows) then
            recvbuffersize = localncols
            allocate(recvbuffer(recvbuffersize))
            call MPI_RECV(recvbuffer, recvbuffersize, MPI_DOUBLE_PRECISION, myid+1, myid+1, &!
                MPI_COMM_WORLD, status, ierr)
            Pw(localnrows+1,1:localncols) = recvbuffer(1:localncols)
            deallocate(recvbuffer)
        end if

        if(pcol /= 1) then
            call MPI_WAIT(requestl, status, ierr)
        end if
        if(pcol /= pncols) then
            call MPI_WAIT(requestr, status, ierr)
        end if
        if(prow /= pnrows) then
            call MPI_WAIT(requestu, status, ierr)
        end if
        if(prow /= 1) then
            call MPI_WAIT(requestd, status, ierr)
        end if

        call MPI_BUFFER_DETACH(buffer,buffer_size,ierr)

        commtimefinish = MPI_Wtime()
        commtime = commtime + commtimefinish - commtimestart

    end subroutine computePres

    subroutine computeVel()

        implicit none
       
        integer :: i, j, ierr

        ! compute the total new velocities in x direction
        if(pncols==1) then
            do i = 1, localnrows
                if(isDiriX(1,ylower+i-1) == 1) then
                    Uwx(i,1) = -lambdawx(i,1)*((Pw(i,1)-Pw(i,0))*2/(xs(2)-xs(1)) - densiWbarx(i,1)*gravX)
                    Unx(i,1) = -lambdanx(i,1)*((Pw(i,1)-Pw(i,0))*2/(xs(2)-xs(1)) - densiNbarx(i,1)*gravX)
                else
                    Uwx(i,1) = UwBdryX(1,ylower+i-1)
                    Unx(i,1) = UnBdryX(1,ylower+i-1)
                end if
            end do
            do i = 1, localnrows
                if(isDiriX(2,ylower+i-1) == 1) then
                    Uwx(i,localncols+1) = -lambdawx(i,localncols+1)*((Pw(i,localncols+1)-Pw(i,localncols))*2/(xs(nx+1)-xs(nx)) &!
                        - densiWbarx(i,localncols+1)*gravX)
                    Unx(i,localncols+1) = -lambdanx(i,localncols+1)*((Pw(i,localncols+1)-Pw(i,localncols))*2/(xs(nx+1)-xs(nx)) &!
                        - densiNbarx(i,localncols+1)*gravX)
                    if(Uwx(i,localncols+1) < 0) then
                        Uwx(i,localncols+1) = 0
                    end if
                    if(Unx(i,localncols+1) < 0) then
                        Unx(i,localncols+1) = 0
                    end if
                else
                    Uwx(i,localncols+1) = UwBdryX(2,ylower+i-1)
                    Unx(i,localncols+1) = UnBdryX(2,ylower+i-1)
                end if
            end do
            do j = 2, localncols
                do i = 1, localnrows
                    Uwx(i,j) = -lambdawx(i,j)*((Pw(i,j)-Pw(i,j-1))*2/(xs(xlower+j)-xs(xlower+j-2)) &!
                        - densiWbarx(i,j)*gravX)
                    Unx(i,j) = -lambdanx(i,j)*((Pw(i,j)-Pw(i,j-1))*2/(xs(xlower+j)-xs(xlower+j-2)) &!
                        - densiNbarx(i,j)*gravX)
                end do
            end do
        elseif(pcol==1) then
            do i = 1, localnrows
                if(isDiriX(1,ylower+i-1) == 1) then
                    Uwx(i,1) = -lambdawx(i,1)*((Pw(i,1)-Pw(i,0))*2/(xs(2)-xs(1)) - densiWbarx(i,1)*gravX)
                    Unx(i,1) = -lambdanx(i,1)*((Pw(i,1)-Pw(i,0))*2/(xs(2)-xs(1)) - densiNbarx(i,1)*gravX)
                else
                    Uwx(i,1) = UwBdryX(1,ylower+i-1)
                    Unx(i,1) = UnBdryX(1,ylower+i-1)
                end if
            end do
            do j = 2, localncols+1
                do i = 1, localnrows
                    Uwx(i,j) = -lambdawx(i,j)*((Pw(i,j)-Pw(i,j-1))*2/(xs(xlower+j)-xs(xlower+j-2)) &!
                        - densiWbarx(i,j)*gravX)
                    Unx(i,j) = -lambdanx(i,j)*((Pw(i,j)-Pw(i,j-1))*2/(xs(xlower+j)-xs(xlower+j-2)) &!
                        - densiNbarx(i,j)*gravX)
                end do
            end do
        elseif(pcol==pncols) then
            do i = 1, localnrows
                if(isDiriX(2,ylower+i-1) == 1) then
                    Uwx(i,localncols+1) = -lambdawx(i,localncols+1)*((Pw(i,localncols+1)-Pw(i,localncols))*2/(xs(nx+1)-xs(nx)) &!
                        - densiWbarx(i,localncols+1)*gravX)
                    Unx(i,localncols+1) = -lambdanx(i,localncols+1)*((Pw(i,localncols+1)-Pw(i,localncols))*2/(xs(nx+1)-xs(nx)) &!
                        - densiNbarx(i,localncols+1)*gravX)
                    if(Uwx(i,localncols+1) < 0) then
                        Uwx(i,localncols+1) = 0
                    end if
                    if(Unx(i,localncols+1) < 0) then
                        Unx(i,localncols+1) = 0
                    end if
                else
                    Uwx(i,localncols+1) = UwBdryX(2,ylower+i-1)
                    Unx(i,localncols+1) = UnBdryX(2,ylower+i-1)
                end if
            end do
            do j = 1, localncols
                do i = 1, localnrows
                    Uwx(i,j) = -lambdawx(i,j)*((Pw(i,j)-Pw(i,j-1))*2/(xs(xlower+j)-xs(xlower+j-2)) &!
                        - densiWbarx(i,j)*gravX)
                    Unx(i,j) = -lambdanx(i,j)*((Pw(i,j)-Pw(i,j-1))*2/(xs(xlower+j)-xs(xlower+j-2)) &!
                        - densiNbarx(i,j)*gravX)
                end do
            end do
        else
            do j = 1, localncols+1
                do i = 1, localnrows
                    Uwx(i,j) = -lambdawx(i,j)*((Pw(i,j)-Pw(i,j-1))*2/(xs(xlower+j)-xs(xlower+j-2)) &!
                        - densiWbarx(i,j)*gravX)
                    Unx(i,j) = -lambdanx(i,j)*((Pw(i,j)-Pw(i,j-1))*2/(xs(xlower+j)-xs(xlower+j-2)) &!
                        - densiNbarx(i,j)*gravX)
                end do
            end do
        end if

        ! compute the total new velocities in y direction
        if(pnrows==1) then
            do j = 1, localncols
                if(isDiriY(xlower+j-1,1) == 1) then
                    Uwy(1,j) = -lambdawy(1,j)*((Pw(1,j)-Pw(0,j))*2/(ys(2)-ys(1)) - densiWbary(1,j)*gravY)
                    Uny(1,j) = -lambdany(1,j)*((Pw(1,j)-Pw(0,j))*2/(ys(2)-ys(1)) - densiNbary(1,j)*gravY)
                else
                    Uwy(1,j) = UwBdryY(xlower+j-1,1)
                    Uny(1,j) = UnBdryY(xlower+j-1,1)
                end if
            end do
            do j = 1, localncols
                if(isDiriY(xlower+j-1,2) == 1) then
                    Uwy(localnrows+1,j) = -lambdawy(localnrows+1,j)*((Pw(localnrows+1,j)-Pw(localnrows,j))*2/(ys(ny+1)-ys(ny)) &!
                        - densiWbary(localnrows+1,j)*gravY)
                    Uny(localnrows+1,j) = -lambdany(localnrows+1,j)*((Pw(localnrows+1,j)-Pw(localnrows,j))*2/(ys(ny+1)-ys(ny)) &!
                        - densiNbary(localnrows+1,j)*gravY)
                else
                    Uwy(localnrows+1,j) = UwBdryY(xlower+j-1,2)
                    Uny(localnrows+1,j) = UnBdryY(xlower+j-1,2)
                end if
            end do
            do j = 1, localncols
                do i = 2, localnrows
                    Uwy(i,j) = -lambdawy(i,j)*((Pw(i,j)-Pw(i-1,j))*2/(ys(ylower+i)-ys(ylower+i-2)) &!
                        - densiWbary(i,j)*gravY)
                    Uny(i,j) = -lambdany(i,j)*((Pw(i,j)-Pw(i-1,j))*2/(ys(ylower+i)-ys(ylower+i-2)) &!
                        - densiNbary(i,j)*gravY)
                end do
            end do
        elseif(prow==1) then
            do j = 1, localncols
                if(isDiriY(xlower+j-1,1) == 1) then
                    Uwy(1,j) = -lambdawy(1,j)*((Pw(1,j)-Pw(0,j))*2/(ys(2)-ys(1)) - densiWbary(1,j)*gravY)
                    Uny(1,j) = -lambdany(1,j)*((Pw(1,j)-Pw(0,j))*2/(ys(2)-ys(1)) - densiNbary(1,j)*gravY)
                else
                    Uwy(1,j) = UwBdryY(xlower+j-1,1)
                    Uny(1,j) = UnBdryY(xlower+j-1,1)
                end if
            end do
            do j = 1, localncols
                do i = 2, localnrows+1
                    Uwy(i,j) = -lambdawy(i,j)*((Pw(i,j)-Pw(i-1,j))*2/(ys(ylower+i)-ys(ylower+i-2)) &!
                        - densiWbary(i,j)*gravY)
                    Uny(i,j) = -lambdany(i,j)*((Pw(i,j)-Pw(i-1,j))*2/(ys(ylower+i)-ys(ylower+i-2)) &!
                        - densiNbary(i,j)*gravY)
                end do
            end do
        elseif(prow==pnrows) then
            do j = 1, localncols
                if(isDiriY(xlower+j-1,2) == 1) then
                    Uwy(localnrows+1,j) = -lambdawy(localnrows+1,j)*((Pw(localnrows+1,j)-Pw(localnrows,j))*2/(ys(ny+1)-ys(ny)) &!
                        - densiWbary(localnrows+1,j)*gravY)
                    Uny(localnrows+1,j) = -lambdany(localnrows+1,j)*((Pw(localnrows+1,j)-Pw(localnrows,j))*2/(ys(ny+1)-ys(ny)) &!
                        - densiNbary(localnrows+1,j)*gravY)
                else
                    Uwy(localnrows+1,j) = UwBdryY(xlower+j-1,2)
                    Uny(localnrows+1,j) = UnBdryY(xlower+j-1,2)
                end if
            end do
            do j = 1, localncols
                do i = 1, localnrows
                    Uwy(i,j) = -lambdawy(i,j)*((Pw(i,j)-Pw(i-1,j))*2/(ys(ylower+i)-ys(ylower+i-2)) &!
                        - densiWbary(i,j)*gravY)
                    Uny(i,j) = -lambdany(i,j)*((Pw(i,j)-Pw(i-1,j))*2/(ys(ylower+i)-ys(ylower+i-2)) &!
                        - densiNbary(i,j)*gravY)
                end do
            end do
        else
            do j = 1, localncols
                do i = 1, localnrows+1
                    Uwy(i,j) = -lambdawy(i,j)*((Pw(i,j)-Pw(i-1,j))*2/(ys(ylower+i)-ys(ylower+i-2)) &!
                        - densiWbary(i,j)*gravY)
                    Uny(i,j) = -lambdany(i,j)*((Pw(i,j)-Pw(i-1,j))*2/(ys(ylower+i)-ys(ylower+i-2)) &!
                        - densiNbary(i,j)*gravY)
                end do
            end do
        end if

    end subroutine computeVel

    subroutine computeMoleFrac()

        implicit none
        real(kind=8), dimension(:), pointer :: xic
        real(kind=8) :: right, left, up, down, div
        real(kind=8) :: xinew
        real(kind=8), dimension(:), pointer :: zsent, recvbuffer
        integer :: recvbuffersize
        integer :: requestl, requestr, requestu, requestd
        integer :: status(MPI_STATUS_SIZE)
        integer :: i, j, m, k

        integer :: ierr, errorcode
        integer :: buffer_size = MAX_BUF
        real(kind=8) :: buffer(MAX_BUF)
        real(kind=8) :: commtimestart, commtimefinish

        allocate(xic(Nc))

        do j = 1, localncols
            do i = 1, localnrows
                do m = 1, Nc
                    right = Uwx(i,j+1)*xiWbarx(i,j+1)*xWbarx(m,i,j+1) + Unx(i,j+1)* &!
                        xiNbarx(i,j+1)*xNbarx(m,i,j+1)
                    left = Uwx(i,j)*xiWbarx(i,j)*xWbarx(m,i,j) + &!
                        Unx(i,j)*xiNbarx(i,j)*xNbarx(m,i,j)
                    up = Uwy(i+1,j)*xiWbary(i+1,j)*xWbary(m,i+1,j) + Uny(i+1,j)* &!
                        xiNbary(i+1,j)*xNbary(m,i+1,j)
                    down = Uwy(i,j)*xiWbary(i,j)*xWbary(m,i,j) + &!
                        Uny(i,j)*xiNbary(i,j)*xNbary(m,i,j)
                    div = (right-left)/(xs(xlower+j)-xs(xlower+j-1)) + (up-down)/(ys(ylower+i)-ys(ylower+i-1))
                    xic(m) = (src(m,xlower+j-1,ylower+i-1)-div)*(timeEnd/nt)/poro(xlower+j-1,ylower+i-1) &!
                        + z(m,i,j)*(Sw(i,j)*xiW(i,j)+(1-Sw(i,j))*xiN(i,j))
                    if(xic(m)<0) then
                        print *, 'Please tune the time step.', xic(m),myid,m,i,j,Uwy(i+1,j),Uny(i+1,j), &!
                            Uwy(i,j),Uny(i,j), Uwx(i+1,j),Unx(i+1,j),Uwx(i,j),Unx(i,j),xiWbarx(i,j+1), &!
                            xiWbarx(i,j),xiWbary(i,j+1),xiWbary(i,j),xiNbarx(i,j+1), xiNbarx(i,j),xiNbary(i,j+1), &!
                            xiNbary(i,j)
                        call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
                    end if
                end do
                xinew = 0
                do m = 1, Nc
                    xinew = xinew + xic(m)
                end do
                z(1:Nc,i,j) = xic(1:Nc)/xinew
                do m = 1, Nc
                    if(z(m,i,j) < 1.D-99) then
                        z(m,i,j) = 0.0
                    end if
                end do
            end do
        end do

        commtimestart = MPI_Wtime()
        call MPI_BUFFER_ATTACH(buffer,buffer_size,ierr)

        if(pcol /= 1) then
            allocate(zsent(Nc*localnrows))
            k = 1
            do i = 1, localnrows
                do m = 1, Nc
                    zsent(k) = z(m,i,1)
                    k = k + 1
                end do
            end do
            call MPI_IBSEND(zsent, Nc*localnrows, MPI_DOUBLE_PRECISION, myid-pnrows, myid, MPI_COMM_WORLD, requestl, ierr)
            deallocate(zsent)
        end if

        if(pcol /= pncols) then
            allocate(zsent(Nc*localnrows))
            k = 1
            do i = 1, localnrows
                do m = 1, Nc
                    zsent(k) = z(m,i,localncols)
                    k = k + 1
                end do
            end do
            call MPI_IBSEND(zsent, Nc*localnrows, MPI_DOUBLE_PRECISION, myid+pnrows, myid, MPI_COMM_WORLD, requestr, ierr)
            deallocate(zsent)
        end if

        if(prow /= 1) then
            allocate(zsent(Nc*localncols))
            k = 1
            do i = 1, localncols
                do m = 1, Nc
                    zsent(k) = z(m,1,i)
                    k = k + 1
                end do
            end do
            call MPI_IBSEND(zsent, Nc*localncols, MPI_DOUBLE_PRECISION, myid-1, myid, MPI_COMM_WORLD, requestd, ierr)
            deallocate(zsent)
        end if

        if(prow /= pnrows) then
            allocate(zsent(Nc*localncols))
            k = 1
            do i = 1, localncols
                do m = 1, Nc
                    zsent(k) = z(m,localnrows,i)
                    k = k + 1
                end do
            end do
            call MPI_IBSEND(zsent, Nc*localncols, MPI_DOUBLE_PRECISION, myid+1, myid, MPI_COMM_WORLD, requestu, ierr)
            deallocate(zsent)
        end if

        if(pcol /= 1) then
            recvbuffersize = Nc*localnrows
            allocate(recvbuffer(recvbuffersize))
            call MPI_RECV(recvbuffer, recvbuffersize, MPI_DOUBLE_PRECISION, myid-pnrows, myid-pnrows, &!
                MPI_COMM_WORLD, status, ierr)
            k = 1
            do i = 1, localnrows
                do m = 1, Nc
                    z(m,i,0) = recvbuffer(k)
                    k = k + 1
                end do
            end do
            deallocate(recvbuffer)
        end if

        if(pcol /= pncols) then
            recvbuffersize = Nc*localnrows
            allocate(recvbuffer(recvbuffersize))
            call MPI_RECV(recvbuffer, recvbuffersize, MPI_DOUBLE_PRECISION, myid+pnrows, myid+pnrows, &!
                MPI_COMM_WORLD, status, ierr)
            k = 1
            do i = 1, localnrows
                do m = 1, Nc
                    z(m,i,localncols+1) = recvbuffer(k)
                    k = k + 1
                end do
            end do
            deallocate(recvbuffer)
        end if

        if(prow /= 1) then
            recvbuffersize = Nc*localncols
            allocate(recvbuffer(recvbuffersize))
            call MPI_RECV(recvbuffer, recvbuffersize, MPI_DOUBLE_PRECISION, myid-1, myid-1, &!
                MPI_COMM_WORLD, status, ierr)
            k = 1
            do i = 1, localncols
                do m = 1, Nc
                    z(m,0,i) = recvbuffer(k)
                    k = k + 1
                end do
            end do
            deallocate(recvbuffer)
        end if

        if(prow /= pnrows) then
            recvbuffersize = Nc*localncols
            allocate(recvbuffer(recvbuffersize))
            call MPI_RECV(recvbuffer, recvbuffersize, MPI_DOUBLE_PRECISION, myid+1, myid+1, &!
                MPI_COMM_WORLD, status, ierr)
            k = 1
            do i = 1, localncols
                do m = 1, Nc
                    z(m,localnrows+1,i) = recvbuffer(k)
                    k = k + 1
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
        if(prow /= pnrows) then
            call MPI_WAIT(requestu, status, ierr)
        end if
        if(prow /= 1) then
            call MPI_WAIT(requestd, status, ierr)
        end if

        call MPI_BUFFER_DETACH(buffer,buffer_size,ierr)

        commtimefinish = MPI_Wtime()
        commtime = commtime + commtimefinish - commtimestart

        deallocate(xic)

    end subroutine computeMoleFrac

    subroutine finalize()

        implicit none

        integer :: m, ierr
        real(kind=8) :: timefinish

        deallocate(xs)
        deallocate(ys)
        deallocate(ts)
        deallocate(Kxx)
        deallocate(Kyy)
        deallocate(poro)
        deallocate(src)
        deallocate(isDiriX)
        deallocate(isDiriY)
        deallocate(PwBdryX)
        deallocate(PwBdryY)
        deallocate(PwInit)
        deallocate(zBdryX)
        deallocate(zBdryY)
        deallocate(zInit)
        deallocate(UwBdryX)
        deallocate(UwBdryY)
        deallocate(UnBdryX)
        deallocate(UnBdryY)
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
        deallocate(Unx)
        deallocate(Uny)
        deallocate(Sw)
        deallocate(lambdawx)
        deallocate(lambdawy)
        deallocate(lambdanx)
        deallocate(lambdany)
        deallocate(Kxxbar)
        deallocate(Kyybar)
        deallocate(z)
        deallocate(densiW)
        deallocate(densiN)
        deallocate(densiWbarx)
        deallocate(densiWbary)
        deallocate(densiNbarx)
        deallocate(densiNbary)
        deallocate(xW)
        deallocate(xN)
        deallocate(xWbarx)
        deallocate(xWbary)
        deallocate(xNbarx)
        deallocate(xNbary)
        deallocate(xiW)
        deallocate(xiN)
        deallocate(xiWbarx)
        deallocate(xiWbary)
        deallocate(xiNbarx)
        deallocate(xiNbary)
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
            do m = 1, Nc
                close(40+m)
            end do
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
