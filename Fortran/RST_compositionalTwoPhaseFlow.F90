
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
    use RST_flashcalculation_sparsegrid
#elif FULL_2c
    use RST_flashcalculation_fullgrid_2c
#elif FULL_3c
    use RST_flashcalculation_fullgrid_3c
#elif NN
    use RST_flashcalculation_nn
#else
    use RST_flashcalculation
#endif
    use RST_viscosity

    implicit none

contains

    subroutine initialize(modelCase)

        type(model), intent(in out) :: modelCase

        logical :: alive
        character(len=50) :: fmhtxt, fmrtxt, ftimetxt
        integer :: i, j, m

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

        allocate(Pw(ny+2, nx+2))
        allocate(Uwx(ny, nx+1))
        Uwx = 0.D0
        allocate(Uwy(ny+1, nx))
        Uwy = 0.D0
        allocate(Unx(ny, nx+1))
        Unx = 0.D0
        allocate(Uny(ny+1, nx))
        Uny = 0.D0
        allocate(Sw(ny, nx))
        allocate(lambdawx(ny, nx+1))
        allocate(lambdawy(ny+1, nx))
        allocate(lambdanx(ny, nx+1))
        allocate(lambdany(ny+1, nx))
        allocate(Kxxbar(ny, nx+1))
        allocate(Kyybar(ny+1, nx))
        allocate(z(Nc, ny, nx))
        allocate(densiW(ny, nx))
        allocate(densiN(ny, nx))
        allocate(densiWbarx(ny, nx+1))
        allocate(densiWbary(ny+1, nx))
        allocate(densiNbarx(ny, nx+1))
        allocate(densiNbary(ny+1, nx))
        allocate(xW(Nc, ny, nx))
        allocate(xN(Nc, ny, nx))
        allocate(xWbarx(Nc, ny, nx+1))
        allocate(xWbary(Nc, ny+1, nx))
        allocate(xNbarx(Nc, ny, nx+1))
        allocate(xNbary(Nc, ny+1, nx))
        allocate(xiW(ny, nx))
        allocate(xiN(ny, nx))
        allocate(xiWbarx(ny, nx+1))
        allocate(xiWbary(ny+1, nx))
        allocate(xiNbarx(ny, nx+1))
        allocate(xiNbary(ny+1, nx))
        allocate(v(Nc, ny, nx))
        allocate(Cf(ny, nx))
        allocate(viscW(ny, nx))
        allocate(viscN(ny, nx))

#if defined(FULL_2c) || defined(FULL_3c)
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
#elif SPARSE
        firstSG = .true.
#elif NN
        isFirstNN = .true.
#endif

        Pw(1, 2:nx+1) = PwBdryY(1:nx, 1)
        Pw(ny+2, 2:nx+1) = PwBdryY(1:nx, 2)
        Pw(2:ny+1, 1) = PwBdryX(1, 1:ny)
        Pw(2:ny+1, nx+2) = PwBdryX(2, 1:ny)
        Pw(2:ny+1,2:nx+1) = transpose(PwInit(1:nx, 1:ny))

        UwBdryX(1, 1:ny) = -UwBdryX(1, 1:ny)
        UwBdryY(1:nx, 1) = -UwBdryY(1:nx, 1)
        UnBdryX(1, 1:ny) = -UnBdryX(1, 1:ny)
        UnBdryY(1:nx, 1) = -UnBdryY(1:nx, 1)

        Uwx(1:ny, 1) = UwBdryX(1, 1:ny)
        Uwx(1:ny, nx+1) = UwBdryX(2, 1:ny)
        Uwy(1, 1:nx) = UwBdryY(1:nx, 1)
        Uwy(ny+1, 1:nx) = UwBdryY(1:nx, 2)

        Unx(1:ny, 1) = UnBdryX(1, 1:ny)
        Unx(1:ny, nx+1) = UnBdryX(2, 1:ny)
        Uny(1, 1:nx) = UnBdryY(1:nx, 1)
        Uny(ny+1, 1:nx) = UnBdryY(1:nx, 2)

        do i = 1, ny
            Kxxbar(i,1) = Kxx(1,i)
            Kxxbar(i,nx+1) = Kxx(nx,i)
        end do
        do i = 1, ny
            do j = 2, nx
                Kxxbar(i,j) = (xs(j+1)-xs(j-1)) / ((xs(j)-xs(j-1))/Kxx(j-1,i)+(xs(j+1)-xs(j))/Kxx(j,i))
            end do
        end do
        do j = 1, nx
            Kyybar(1,j) = Kyy(j,1)
            Kyybar(ny+1,j) = Kyy(j, ny)
        end do
        do i = 2, ny
            do j = 1, nx
                Kyybar(i,j) = (ys(i+1)-ys(i-1)) / ((ys(i)-ys(i-1))/Kyy(j,i-1)+(ys(i+1)-ys(i))/Kyy(j,i))
            end do
        end do

        do m = 1, Nc
            do i = 1, ny
                do j = 1, nx
                    z(m,i,j) = zInit(m,j,i)
                end do
            end do
        end do

        inquire(file = soludoc, exist = alive)
        if(.not.alive) then
            call system("mkdir "//trim(adjustl(soludoc)))
        end if

        fmhtxt = trim(adjustl(soludoc))//"/soln_cmp2PhFlw_moleHistory.txt"
        fmrtxt = trim(adjustl(soludoc))//"/soln_cmp2PhFlw_moleRatioHistory.txt"
        ftimetxt = trim(adjustl(soludoc))//"/soln_cmp2PhFlw_time.txt"

        open(unit=70, file=fmhtxt, status='replace')
        open(unit=80, file=fmrtxt, status='replace')
        open(unit=90, file=ftimetxt, status='replace')

#ifdef FULL_2c
        open(unit=120, file="Fullgrid_2c/xW1.txt", status='old')
        open(unit=130, file="Fullgrid_2c/xW2.txt", status='old')
        open(unit=140, file="Fullgrid_2c/xN1.txt", status='old')
        open(unit=150, file="Fullgrid_2c/xN2.txt", status='old')
        open(unit=160, file="Fullgrid_2c/xiW.txt", status='old')
        open(unit=170, file="Fullgrid_2c/xiN.txt", status='old')
        open(unit=180, file="Fullgrid_2c/densiW.txt", status='old')
        open(unit=190, file="Fullgrid_2c/densiN.txt", status='old')
        open(unit=200, file="Fullgrid_2c/sW.txt", status='old')
        open(unit=210, file="Fullgrid_2c/v1.txt", status='old')
        open(unit=220, file="Fullgrid_2c/v2.txt", status='old')
        open(unit=230, file="Fullgrid_2c/Cf.txt", status='old')
        open(unit=270, file="Fullgrid_2c/liquid.txt", status='old')
        open(unit=280, file="Fullgrid_2c/gas.txt", status='old')

        read(120,*) xtable(1,1:TABLESIZE)
        read(130,*) xtable(2,1:TABLESIZE)
        read(140,*) ytable(1,1:TABLESIZE)
        read(150,*) ytable(2,1:TABLESIZE)
        read(160,*) xiLtable(1:TABLESIZE)
        read(170,*) xiGtable(1:TABLESIZE)
        read(180,*) rhoLtable(1:TABLESIZE)
        read(190,*) rhoGtable(1:TABLESIZE)
        read(200,*) sLtable(1:TABLESIZE)
        read(210,*) vtable(1,1:TABLESIZE)
        read(220,*) vtable(2,1:TABLESIZE)
        read(230,*) Cftable(1:TABLESIZE)
        read(270,*) isWtable(1:TABLESIZE)
        read(280,*) isNtable(1:TABLESIZE)
#elif FULL_3c
        open(unit=120, file="Fullgrid_3c/xW1.txt", status='old')
        open(unit=130, file="Fullgrid_3c/xW2.txt", status='old')
        open(unit=140, file="Fullgrid_3c/xN1.txt", status='old')
        open(unit=150, file="Fullgrid_3c/xN2.txt", status='old')
        open(unit=160, file="Fullgrid_3c/xiW.txt", status='old')
        open(unit=170, file="Fullgrid_3c/xiN.txt", status='old')
        open(unit=180, file="Fullgrid_3c/densiW.txt", status='old')
        open(unit=190, file="Fullgrid_3c/densiN.txt", status='old')
        open(unit=200, file="Fullgrid_3c/sW.txt", status='old')
        open(unit=210, file="Fullgrid_3c/v1.txt", status='old')
        open(unit=220, file="Fullgrid_3c/v2.txt", status='old')
        open(unit=230, file="Fullgrid_3c/Cf.txt", status='old')
        open(unit=270, file="Fullgrid_3c/liquid.txt", status='old')
        open(unit=280, file="Fullgrid_3c/gas.txt", status='old')

        open(unit=240, file="Fullgrid_3c/xW3.txt", status='old')
        open(unit=250, file="Fullgrid_3c/xN3.txt", status='old')
        open(unit=260, file="Fullgrid_3c/v3.txt", status='old')

        read(120,*) xtable(1,1:TABLESIZE)
        read(130,*) xtable(2,1:TABLESIZE)
        read(140,*) ytable(1,1:TABLESIZE)
        read(150,*) ytable(2,1:TABLESIZE)
        read(160,*) xiLtable(1:TABLESIZE)
        read(170,*) xiGtable(1:TABLESIZE)
        read(180,*) rhoLtable(1:TABLESIZE)
        read(190,*) rhoGtable(1:TABLESIZE)
        read(200,*) sLtable(1:TABLESIZE)
        read(210,*) vtable(1,1:TABLESIZE)
        read(220,*) vtable(2,1:TABLESIZE)
        read(230,*) Cftable(1:TABLESIZE)
        read(270,*) isWtable(1:TABLESIZE)
        read(280,*) isNtable(1:TABLESIZE)

        read(240,*) xtable(3,1:TABLESIZE)
        read(250,*) ytable(3,1:TABLESIZE)
        read(260,*) vtable(3,1:TABLESIZE)
#endif

        call CPU_TIME(timebegin)

        totalmole = 0.D0
        t = 2

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

        real(kind=8), dimension(:), pointer :: xWtemp, xNtemp, ztemp, vtemp
        real(kind=8) :: ignore1, Swtemp, viscWtemp, viscNtemp
        real(kind=8), dimension(:,:,:), pointer :: moleincell
        real(kind=8), dimension(:), pointer :: leftmole
        real(kind=8) :: totaldesiredleftmole
        logical :: isW, isN, isRea
        integer :: i, j, m

        allocate(ztemp(Nc))
        allocate(xWtemp(Nc))
        allocate(xNtemp(Nc))
        allocate(vtemp(Nc))

        do j = 1, nx
            do i = 1, ny

                ztemp(1:Nc) = z(1:Nc,i,j)

#ifdef FULL_2c
                call flashcalculation_fullgrid_2c(Pw(i+1,j+1), ztemp, xWtemp, xNtemp, xiW(i,j), xiN(i,j), &!
                    densiW(i,j), densiN(i,j), Sw(i,j), vtemp, Cf(i,j), isW, isN)
#elif FULL_3c
                call flashcalculation_fullgrid_3c(Pw(i+1,j+1), ztemp, xWtemp, xNtemp, xiW(i,j), xiN(i,j), &!
                    densiW(i,j), densiN(i,j), Sw(i,j), vtemp, Cf(i,j), isW, isN)
#elif SPARSE
                call flashcalculation_sparsegrid( Pw(i+1,j+1), ztemp, xWtemp, xNtemp, xiW(i,j), xiN(i,j), &!
                    densiW(i,j), densiN(i,j), Sw(i,j), vtemp, Cf(i,j), isW, isN )
#elif NN
                call flashcalculation_nn( Pw(i+1,j+1), ztemp, xWtemp, xNtemp, xiW(i,j), xiN(i,j), &!
                    densiW(i,j), densiN(i,j), Sw(i,j), vtemp, Cf(i,j), isW, isN )
#else
                call flashcalculation( Pw(i+1,j+1), ztemp, xWtemp, xNtemp, xiW(i,j), xiN(i,j), &!
                    densiW(i,j), densiN(i,j), Sw(i,j), vtemp, Cf(i,j), isW, isN, isRea )
#endif

                xW(1:Nc,i,j) = xWtemp(1:Nc)
                xN(1:Nc,i,j) = xNtemp(1:Nc)
                v(1:Nc,i,j) = vtemp(1:Nc)

                if(isW) then
                    viscW(i,j) = viscosity( xWtemp, xiW(i,j), Pw(i+1,j+1), 'l' )
                else
                    viscW(i,j) = 1.D12
                end if

                if(isN) then
                    viscN(i,j) = viscosity( xNtemp, xiN(i,j), Pw(i+1,j+1), 'g' )
                else
                    viscN(i,j) = 1.D12
                end if
            end do
        end do

        allocate(moleincell(Nc, ny, nx))
        do j = 1, nx
            do i = 1, ny
                do m = 1, Nc
                    moleincell(m,i,j) = z(m,i,j)*(xiW(i,j)*Sw(i,j) + xiN(i,j)*(1-Sw(i,j)))&!
                        *(xs(j+1)-xs(j))*(ys(i+1)-ys(i))*poro(j,i)
                end do
            end do
        end do

        allocate(leftmole(Nc))
        leftmole = 0
        do m = 1, Nc
            do j = 1, nx
                do i = 1, ny
                    leftmole(m) = leftmole(m) + moleincell(m,i,j)
                end do
            end do
        end do

        totaldesiredleftmole = 0.D0
        do m = 2, Nc
            totaldesiredleftmole = totaldesiredleftmole + leftmole(m)
        end do

        if(t == 2) then
            totalmole = totaldesiredleftmole
        end if

        write(70, fmt="(es15.8)") (totalmole-totaldesiredleftmole)/totalmole
        print *, '% = ', (totalmole-totaldesiredleftmole)/totalmole
        write(80, fmt="(es15.8)") totaldesiredleftmole/leftmole(1)

        deallocate(leftmole)
        deallocate(moleincell)

        do i = 1, ny
            if((Uwx(i,1) > 0).or.(Unx(i,1) > 0)) then
                ztemp(1:Nc) = zBdryX(1:Nc,1,i)
#ifdef FULL_2c
                call flashcalculation_fullgrid_2c(Pw(i+1,2), ztemp, xWtemp, xNtemp, xiWbarx(i,1), xiNbarx(i,1), &!
                    densiWbarx(i,1), densiNbarx(i,1), Swtemp, vtemp, ignore1, isW, isN)
#elif FULL_3c
                call flashcalculation_fullgrid_3c(Pw(i+1,2), ztemp, xWtemp, xNtemp, xiWbarx(i,1), xiNbarx(i,1), &!
                    densiWbarx(i,1), densiNbarx(i,1), Swtemp, vtemp, ignore1, isW, isN)
#elif SPARSE
                call flashcalculation_sparsegrid(Pw(i+1,2), ztemp, xWtemp, xNtemp, xiWbarx(i,1), xiNbarx(i,1), &!
                    densiWbarx(i,1), densiNbarx(i,1), Swtemp, vtemp, ignore1, isW, isN)
#elif NN
                call flashcalculation_nn(Pw(i+1,2), ztemp, xWtemp, xNtemp, xiWbarx(i,1), xiNbarx(i,1), &!
                    densiWbarx(i,1), densiNbarx(i,1), Swtemp, vtemp, ignore1, isW, isN)
#else
                call flashcalculation(Pw(i+1,2), ztemp, xWtemp, xNtemp, xiWbarx(i,1), xiNbarx(i,1), &!
                    densiWbarx(i,1), densiNbarx(i,1), Swtemp, vtemp, ignore1, isW, isN, isRea)
#endif

                xWbarx(1:Nc,i,1) = xWtemp(1:Nc)
                xNbarx(1:Nc,i,1) = xNtemp(1:Nc)
                if(isW) then
                    viscWtemp = viscosity( xWtemp, xiWbarx(i,1), Pw(i+1,2), 'l' )
                end if
                if(isN) then
                    viscNtemp = viscosity( xNtemp, xiNbarx(i,1), Pw(i+1,2), 'g' )
                end if
            end if
            if(Uwx(i,1) > 0) then
                lambdawx(i,1) = Kxxbar(i,1)*computekr_W(Swtemp)/viscWtemp
            elseif(Uwx(i,1) < 0) then
                lambdawx(i,1) = Kxxbar(i,1)*computekr_W(Sw(i,1))/viscW(i,1)
                xiWbarx(i,1) = xiW(i,1)
                densiWbarx(i,1) = densiW(i,1)
                xWbarx(1:Nc,i,1) = xW(1:Nc,i,1)
            else
                if(isDiriX(1,i) == 1) then
                    lambdawx(i,1) = Kxxbar(i,1)*computekr_W(Sw(i,1))/viscW(i,1)
                    xiWbarx(i,1) = xiW(i,1)
                    densiWbarx(i,1) = densiW(i,1)
                    xWbarx(1:Nc,i,1) = xW(1:Nc,i,1)
                end if
            end if
            if (Unx(i,1) > 0) then
                lambdanx(i,1) = Kxxbar(i,1)*computekr_N(Swtemp)/viscNtemp
            elseif (Unx(i,1) < 0) then
                lambdanx(i,1) = Kxxbar(i,1)*computekr_N(Sw(i,1))/viscN(i,1)
                xiNbarx(i,1) = xiN(i,1)
                densiNbarx(i,1) = densiN(i,1)
                xNbarx(1:Nc,i,1) = xN(1:Nc,i,1)
            else
                if(isDiriX(1,i) == 1) then
                    lambdanx(i,1) = Kxxbar(i,1)*computekr_N(Sw(i,1))/viscN(i,1)
                    xiNbarx(i,1) = xiN(i,1)
                    densiNbarx(i,1) = densiN(i,1)
                    xNbarx(1:Nc,i,1) = xN(1:Nc,i,1)
                end if
            end if
        end do

        do i = 1, ny
            if((Uwx(i,nx+1) < 0).or.(Unx(i,nx+1) < 0)) then
                ztemp(1:Nc) = zBdryX(1:Nc,2,i)
#ifdef FULL_2c
                call flashcalculation_fullgrid_2c(Pw(i+1,nx+1), ztemp, xWtemp, xNtemp, xiWbarx(i,nx+1), xiNbarx(i,nx+1), &!
                    densiWbarx(i,nx+1), densiNbarx(i,nx+1), Swtemp, vtemp, ignore1, isW, isN)
#elif FULL_3c
                call flashcalculation_fullgrid_3c(Pw(i+1,nx+1), ztemp, xWtemp, xNtemp, xiWbarx(i,nx+1), xiNbarx(i,nx+1), &!
                    densiWbarx(i,nx+1), densiNbarx(i,nx+1), Swtemp, vtemp, ignore1, isW, isN)
#elif SPARSE
                call flashcalculation_sparsegrid(Pw(i+1,nx+1), ztemp, xWtemp, xNtemp, xiWbarx(i,nx+1), &!
                    xiNbarx(i,nx+1), densiWbarx(i,nx+1), densiNbarx(i,nx+1), Swtemp, vtemp, ignore1, isW, isN)
#elif NN
                call flashcalculation_nn(Pw(i+1,nx+1), ztemp, xWtemp, xNtemp, xiWbarx(i,nx+1), &!
                    xiNbarx(i,nx+1), densiWbarx(i,nx+1), densiNbarx(i,nx+1), Swtemp, vtemp, ignore1, isW, isN)
#else
                call flashcalculation(Pw(i+1,nx+1), ztemp, xWtemp, xNtemp, xiWbarx(i,nx+1), &!
                    xiNbarx(i,nx+1), densiWbarx(i,nx+1), densiNbarx(i,nx+1), Swtemp, vtemp, ignore1, isW, isN, isRea)
#endif

                xWbarx(1:Nc,i,nx+1) = xWtemp(1:Nc)
                xNbarx(1:Nc,i,nx+1) = xNtemp(1:Nc)
                if(isW) then
                    viscWtemp = viscosity( xWtemp, xiWbarx(i,nx+1), Pw(i+1,nx+1), 'l' )
                end if
                if(isN) then
                    viscNtemp = viscosity( xNtemp, xiNbarx(i,nx+1), Pw(i+1,nx+1), 'g' )
                end if
            end if
            if(Uwx(i,nx+1) < 0) then
                lambdawx(i,nx+1) = Kxxbar(i,nx+1)*computekr_W(Swtemp)/viscWtemp
            elseif(Uwx(i,nx+1) > 0) then
                lambdawx(i,nx+1) = Kxxbar(i,nx+1)*computekr_W(Sw(i,nx))/viscW(i,nx)
                xiWbarx(i,nx+1) = xiW(i,nx)
                densiWbarx(i,nx+1) = densiW(i,nx)
                xWbarx(1:Nc,i,nx+1) = xW(1:Nc,i,nx)
            else
                if(isDiriX(2,i) == 1) then
                    lambdawx(i,nx+1) = Kxxbar(i,nx+1)*computekr_W(Sw(i,nx))/viscW(i,nx)
                    xiWbarx(i,nx+1) = xiW(i,nx)
                    densiWbarx(i,nx+1) = densiW(i,nx)
                    xWbarx(1:Nc,i,nx+1) = xW(1:Nc,i,nx)
                end if
            end if
            if(Unx(i,nx+1) < 0) then
                lambdanx(i,nx+1) = Kxxbar(i,nx+1)*computekr_N(Swtemp)/viscNtemp
            elseif(Unx(i,nx+1) > 0) then
                lambdanx(i,nx+1) = Kxxbar(i,nx+1)*computekr_N(Sw(i,nx))/viscN(i,nx)
                xiNbarx(i,nx+1) = xiN(i,nx)
                densiNbarx(i,nx+1) = densiN(i,nx)
                xNbarx(1:Nc,i,nx+1) = xN(1:Nc,i,nx)
            else
                if(isDiriX(2,i) == 1) then
                    lambdanx(i,nx+1) = Kxxbar(i,nx+1)*computekr_N(Sw(i,nx))/viscN(i,nx)
                    xiNbarx(i,nx+1) = xiN(i,nx)
                    densiNbarx(i,nx+1) = densiN(i,nx)
                    xNbarx(1:Nc,i,nx+1) = xN(1:Nc,i,nx)
                end if
            end if
        end do

        do j = 2, nx
            do i = 1, ny
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
                if((Unx(i,j) > 0)) then
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

        do i = 1, nx
            if((Uwy(1,i) > 0).or.(Uny(1,i) > 0)) then
                ztemp(1:Nc) = zBdryY(1:Nc,i,1)
#ifdef FULL_2c
                call flashcalculation_fullgrid_2c(Pw(2,i+1), ztemp, xWtemp, xNtemp, xiWbary(1,i), xiNbary(1,i), &!
                    densiWbary(1,i), densiNbary(1,i), Swtemp, vtemp, ignore1, isW, isN)
#elif FULL_3c
                call flashcalculation_fullgrid_3c(Pw(2,i+1), ztemp, xWtemp, xNtemp, xiWbary(1,i), xiNbary(1,i), &!
                    densiWbary(1,i), densiNbary(1,i), Swtemp, vtemp, ignore1, isW, isN)
#elif SPARSE
                call flashcalculation_sparsegrid(Pw(2,i+1), ztemp, xWtemp, xNtemp, xiWbary(1,i), xiNbary(1,i), &!
                    densiWbary(1,i), densiNbary(1,i), Swtemp, vtemp, ignore1, isW, isN)
#elif NN
                call flashcalculation_nn(Pw(2,i+1), ztemp, xWtemp, xNtemp, xiWbary(1,i), xiNbary(1,i), &!
                    densiWbary(1,i), densiNbary(1,i), Swtemp, vtemp, ignore1, isW, isN)
#else
                call flashcalculation(Pw(2,i+1), ztemp, xWtemp, xNtemp, xiWbary(1,i), xiNbary(1,i), &!
                    densiWbary(1,i), densiNbary(1,i), Swtemp, vtemp, ignore1, isW, isN, isRea)
#endif

                xWbary(1:Nc,1,i) = xWtemp(1:Nc)
                xNbary(1:Nc,1,i) = xNtemp(1:Nc)
                if(isW) then
                    viscWtemp = viscosity( xWtemp, xiWbary(1,i), Pw(2,i+1), 'l' )
                end if
                if(isN) then
                    viscNtemp = viscosity( xNtemp, xiNbary(1,i), Pw(2,i+1), 'g' )
                end if
            end if
            if(Uwy(1,i) > 0) then
                lambdawy(1,i) = Kyybar(1,i)*computekr_W(Swtemp)/viscWtemp
            elseif(Uwy(1,i) < 0) then
                lambdawy(1,i) = Kyybar(1,i)*computekr_W(Sw(1,i))/viscW(1,i)
                xiWbary(1,i) = xiW(1,i)
                densiWbary(1,i) = densiW(1,i)
                xWbary(1:Nc,1,i) = xW(1:Nc,1,i)
            else
                if(isDiriY(i,1) == 1) then
                    lambdawy(1,i) = Kyybar(1,i)*computekr_W(Sw(1,i))/viscW(1,i)
                    xiWbary(1,i) = xiW(1,i)
                    densiWbary(1,i) = densiW(1,i)
                    xWbary(1:Nc,1,i) = xW(1:Nc,1,i)
                end if
            end if
            if(Uny(1,i) > 0) then
                lambdany(1,i) = Kyybar(1,i)*computekr_N(Swtemp)/viscNtemp
            elseif(Uny(1,i) < 0) then
                lambdany(1,i) = Kyybar(1,i)*computekr_N(Sw(1,i))/viscN(1,i)
                xiNbary(1,i) = xiN(1,i)
                densiNbary(1,i) = densiN(1,i)
                xNbary(1:Nc,1,i) = xN(1:Nc,1,i)
            else
                if(isDiriY(i,1) == 1) then
                    lambdany(1,i) = Kyybar(1,i)*computekr_N(Sw(1,i))/viscN(1,i)
                    xiNbary(1,i) = xiN(1,i)
                    densiNbary(1,i) = densiN(1,i)
                    xNbary(1:Nc,1,i) = xN(1:Nc,1,i)
                end if
            end if
        end do

        do i = 1, nx
            if((Uwy(ny+1,i) < 0).or.(Uny(ny+1,i) < 0)) then
                ztemp(1:Nc) = zBdryY(1:Nc,i,2)
#ifdef FULL_2c
                call flashcalculation_fullgrid_2c(Pw(ny+1,i+1), ztemp, xWtemp, xNtemp, xiWbary(ny+1,i), &!
                    xiNbary(ny+1,i), densiWbary(ny+1,i), densiNbary(ny+1,i), Swtemp, vtemp, ignore1, isW, isN)
#elif FULL_3c
                call flashcalculation_fullgrid_3c(Pw(ny+1,i+1), ztemp, xWtemp, xNtemp, xiWbary(ny+1,i), &!
                    xiNbary(ny+1,i), densiWbary(ny+1,i), densiNbary(ny+1,i), Swtemp, vtemp, ignore1, isW, isN)
#elif SPARSE
                call flashcalculation_sparsegrid(Pw(ny+1,i+1), ztemp, xWtemp, xNtemp, xiWbary(ny+1,i), &!
                    xiNbary(ny+1,i), densiWbary(ny+1,i), densiNbary(ny+1,i), Swtemp, vtemp, ignore1, isW, isN)
#elif NN
                call flashcalculation_nn(Pw(ny+1,i+1), ztemp, xWtemp, xNtemp, xiWbary(ny+1,i), &!
                    xiNbary(ny+1,i), densiWbary(ny+1,i), densiNbary(ny+1,i), Swtemp, vtemp, ignore1, isW, isN)
#else
                call flashcalculation(Pw(ny+1,i+1), ztemp, xWtemp, xNtemp, xiWbary(ny+1,i), &!
                    xiNbary(ny+1,i), densiWbary(ny+1,i), densiNbary(ny+1,i), Swtemp, vtemp, ignore1, isW, isN, isRea)
#endif

                xWbary(1:Nc,ny+1,i) = xWtemp(1:Nc)
                xNbary(1:Nc,ny+1,i) = xNtemp(1:Nc)
                if(isW) then
                    viscWtemp = viscosity( xWtemp, xiWbary(ny+1,i), Pw(ny+1,i+1), 'l' )
                end if
                if(isN) then
                    viscNtemp = viscosity( xNtemp, xiNbary(ny+1,i), Pw(ny+1,i+1), 'g' )
                end if
            end if
            if(Uwy(ny+1,i) < 0) then
                lambdawy(ny+1,i) = Kyybar(ny+1,i)*computekr_W(Swtemp)/viscWtemp
            elseif(Uwy(ny+1,i) > 0) then
                lambdawy(ny+1,i) = Kyybar(ny+1,i)*computekr_W(Sw(ny,i))/viscW(ny,i)
                xiWbary(ny+1,i) = xiW(ny,i)
                densiWbary(ny+1,i) = densiW(ny,i)
                xWbary(1:Nc,ny+1,i) = xW(1:Nc,ny,i)
            else
                if(isDiriY(i,2) == 1) then
                    lambdawy(ny+1,i) = Kyybar(ny+1,i)*computekr_W(Sw(ny,i))/viscW(ny,i)
                    xiWbary(ny+1,i) = xiW(ny,i)
                    densiWbary(ny+1,i) = densiW(ny,i)
                    xWbary(1:Nc,ny+1,i) = xW(1:Nc,ny,i)
                end if
            end if
            if(Uny(ny+1,i) < 0) then
                lambdany(ny+1,i) = Kyybar(ny+1,i)*computekr_N(Swtemp)/viscNtemp
            elseif(Uny(ny+1,i) > 0) then
                lambdany(ny+1,i) = Kyybar(ny+1,i)*computekr_N(Sw(ny,i))/viscN(ny,i)
                xiNbary(ny+1,i) = xiN(ny,i)
                densiNbary(ny+1,i) = densiN(ny,i)
                xNbary(1:Nc,ny+1,i) = xN(1:Nc,ny,i)
            else
                if(isDiriY(i,2) == 1) then
                    lambdany(ny+1,i) = Kyybar(ny+1,i)*computekr_N(Sw(ny,i))/viscN(ny,i)
                    xiNbary(ny+1,i) = xiN(ny,i)
                    densiNbary(ny+1,i) = densiN(ny,i)
                    xNbary(1:Nc,ny+1,i) = xN(1:Nc,ny,i)
                end if
            end if
        end do

        do j = 1, nx
            do i = 2, ny
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

        deallocate(ztemp)
        deallocate(xWtemp)
        deallocate(xNtemp)
        deallocate(vtemp)

    end subroutine computeParameters

    subroutine computePres()

        real(kind=8), dimension(:,:), pointer :: A
        real(kind=8), dimension(:), pointer :: b
        real(kind=8) :: xedge, yedge, ledge, redge, uedge, dedge
        real(kind=8) :: cotwx2, cotnx2, cotwx1, cotnx1, cotwy1, cotny1, cotwy2, cotny2
        real(kind=8) :: up, down, left, right
        real(kind=8) :: sum1, sum2
        integer :: i, j, m, r
        integer :: INFO
        integer, dimension(:), pointer :: IPIV

        allocate(A(nx*ny,nx*ny))
        allocate(b(nx*ny))
        allocate(IPIV(nx*ny))
        A(:,:) = 0.D0
        b(:) = 0.D0
        IPIV(:) = 0.D0

        do i = 2, ny+1
            do j = 2, nx+1
                r = (i-2)*nx + (j-1)

                xedge = xs(j) - xs(j-1)
                yedge = ys(i) - ys(i-1)

                if(j /= 2) then
                    ledge = xs(j-1) - xs(j-2)
                else
                    ledge = 0.D0
                end if

                if(j /= nx+1) then
                    redge = xs(j+1) - xs(j)
                else
                    redge = 0.D0
                end if

                if(i /= ny+1) then
                    uedge = ys(i+1) - ys(i)
                else
                    uedge = 0.D0
                end if

                if(i /= 2) then
                    dedge = ys(i-1) - ys(i-2)
                else
                    dedge = 0.D0
                end if

                A(r,r) = poro(j-1,i-1)*Cf(i-1,j-1)/(timeEnd/nt)
                b(r) = poro(j-1,i-1)*Cf(i-1,j-1)/(timeEnd/nt)*Pw(i,j)
                do m = 1, Nc
                    b(r) = b(r) + v(m,i-1,j-1)*src(m,j-1,i-1)
                end do

                cotwx2 = 0.D0
                cotnx2 = 0.D0
                cotwx1 = 0.D0
                cotnx1 = 0.D0
                cotwy1 = 0.D0
                cotny1 = 0.D0
                cotwy2 = 0.D0
                cotny2 = 0.D0
                do m = 1, Nc
                    cotwx2 = cotwx2 + lambdawx(i-1,j-1)*xiWbarx(i-1,j-1)*xWbarx(m,i-1,j-1) &!
                        *v(m,i-1,j-1)
                    cotnx2 = cotnx2 + lambdanx(i-1,j-1)*xiNbarx(i-1,j-1)*xNbarx(m,i-1,j-1) &!
                        *v(m,i-1,j-1)
                    cotwx1 = cotwx1 + lambdawx(i-1,j)*xiWbarx(i-1,j)*xWbarx(m,i-1,j)*v(m,i-1,j-1)
                    cotnx1 = cotnx1 + lambdanx(i-1,j)*xiNbarx(i-1,j)*xNbarx(m,i-1,j)*v(m,i-1,j-1)
                    cotwy1 = cotwy1 + lambdawy(i,j-1)*xiWbary(i,j-1)*xWbary(m,i,j-1)*v(m,i-1,j-1)
                    cotny1 = cotny1 + lambdany(i,j-1)*xiNbary(i,j-1)*xNbary(m,i,j-1)*v(m,i-1,j-1)
                    cotwy2 = cotwy2 + lambdawy(i-1,j-1)*xiWbary(i-1,j-1)*xWbary(m,i-1,j-1) &!
                        *v(m,i-1,j-1)
                    cotny2 = cotny2 + lambdany(i-1,j-1)*xiNbary(i-1,j-1)*xNbary(m,i-1,j-1) &!
                        *v(m,i-1,j-1)
                end do

                up = -2.D0*(cotwy1+cotny1)/yedge/(yedge+uedge)
                down = -2.D0*(cotwy2+cotny2)/yedge/(yedge+dedge)
                left = -2.D0*(cotwx2+cotnx2)/xedge/(xedge+ledge)
                right = -2.D0*(cotwx1+cotnx1)/xedge/(xedge+redge)

                if((i == ny+1).and.(isDiriY(j-1,2) == 0)) then
                    sum1 = 0.D0
                    sum2 = 0.D0
                    do m = 1, Nc
                        sum1 = sum1 + xWbary(m,i,j-1)*v(m,i-1,j-1)
                        sum2 = sum2 + xNbary(m,i,j-1)*v(m,i-1,j-1)
                    end do
                    b(r) = b(r) - UwBdryY(j-1,2)*xiWbary(i,j-1)*sum1/yedge - &!
                        UnBdryY(j-1,2)*xiNbary(i,j-1)*sum2/yedge
                else if((i == ny+1).and.(isDiriY(j-1,2) == 1)) then
                    b(r) = b(r) - up*Pw(i+1,j) - cotwy1*densiWbary(i,j-1)*gravY/yedge - &!
                        cotny1*densiNbary(i,j-1)*gravY/yedge
                    A(r,r) = A(r,r)-up
                else
                    A(r,r+nx) = up
                    A(r,r) = A(r,r)-up
                    b(r) = b(r) - cotwy1*densiWbary(i,j-1)*gravY/yedge - &!
                        cotny1*densiNbary(i,j-1)*gravY/yedge
                end if

                if((i == 2).and.(isDiriY(j-1,1) == 0)) then
                    sum1 = 0.D0
                    sum2 = 0.D0
                    do m = 1, Nc
                        sum1 = sum1 + xWbary(m,i-1,j-1)*v(m,i-1,j-1)
                        sum2 = sum2 + xNbary(m,i-1,j-1)*v(m,i-1,j-1)
                    end do
                    b(r) = b(r) + UwBdryY(j-1,1)*xiWbary(i-1,j-1)*sum1/yedge + &!
                        UnBdryY(j-1,1)*xiNbary(i-1,j-1)*sum2/yedge
                else if((i == 2).and.(isDiriY(j-1,1) == 1)) then
                    b(r) = b(r) - down*Pw(i-1,j) + cotwy2*densiWbary(i-1,j-1)*gravY/yedge + &!
                        cotny2*densiNbary(i-1,j-1)*gravY/yedge
                    A(r,r) = A(r,r)-down
                else
                    A(r,r-nx) = down
                    A(r,r) = A(r,r)-down
                    b(r) = b(r) + cotwy2*densiWbary(i-1,j-1)*gravY/yedge + &!
                        cotny2*densiNbary(i-1,j-1)*gravY/yedge
                end if

                if((j == 2).and.(isDiriX(1,i-1) == 0)) then
                    sum1 = 0.D0
                    sum2 = 0.D0
                    do m = 1, Nc
                        sum1 = sum1 + xWbarx(m,i-1,j-1)*v(m,i-1,j-1)
                        sum2 = sum2 + xNbarx(m,i-1,j-1)*v(m,i-1,j-1)
                    end do
                    b(r) = b(r) + UwBdryX(1,i-1)*xiWbarx(i-1,j-1)*sum1/xedge + &!
                        UnBdryX(1,i-1)*xiNbarx(i-1,j-1)*sum2/xedge
                else if((j == 2).and.(isDiriX(1,i-1) == 1)) then
                    b(r) = b(r) - left*Pw(i,j-1) + cotwx2*densiWbarx(i-1,j-1)*gravX/xedge + &!
                        cotnx2*densiNbarx(i-1,j-1)*gravX/xedge
                    A(r,r) = A(r,r)-left
                else
                    A(r,r-1) = left
                    A(r,r) = A(r,r)-left
                    b(r) = b(r) + cotwx2*densiWbarx(i-1,j-1)*gravX/xedge + &!
                        cotnx2*densiNbarx(i-1,j-1)*gravX/xedge
                end if

                if((j == nx+1).and.(isDiriX(2,i-1) == 0)) then
                    sum1 = 0.D0
                    sum2 = 0.D0
                    do m = 1, Nc
                        sum1 = sum1 + xWbarx(m,i-1,j)*v(m,i-1,j-1)
                        sum2 = sum2 + xNbarx(m,i-1,j)*v(m,i-1,j-1)
                    end do
                    b(r) = b(r) - UwBdryX(2,i-1)*xiWbarx(i-1,j)*sum1/xedge - &!
                        UnBdryX(2,i-1)*xiNbarx(i-1,j)*sum2/xedge
                else if((j == nx+1).and.(isDiriX(2,i-1) == 1)) then
                    b(r) = b(r) - right*Pw(i,j+1) - cotwx1*densiWbarx(i-1,j)*gravX/xedge - &!
                        cotnx1*densiNbarx(i-1,j)*gravX/xedge
                    A(r,r) = A(r,r)-right
                else
                    A(r,r+1) = right
                    A(r,r) = A(r,r)-right
                    b(r) = b(r) - cotwx1*densiWbarx(i-1,j)*gravX/xedge - &!
                        cotnx1*densiNbarx(i-1,j)*gravX/xedge
                end if
            end do
        end do

        call dgesv(nx*ny, 1, A, nx*ny, IPIV, b, nx*ny, INFO)

        do i = 2, ny+1
            do j = 2, nx+1
                Pw(i,j) = b((i-2)*nx+(j-1))
            end do
        end do

        deallocate(A)
        deallocate(b)
        deallocate(IPIV)

    end subroutine computePres

    subroutine computeMoleFrac()

        real(kind=8), dimension(:), pointer :: xic
        real(kind=8) :: right, left, up, down, div
        real(kind=8) :: xinew
        real(kind=8) :: timefinish
        integer :: i, j, m

        allocate(xic(Nc))

        do j = 1, nx
            do i = 1, ny
                do m = 1, Nc
                    right = Uwx(i,j+1)*xiWbarx(i,j+1)*xWbarx(m,i,j+1) + Unx(i,j+1)* &!
                        xiNbarx(i,j+1)*xNbarx(m,i,j+1)
                    left = Uwx(i,j)*xiWbarx(i,j)*xWbarx(m,i,j) + &!
                        Unx(i,j)*xiNbarx(i,j)*xNbarx(m,i,j)
                    up = Uwy(i+1,j)*xiWbary(i+1,j)*xWbary(m,i+1,j) + Uny(i+1,j)* &!
                        xiNbary(i+1,j)*xNbary(m,i+1,j)
                    down = Uwy(i,j)*xiWbary(i,j)*xWbary(m,i,j) + &!
                        Uny(i,j)*xiNbary(i,j)*xNbary(m,i,j)
                    div = (right-left)/(xs(j+1)-xs(j)) + (up-down)/(ys(i+1)-ys(i))
                    xic(m) = (src(m,j,i)-div)*(timeEnd/nt)/poro(j,i) + z(m,i,j)*(Sw(i,j)*xiW(i,j)+(1-Sw(i,j))*xiN(i,j))
                    if(xic(m) < 0) then
                        print *, 'Please tune the time step.', xic(m), i, j, m, Uwx(i,j+1), &!
                            xiWbarx(i,j+1), xWbarx(m,i,j+1),Unx(i,j+1), xiNbarx(i,j+1),xNbarx(m,i,j+1)
                        stop
                    end if
                end do
                xinew = 0.D0
                do m = 1, Nc
                    xinew = xinew + xic(m)
                end do
                z(1:Nc,i,j) = xic(1:Nc)/xinew
                do m = 1, Nc
                    if(z(m,i,j) < 1.D-12) then
                        z(m,i,j) = 0.D0
                    end if
                end do
            end do
        end do

        call CPU_TIME(timefinish)
        write(90, fmt="(es15.8)") timefinish-timebegin

        deallocate(xic)

    end subroutine computeMoleFrac

    subroutine computeVel()

        integer :: i, j

        ! compute the total new velocities in x direction
        do i = 1, ny
            if(isDiriX(1,i) == 1) then
                Uwx(i,1) = -lambdawx(i,1)*((Pw(i+1,2)-Pw(i+1,1))*2.D0/(xs(2)-xs(1)) - densiWbarx(i,1)*gravX)
                Unx(i,1) = -lambdanx(i,1)*((Pw(i+1,2)-Pw(i+1,1))*2.D0/(xs(2)-xs(1)) - densiNbarx(i,1)*gravX)
            else
                Uwx(i,1) = UwBdryX(1,i)
                Unx(i,1) = UnBdryX(1,i)
            end if
        end do

        do i = 1, ny
            if(isDiriX(2,i) == 1) then
                Uwx(i,nx+1) = -lambdawx(i,nx+1)*((Pw(i+1,nx+2)-Pw(i+1,nx+1))*2.D0/(xs(nx+1)-xs(nx)) - densiWbarx(i,nx+1)*gravX)
                Unx(i,nx+1) = -lambdanx(i,nx+1)*((Pw(i+1,nx+2)-Pw(i+1,nx+1))*2.D0/(xs(nx+1)-xs(nx)) - densiNbarx(i,nx+1)*gravX)
                if(Uwx(i,nx+1) < 0) then
                    Uwx(i,nx+1) = 0.D0
                end if
                if(Unx(i,nx+1) < 0) then
                    Unx(i,nx+1) = 0.D0
                end if
            else
                Uwx(i,nx+1) = UwBdryX(2,i)
                Unx(i,nx+1) = UnBdryX(2,i)
            end if
        end do

        do j = 2, nx
            do i = 1, ny
                Uwx(i,j) = -lambdawx(i,j)*((Pw(i+1,j+1)-Pw(i+1,j))*2.D0/(xs(j+1)-xs(j-1)) - densiWbarx(i,j)*gravX)
                Unx(i,j) = -lambdanx(i,j)*((Pw(i+1,j+1)-Pw(i+1,j))*2.D0/(xs(j+1)-xs(j-1)) - densiNbarx(i,j)*gravX)
            end do
        end do

        ! compute the total new velocities in y direction
        do j = 1, nx
            if(isDiriY(j,1) == 1) then
                Uwy(1,j) = -lambdawy(1,j)*((Pw(2,j+1)-Pw(1,j+1))*2.D0/(ys(2)-ys(1)) - densiWbary(1,j)*gravY)
                Uny(1,j) = -lambdany(1,j)*((Pw(2,j+1)-Pw(1,j+1))*2.D0/(ys(2)-ys(1)) - densiNbary(1,j)*gravY)
                if(Uwy(1,j) > 0) then
                    Uwy(1,j) = 0.D0
                end if
                if(Uny(1,j) > 0) then
                    Uny(1,j) = 0.D0
                end if
            else
                Uwy(1,j) = UwBdryY(j,1)
                Uny(1,j) = UnBdryY(j,1)
            end if
        end do

        do j = 1, nx
            if(isDiriY(j,2) == 1) then
                Uwy(ny+1,j) = -lambdawy(ny+1,j)*((Pw(ny+2,j+1)-Pw(ny+1,j+1))*2.D0/(ys(ny+1)-ys(ny)) - densiWbary(ny+1,j)*gravY)
                Uny(ny+1,j) = -lambdany(ny+1,j)*((Pw(ny+2,j+1)-Pw(ny+1,j+1))*2.D0/(ys(ny+1)-ys(ny)) - densiNbary(ny+1,j)*gravY)
            else
                Uwy(ny+1,j) = UwBdryY(j,2)
                Uny(ny+1,j) = UnBdryY(j,2)
            end if
        end do

        do j = 1, nx
            do i = 2, ny
                Uwy(i,j) = -lambdawy(i,j)*((Pw(i+1,j+1)-Pw(i,j+1))*2.D0/(ys(i+1)-ys(i-1)) - densiWbary(i,j)*gravY)
                Uny(i,j) = -lambdany(i,j)*((Pw(i+1,j+1)-Pw(i,j+1))*2.D0/(ys(i+1)-ys(i-1)) - densiNbary(i,j)*gravY)
            end do
        end do

    end subroutine computeVel

    subroutine finalize()

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

#if defined(FULL_2c) || defined(FULL_3c)
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

        close(70)
        close(80)
        close(90)

#if defined(FULL_2c) || defined(FULL_3c)
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
        if(Nc == 3) then
            close(240)
            close(250)
            close(260)
        end if
        close(270)
        close(280)
#endif

    end subroutine finalize

end module RST_compositionalTwoPhaseFlow
