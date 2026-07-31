
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com
!!$
!!$ This is a two-phase case. Inject to 1 and output from 2.

program infile_RSTi_compositionalTwoPhaseFlow

    use RST_model
    use RST_compositionalTwoPhaseFlowDriver
    use RST_proceAlloc

    implicit none

    type(model) :: modelCase
    integer :: i, j, k
    real(kind=8), parameter :: UCONST = 2.D-6
    real(kind=8), parameter :: PCONST = 2012018
    real(kind=8), dimension(:), pointer :: initPres
    integer :: times !!!!!

    times = 1!!!!!
    modelCase%Nc = 2
    modelCase%Temp = 480.0
    modelCase%Lx = 4.0
    modelCase%Ly = 4.0
    modelCase%Lz = 4.0
    modelCase%timeEnd = 0.1*365*24*3600.0
    modelCase%nx = 40*times!!!!
    modelCase%ny = 40*times!!!!
    modelCase%nz = 40*times!!!!
    modelCase%nt = 0.1*365*24*1!!!!!

    call call proceAlloc(8, modelCase%nx, modelCase%ny, modelCase%nz, modelCase%pncols, modelCase%pnrows, modelCase%pnlays)

    allocate(modelCase%xs(modelCase%nx+1))
    do i = 1, modelCase%nx+1
        modelCase%xs(i) = (i-1)*modelCase%Lx/modelCase%nx
    end do
    allocate(modelCase%ys(modelCase%ny+1))
    do i = 1, modelCase%ny+1
        modelCase%ys(i) = (i-1)*modelCase%Ly/modelCase%ny
    end do
    allocate(modelCase%zs(modelCase%nz+1))
    do i = 1, modelCase%nz+1
        modelCase%zs(i) = (i-1)*modelCase%Lz/modelCase%nz
    end do
    allocate(modelCase%ts(modelCase%nt+1))
    do i = 1, modelCase%nt+1
        modelCase%ts(i) = (i-1)*modelCase%timeEnd/modelCase%nt
    end do

    modelCase%gravX = 0.0
    modelCase%gravY = -9.807
    modelCase%gravZ = 0.0

    allocate(modelCase%Kxx(modelCase%nx,modelCase%ny,modelCase%nz))
    modelCase%Kxx(:,:,:) = 9.869233*1.D-15
    allocate(modelCase%Kyy(modelCase%nx,modelCase%ny,modelCase%nz))
    modelCase%Kyy = modelCase%Kxx
    allocate(modelCase%Kzz(modelCase%nx,modelCase%ny,modelCase%nz))
    modelCase%Kzz = modelCase%Kxx

    allocate(modelCase%poro(modelCase%nx,modelCase%ny,modelCase%nz))
    modelCase%poro(:,:,:) = 0.2

    allocate(modelCase%src(modelCase%Nc,modelCase%nx,modelCase%ny,modelCase%nz))
    modelCase%src(:,:,:,:) = 0.0
    modelCase%src(1,1,1,1) = 0.02

    allocate(modelCase%isDiriX(2,modelCase%ny,modelCase%nz))
    modelCase%isDiriX(:,:,:) = 0
    modelCase%isDiriX(2,1:times,1:times) = 1!!!!!
    allocate(modelCase%isDiriY(modelCase%nx,2,modelCase%nz))
    modelCase%isDiriY(:,:,:) = 0
    modelCase%isDiriY(modelCase%nx-times+1:modelCase%nx,1,1:times) = 1!!!!!
    allocate(modelCase%isDiriZ(modelCase%nx,modelCase%ny,2))
    modelCase%isDiriZ(:,:,:) = 0
    modelCase%isDiriZ(modelCase%nx-times+1:modelCase%nx,1:times,1) = 1!!!!!

    allocate(initPres(modelCase%ny))
    open(90,file='../InitialReservoirData/P4.txt')
    do i = 1, modelCase%ny
        read(90,*) initPres(i)
    end do

    allocate(modelCase%PwBdryX(2,modelCase%ny,modelCase%nz))
    modelCase%PwBdryX(:,:,:) = 0.0
    do j = 1, times
        do k = 1, times
            modelCase%PwBdryX(2,j,k) = initPres(j)
        end do
    end do
    allocate(modelCase%PwBdryY(modelCase%nx,2,modelCase%nz))
    modelCase%PwBdryY(:,:,:) = 0.0
    modelCase%PwBdryY(modelCase%nx-times+1:modelCase%nx,1,1:times) = initPres(1)
    allocate(modelCase%PwBdryZ(modelCase%nx,modelCase%ny,2))
    modelCase%PwBdryZ(:,:,:) = 0.0
    do i = modelCase%nx-times+1, modelCase%nx
        do j = 1, times
            modelCase%PwBdryZ(i,j,1) = initPres(j)
        end do
    end do

    allocate(modelCase%PwInit(modelCase%nx,modelCase%ny,modelCase%nz))
    do i = 1, modelCase%nx
        do j = 1, modelCase%ny
            do k = 1, modelCase%nz
                modelCase%PwInit(i, j, k) = initPres(j)
            end do
        end do
    end do

    allocate(modelCase%zBdryX(modelCase%Nc,2,modelCase%ny,modelCase%nz))
    modelCase%zBdryX(:,:,:,:) = 0.0
    allocate(modelCase%zBdryY(modelCase%Nc,modelCase%nx,2,modelCase%nz))
    modelCase%zBdryY(:,:,:,:) = 0.0
    allocate(modelCase%zBdryZ(modelCase%Nc,modelCase%nx,modelCase%ny,2))
    modelCase%zBdryZ(:,:,:,:) = 0.0

    allocate(modelCase%zInit(modelCase%Nc,modelCase%nx,modelCase%ny,modelCase%nz))
    modelCase%zInit(1,:,:,:) = 0.0
    modelCase%zInit(2,:,:,:) = 1.0

    allocate(modelCase%UwBdryX(2,modelCase%ny,modelCase%nz))
    modelCase%UwBdryX(:,:,:) = 0.0
    allocate(modelCase%UwBdryY(modelCase%nx,2,modelCase%nz))
    modelCase%UwBdryY(:,:,:) = 0.0
    allocate(modelCase%UwBdryZ(modelCase%nx,modelCase%ny,2))
    modelCase%UwBdryZ(:,:,:) = 0.0
    allocate(modelCase%UnBdryX(2,modelCase%ny,modelCase%nz))
    modelCase%UnBdryX(:,:,:) = 0.0
    allocate(modelCase%UnBdryY(modelCase%nx,2,modelCase%nz))
    modelCase%UnBdryY(:,:,:) = 0.0
    allocate(modelCase%UnBdryZ(modelCase%nx,modelCase%ny,2))
    modelCase%UnBdryZ(:,:,:) = 0.0

    allocate(modelCase%ct(modelCase%Nc))
    modelCase%ct(1) = 190
    modelCase%ct(2) = 370
    allocate(modelCase%cp(modelCase%Nc))
    modelCase%cp(1) = 4.6*1.D6
    modelCase%cp(2) = 4.2*1.D6
    allocate(modelCase%af(modelCase%Nc))
    modelCase%af(1) = 0.01
    modelCase%af(2) = 0.15
    allocate(modelCase%mw(modelCase%Nc))
    modelCase%mw(1) = 0.016
    modelCase%mw(2) = 0.044
    allocate(modelCase%cv(modelCase%Nc))
    modelCase%cv(1) = 0.0062
    modelCase%cv(2) = 0.0045
    allocate(modelCase%psatA(modelCase%Nc))
    modelCase%psatA(1) = 6.69561
    modelCase%psatA(2) = 6.82973
    allocate(modelCase%psatB(modelCase%Nc))
    modelCase%psatB(1) = 405.420
    modelCase%psatB(2) = 813.2
    allocate(modelCase%psatC(modelCase%Nc))
    modelCase%psatC(1) = 267.777
    modelCase%psatC(2) = 248
    allocate(modelCase%delta(modelCase%Nc,modelCase%Nc))
    modelCase%delta(:,:) = 0.0
    modelCase%delta(1,2) = 0.036
    modelCase%delta(2,1) = modelCase%delta(1,2)

    modelCase%soludoc = 'case2'

    call driver(modelCase)

    deallocate(initPres)

end program infile_RSTi_compositionalTwoPhaseFlow
