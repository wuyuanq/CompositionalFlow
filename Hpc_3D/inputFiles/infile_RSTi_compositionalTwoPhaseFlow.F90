
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com
!!$
!!$ This is a two-phase case. Inject to one corner above and output from one corner below.
!!$ Stipulate that the injected component is always numbered as 1, and the components inside
!!$ the reservoir are numbered 2 and so on.

program infile_RSTi_compositionalTwoPhaseFlow

    use RST_model
    use RST_compositionalTwoPhaseFlowDriver
    use RST_proceAlloc
    implicit none

    type(model) :: modelCase
    integer :: i, j, k
    real(kind=8), parameter :: UCONST = 2.D-6
    real(kind=8), parameter :: PCONST = 2012018
    integer :: times !!!

    times = 1!!!
    modelCase%Nc = 3
    modelCase%Temp = 220.0
    modelCase%Lx = 0.2
    modelCase%Ly = 0.1
    modelCase%Lz = 0.4
    modelCase%timeEnd = 0.1*365*24*3600.0
    modelCase%nx = 8*times!!!!
    modelCase%ny = 4*times!!!!
    modelCase%nz = 16*times!!!!
    modelCase%nt = 0.1*365*24*3600

    call proceAlloc(1, modelCase%nx, modelCase%ny, modelCase%nz, modelCase%pncols, modelCase%pnrows, modelCase%pnlays)

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

    ! the input point is at the left, up and back corner
    allocate(modelCase%src(modelCase%Nc,modelCase%nx,modelCase%ny,modelCase%nz))
    modelCase%src(:,:,:,:) = 0.0
    modelCase%src(1,1,modelCase%ny,modelCase%nz) = 0.02

    allocate(modelCase%isDiriX(2,modelCase%ny,modelCase%nz))
    modelCase%isDiriX(:,:,:) = 0
    modelCase%isDiriX(2,1:times,1:times) = 1!!!!!
    allocate(modelCase%isDiriY(modelCase%nx,2,modelCase%nz))
    modelCase%isDiriY(:,:,:) = 0
    modelCase%isDiriY(modelCase%nx-times+1:modelCase%nx,1,1:times) = 1!!!!!
    allocate(modelCase%isDiriZ(modelCase%nx,modelCase%ny,2))
    modelCase%isDiriZ(:,:,:) = 0
    modelCase%isDiriZ(modelCase%nx-times+1:modelCase%nx,1:times,1) = 1!!!!!

    ! the output point is at the right, down and front corner
    allocate(modelCase%PwBdryX(2,modelCase%ny,modelCase%nz))
    modelCase%PwBdryX(:,:,:) = 0.0
    modelCase%PwBdryX(2,1:times,1:times) = PCONST
    allocate(modelCase%PwBdryY(modelCase%nx,2,modelCase%nz))
    modelCase%PwBdryY(:,:,:) = 0.0
    modelCase%PwBdryY(modelCase%nx-times+1:modelCase%nx,1,1:times) = PCONST
    allocate(modelCase%PwBdryZ(modelCase%nx,modelCase%ny,2))
    modelCase%PwBdryZ(:,:,:) = 0.0
    modelCase%PwBdryZ(modelCase%nx-times+1:modelCase%nx,1:times,1) = PCONST

    allocate(modelCase%PwInit(modelCase%nx,modelCase%ny,modelCase%nz))
    do i = 1, modelCase%nx
        do j = 1, modelCase%ny
            do k = 1, modelCase%nz
                modelCase%PwInit(i, j, k) = PCONST
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
    modelCase%zInit(2,:,:,:) = 0.5
    modelCase%zInit(3,:,:,:) = 0.5

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

    ! 1 is CO2, 2 is CH4 and 3 is C3H8
    allocate(modelCase%ct(modelCase%Nc))
    modelCase%ct(1) = 304.4
    modelCase%ct(2) = 190
    modelCase%ct(3) = 370
    allocate(modelCase%cp(modelCase%Nc))
    modelCase%cp(1) = 7.4*1.D6
    modelCase%cp(2) = 4.6*1.D6
    modelCase%cp(3) = 4.2*1.D6
    allocate(modelCase%af(modelCase%Nc))
    modelCase%af(1) = 0.23
    modelCase%af(2) = 0.01
    modelCase%af(3) = 0.15
    allocate(modelCase%mw(modelCase%Nc))
    modelCase%mw(1) = 0.044
    modelCase%mw(2) = 0.016
    modelCase%mw(3) = 0.044
    allocate(modelCase%cv(modelCase%Nc))
    modelCase%cv(1) = 0.0021
    modelCase%cv(2) = 0.0062
    modelCase%cv(3) = 0.0045
    allocate(modelCase%psatA(modelCase%Nc))
    modelCase%psatA(1) = 6.81228
    modelCase%psatA(2) = 6.69561
    modelCase%psatA(3) = 6.82973
    allocate(modelCase%psatB(modelCase%Nc))
    modelCase%psatB(1) = 1301.679
    modelCase%psatB(2) = 405.420
    modelCase%psatB(3) = 813.2
    allocate(modelCase%psatC(modelCase%Nc))
    modelCase%psatC(1) = 269.506
    modelCase%psatC(2) = 267.777
    modelCase%psatC(3) = 248
    allocate(modelCase%delta(modelCase%Nc,modelCase%Nc))
    modelCase%delta(:,:) = 0.0
    modelCase%delta(1,2) = 0.15
    modelCase%delta(2,1) = modelCase%delta(1,2)
    modelCase%delta(1,3) = 0.1239
    modelCase%delta(3,1) = modelCase%delta(1,3)
    modelCase%delta(2,3) = 0.036
    modelCase%delta(3,2) = modelCase%delta(2,3)

    modelCase%soludoc = 'case5'

    call driver(modelCase)

end program infile_RSTi_compositionalTwoPhaseFlow
