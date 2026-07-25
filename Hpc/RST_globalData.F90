
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_globalData

    implicit none

    ! the program parameters
    integer, parameter :: MAX_BUF = 1.D5
    integer, parameter :: MAX_COMMBUF = 1.D6

    ! the model parameters
    integer :: pncols 
    integer :: pnrows
    real(kind=8) :: Lx 
    real(kind=8) :: Ly 
    real(kind=8) :: timeEnd 
    integer :: nx 
    integer :: ny 
    integer :: nt
    real(kind=8) :: gravX
    real(kind=8) :: gravY
    real(kind=8), dimension(:), pointer :: xs
    real(kind=8), dimension(:), pointer :: ys 
    real(kind=8), dimension(:), pointer :: ts 
    real(kind=8), dimension(:,:), pointer :: Kxx 
    real(kind=8), dimension(:,:), pointer :: Kyy 
    real(kind=8), dimension(:,:), pointer :: poro 
    real(kind=8), dimension(:,:,:), pointer :: src 
    integer, dimension(:,:), pointer :: isDiriX 
    integer, dimension(:,:), pointer :: isDiriY 
    real(kind=8), dimension(:,:), pointer :: PwBdryX 
    real(kind=8), dimension(:,:), pointer :: PwBdryY 
    real(kind=8), dimension(:,:), pointer :: PwInit 
    real(kind=8), dimension(:,:,:), pointer :: zBdryX 
    real(kind=8), dimension(:,:,:), pointer :: zBdryY 
    real(kind=8), dimension(:,:,:), pointer :: zInit 
    real(kind=8), dimension(:,:), pointer :: UwBdryX 
    real(kind=8), dimension(:,:), pointer :: UwBdryY 
    real(kind=8), dimension(:,:), pointer :: UnBdryX 
    real(kind=8), dimension(:,:), pointer :: UnBdryY 
    character(len = 10) :: soludoc 

    ! the global variables
    real(kind=8), dimension(:,:), pointer :: Pw
    real(kind=8), dimension(:,:), pointer :: Uwx
    real(kind=8), dimension(:,:), pointer :: Uwy
    real(kind=8), dimension(:,:), pointer :: Unx
    real(kind=8), dimension(:,:), pointer :: Uny
    real(kind=8), dimension(:,:), pointer :: Sw
    real(kind=8), dimension(:,:), pointer :: lambdawx
    real(kind=8), dimension(:,:), pointer :: lambdawy
    real(kind=8), dimension(:,:), pointer :: lambdanx
    real(kind=8), dimension(:,:), pointer :: lambdany
    real(kind=8), dimension(:,:), pointer :: Kxxbar
    real(kind=8), dimension(:,:), pointer :: Kyybar
    real(kind=8), dimension(:,:,:), pointer :: z
    real(kind=8), dimension(:,:), pointer :: densiW
    real(kind=8), dimension(:,:), pointer :: densiN
    real(kind=8), dimension(:,:), pointer :: densiWbarx
    real(kind=8), dimension(:,:), pointer :: densiWbary
    real(kind=8), dimension(:,:), pointer :: densiNbarx
    real(kind=8), dimension(:,:), pointer :: densiNbary
    real(kind=8), dimension(:,:,:), pointer :: xW
    real(kind=8), dimension(:,:,:), pointer :: xN
    real(kind=8), dimension(:,:,:), pointer :: xWbarx
    real(kind=8), dimension(:,:,:), pointer :: xWbary
    real(kind=8), dimension(:,:,:), pointer :: xNbarx
    real(kind=8), dimension(:,:,:), pointer :: xNbary
    real(kind=8), dimension(:,:), pointer :: xiW
    real(kind=8), dimension(:,:), pointer :: xiN
    real(kind=8), dimension(:,:), pointer :: xiWbarx
    real(kind=8), dimension(:,:), pointer :: xiWbary
    real(kind=8), dimension(:,:), pointer :: xiNbarx
    real(kind=8), dimension(:,:), pointer :: xiNbary
    real(kind=8), dimension(:,:,:), pointer :: v
    real(kind=8), dimension(:,:), pointer :: Cf
    real(kind=8), dimension(:,:), pointer :: viscW
    real(kind=8), dimension(:,:), pointer :: viscN
    real(kind=8) :: totalmole
    integer :: t

    ! the HPC variables
    real(kind=8) :: timestart, flashtime, solvertime, commtime
    integer :: num_procs, myid
    integer :: localncols, localnrows
    integer :: pcol, prow
    integer :: xlower, xupper, ylower, yupper
    integer :: local_size
    real(kind=8), dimension(:), pointer :: initial_x_guess
    integer :: stencil_indices(5)
    integer :: offsets(5,2)
    integer :: ilower(2), iupper(2)
    integer(kind=8) :: grid
    integer(kind=8) :: stencil
    integer(kind=8) :: global_A
    integer(kind=8) :: global_b
    integer(kind=8) :: global_x
    integer(kind=8) :: solver
    
end module RST_globalData

