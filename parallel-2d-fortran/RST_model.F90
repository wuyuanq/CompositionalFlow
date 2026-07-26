
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_model

    implicit none

    type :: model
        integer :: pncols ! the number of columns of processes, for the parallel program
        integer :: pnrows ! the number of rows of processes, for the parallel program
        integer :: Nc ! the number of components
        real(kind=8) :: Temp ! temperature
        real(kind=8) :: Lx ! the length of the reservoir in x-direction
        real(kind=8) :: Ly ! the length of the reservoir in y-direction
        real(kind=8) :: timeEnd ! the total simulation time
        integer :: nx ! the number of cells in x-direction
        integer :: ny ! the number of cells in y-direction
        integer :: nt ! the number of time steps
        real(kind=8) :: gravX ! the gravity in x-direction
        real(kind=8) :: gravY ! the gravity in y-direction
        real(kind=8), dimension(:), pointer :: xs ! the grid points in x-direction
        real(kind=8), dimension(:), pointer :: ys ! the grid points in y-direction
        real(kind=8), dimension(:), pointer :: ts ! the time points
        real(kind=8), dimension(:,:), pointer :: Kxx ! the permeability in x-direction
        real(kind=8), dimension(:,:), pointer :: Kyy ! the permeability in y-direction
        real(kind=8), dimension(:,:), pointer :: poro ! the porosity of the medium
        real(kind=8), dimension(:,:,:), pointer :: src ! the source
        integer, dimension(:,:), pointer :: isDiriX ! whether the boundary perpendicular with x-direction is Dirichlet condition, 1 means true, 0 means false
        integer, dimension(:,:), pointer :: isDiriY ! whether the boundary perpendicular with y-direction is Dirichlet condition, 1 means true, 0 means false
        real(kind=8), dimension(:,:), pointer :: PwBdryX ! the pressure on the boundary perpendicular with x-direction
        real(kind=8), dimension(:,:), pointer :: PwBdryY ! the pressure on the boundary perpendicular with y-direction
        real(kind=8), dimension(:,:), pointer :: PwInit ! the initial pressure field
        real(kind=8), dimension(:,:,:), pointer :: zBdryX ! the mole fraction on the boundary perpendicular with x-direction
        real(kind=8), dimension(:,:,:), pointer :: zBdryY ! the mole fraction on the boundary perpendicular with y-direction
        real(kind=8), dimension(:,:,:), pointer :: zInit ! the initial mole fraction
        real(kind=8), dimension(:,:), pointer :: UwBdryX ! the wetting velocities on the boundary perpendicular with x-direction
        real(kind=8), dimension(:,:), pointer :: UwBdryY ! the wetting velocities on the boundary perpendicular with y-direction
        real(kind=8), dimension(:,:), pointer :: UnBdryX ! the nonwetting velocities on the boundary perpendicular with x-direction
        real(kind=8), dimension(:,:), pointer :: UnBdryY ! the nonwetting velocities on the boundary perpendicular with y-direction
        real(kind=8), dimension(:), pointer :: ct ! critical temperature
        real(kind=8), dimension(:), pointer :: cp ! critical pressure
        real(kind=8), dimension(:), pointer :: af ! acentric factor
        real(kind=8), dimension(:), pointer :: mw ! mole weight
        real(kind=8), dimension(:), pointer :: cv ! critical volumn
        real(kind=8), dimension(:), pointer :: psatA ! the coefficients of saturation pressure equation
        real(kind=8), dimension(:), pointer :: psatB
        real(kind=8), dimension(:), pointer :: psatC
        real(kind=8), dimension(:,:), pointer :: delta ! binary interaction coefficient
        character(len = 10) :: soludoc ! the result document
        real(kind=8), dimension(:,:), pointer :: range ! for the Bitstringtree program
    end type model

contains

    function computeKr_W(satW) result(kr_W)

        implicit none
        real(kind=8), intent(in) :: satW
        real(kind=8) :: m4sat = 2.D0
        real(kind=8) :: kr_W

        kr_W = min(satW, 1.D0)**m4sat

    end function computeKr_W

    function computeKr_N(satW) result(kr_N)

        implicit none
        real(kind=8), intent(in) :: satW
        real(kind=8) :: m4sat = 2.D0
        real(kind=8) :: kr_N

        kr_N = (1-min(satW, 1.D0))**m4sat

    end function computeKr_N

end module RST_model
