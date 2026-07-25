
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_compositionalTwoPhaseFlowDriver

    use RST_model
    use RST_globalData
    use RST_compositionalTwoPhaseFlow
    use RST_writeFile
    use RST_genPlotfileMatlab

    implicit none

contains

    subroutine driver(modelCase)

        type(model), intent(in out) :: modelCase

        call initialize(modelCase)

        do t = 2, nt+1

            if(mod(t,100) == 0) then
                print *, t
            end if

            call computeParameters()

            call computePres()

            call computeVel()

            call computeMoleFrac()
            
        end do

        call writeFile()

        call genPlotfileMatlab()

        call finalize()

    end subroutine driver

end module RST_compositionalTwoPhaseFlowDriver
