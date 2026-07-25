
# Author:
#   Yuanqing Wu, KAUST, Saudi Arabia
#
# History:
#   2015-2-9 by Yuanqing Wu
#
# Support:
#   wuyuanq@gmail.com

PROG = CompositionalFlow_fortran
defaults: ${PROG}

FCC = gfortran
#FCC = mpif90

#FCFLAGS = -O3 -fdefault-real-8 -g
FCFLAGS = -g -O3
#FCFLAGS = -O3 -fbounds-check 
#FCFLAGS = -O3

FLASHTYPE = FULL_2c
FLASHTYPE = FULL_3c
FLASHTYPE = SPARSE
FLASHTYPE = NN
#FLASHTYPE =

ifeq ($(FLASHTYPE), FULL_2c)
FLASHFLAGS = -DFULL_2c
else ifeq ($(FLASHTYPE), FULL_3c)
FLASHFLAGS = -DFULL_3c
else ifeq ($(FLASHTYPE), SPARSE)
FLASHFLAGS = -DSPARSE
else ifeq ($(FLASHTYPE), NN)
FLASHFLAGS = -DNN
else
FLASHFLAGS = 
endif

VPATH = ../../FlashCalculation:../../SparseGrids/sg:../../DeepLearning

ifeq ($(FLASHTYPE), FULL_2c)
FLASHSRCS = RST_globalFlashData.F90 
else ifeq ($(FLASHTYPE), FULL_3c)
FLASHSRCS = RST_globalFlashData.F90 
else ifeq ($(FLASHTYPE), SPARSE)
FLASHSRCS = RST_globalFlashData.F90 RST_mathlib.F90 RST_PREOS.F90 RST_comprefac.F90 \
RST_fugacitycoef.F90 RST_pmv.F90 RST_globalSgData.F90 RST_hashTable.F90
else ifeq ($(FLASHTYPE), NN)
FLASHSRCS = RST_globalFlashData.F90 RST_mathlib.F90 RST_PREOS.F90 RST_comprefac.F90 \
RST_fugacitycoef.F90 RST_pmv.F90 RST_getEquResidual.F90
else
FLASHSRCS = RST_globalFlashData.F90 RST_mathlib.F90 RST_PREOS.F90 RST_comprefac.F90 \
RST_fugacitycoef.F90 RST_stability.F90 RST_pmv.F90 RST_flashcalculation.F90 
endif

ifeq ($(FLASHTYPE), FULL_2c)
FSRCS = RST_model.F90 RST_globalData.F90 RST_writeFile.F90 RST_genPlotfileMatlab.F90 \
RST_viscosity.F90 RST_flashcalculation_fullgrid_2c.F90 RST_compositionalTwoPhaseFlow.F90 \
RST_compositionalTwoPhaseFlowDriver.F90 infile_RSTi_compositionalTwoPhaseFlow.F90
else ifeq ($(FLASHTYPE), FULL_3c)
FSRCS = RST_model.F90 RST_globalData.F90 RST_writeFile.F90 RST_genPlotfileMatlab.F90 \
RST_viscosity.F90 RST_flashcalculation_fullgrid_3c.F90 RST_compositionalTwoPhaseFlow.F90 \
RST_compositionalTwoPhaseFlowDriver.F90 infile_RSTi_compositionalTwoPhaseFlow.F90 
else ifeq ($(FLASHTYPE), SPARSE)
FSRCS = RST_model.F90 RST_globalData.F90 RST_writeFile.F90 RST_genPlotfileMatlab.F90 \
RST_viscosity.F90 RST_flashcalculation_sparsegrid.F90 RST_compositionalTwoPhaseFlow.F90 \
RST_compositionalTwoPhaseFlowDriver.F90 infile_RSTi_compositionalTwoPhaseFlow.F90 
else ifeq ($(FLASHTYPE), NN)
FSRCS = RST_model.F90 RST_globalData.F90 RST_writeFile.F90 RST_genPlotfileMatlab.F90 \
RST_viscosity.F90 RST_flashcalculation_nn.F90 RST_compositionalTwoPhaseFlow.F90 \
RST_compositionalTwoPhaseFlowDriver.F90 infile_RSTi_compositionalTwoPhaseFlow.F90 
else
FSRCS = RST_model.F90 RST_globalData.F90 RST_writeFile.F90 RST_genPlotfileMatlab.F90 \
RST_viscosity.F90 RST_compositionalTwoPhaseFlow.F90 RST_compositionalTwoPhaseFlowDriver.F90 \
infile_RSTi_compositionalTwoPhaseFlow.F90
endif
      
%.o: %.F90
	$(FCC) -c $(FCFLAGS) $(FLASHFLAGS) $^ -o $@ 

FLASHOBJS = $(FLASHSRCS:.F90=.o)
FOBJS = $(FSRCS:.F90=.o)

$(PROG): $(FLASHOBJS) $(FOBJS) 
	$(FCC) $(FCFLAGS) $(FLASHFLAGS) $(FLASHOBJS) $(FOBJS) $(LIBS) -o $(PROG)

LIBS = -llapack

all: $(PROG)

.PHONY: clean

clean:
	rm -rf *.mod *.o $(PROG)

run: ${PROG} 
	./${PROG}
