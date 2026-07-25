
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

-----------------------------------------------------------------------------------------

This is the serial program to simulate the 2D two-phase compositional flow. Any number of components is acceptable and what you have to pay attention on is to set your input file correctly. All the files are self-explained. Of course, you can also use the program to simulate the single-phase flow.

-----------------------------------------------------------------------------------------

You can choose the option "FLASHTYPE" in the Makefile to decide whether to use Sparse Grids to do the flash calculations. If "FLASHTYPE = SPARSE", it means that the program will use Sparse Grids, otherwise the program will use true flash calculations. 

If you decide to use Sparse Grids, you have to generate the Sparse Grids surrogate at first and then put the surrogate files into "Fullgrid". The Sparse Grids only supports at most 3 components now. The sample 2-component surrogates are  put into "Fullgrid_2c" and the sample 3-component surrogates are  put into "Fullgrid_3c". Before you run the 2-component case, you have to copy the files in "Fullgrid_2c" to "Fullgrid"; before you run the 3-component case, you have to copy the files in "Fullgrid_3c" to "Fullgrid".

It is also your responsibility to set the variable "TABLESIZE" in the file "RST_compositionalTwoPhaseFlow.F90" correctly. "TABLESIZE" should be equal with the size of the Sparse Grids surrogate. For example, if there are only 2 components and there are 129 points in p direction and 129 points in z1 direction, then TABLESIZE = 129*129. If there are 3 components and there are 17 points in p direction and 257 points in z1 direction and 257 points in z2 direction, then TABLESIZE = 17*(1+257)/2*257.  

If you have 2 components in the input file, you have to rename the file "RST_dummyFlashcalculation-2c.F90" as "RST_dummyFlashcalculation.F90"; If you have 3 components in the input file, you have to rename the file "RST_dummyFlashcalculation-3c.F90" as "RST_dummyFlashcalculation.F90". 

-----------------------------------------------------------------------------------------

Before you run the cases, you have to remove the numbers at the tail of the name of the input files at first. 

-----------------------------------------------------------------------------------------

You can see the results in the document "case*", and "*" is a number which depends on your input file. In "case*", there is a file called "matlabplot.m". You can run the file to generate the matlab figures in the document "matlabplots".
