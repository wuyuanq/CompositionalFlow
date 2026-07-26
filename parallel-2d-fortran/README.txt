
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

----------------------------------------------------------------------------------------

This is the parallel program to simulate the 2D two-phase compositional flow. Any number of components is acceptable and what you have to pay attention on is to set your input file correctly. The number of processors can be set as 1 if you don't want to use the parallel advantage. All the files are self-explained. Of course, you can also use the program to simulate the single-phase flow.

----------------------------------------------------------------------------------------

1. When you set the parameters of your case, you have to keep the number of the processors appearing in the following three places the same: 
In the input file:
	call proceAlloc(4, model1%nx, model1%ny, model1%pncols, model1%pnrows)
In the makefile:
	run: ${PROG} 
		mpirun -np 4 ./${PROG}
In the Shaheen and Neser shell files:
	/bgsys/drivers/ppcfloor/bin/mpirun -mode VN -np 4 ./CompositionalFlow_fortran_hpc

----------------------------------------------------------------------------------------

Before you run the program on Shaheen, there are some notices you have to pay attention to:

1. You have to set the module environment on Shaheen correctly. You can use the following commands to download proper modules:

module load bluegene

2. You have to create the directory to store the results of your program at first, since there is no system call from FORTRAN on Shaheen. For example, before you run infile1, you can use "mkdir case1" in the terminal at first.

-----------------------------------------------------------------------------------------

Before you run the program on Neser, there are some notices you have to pay attention to:

1. You have to set the module environment on Neser correctly. You can use the following commands to download proper modules:

module unload openmpi/1.5.4/gcc
module load intel-compilers/11.1
module load openmpi/1.6.4/intel

2. You have to create the directory to store the results of your program at first, since there is no system call from FORTRAN on Shaheen. For example, before you run infile1, you can use "mkdir case1" in the terminal at first.

3. When your program size is large, for example, 160*160 cells, you have to consider the issue of "insufficient virtual memory" on each Neser node. In such condition, you can allocate more nodes to the program and reduce the number of processes on each node.

-----------------------------------------------------------------------------------------

You can choose the option "FLASHTYPE" in the Makefile to decide whether to use Sparse Grids to do the flash calculations. If "FLASHTYPE = SPARSE", it means that the program will use Sparse Grids, otherwise the program will use true flash calculations. 

If you decide to use Sparse Grids, you have to generate the Sparse Grids surrogate at first and then put the surrogate files into "CompositionalFlow_fortran/Fullgrid". The Sparse Grids only supports at most 3 components now. The sample 2-component surrogates are  put into "CompositionalFlow_fortran/Fullgrid_2c" and the sample 3-component surrogates are  put into "CompositionalFlow_fortran/Fullgrid_3c". Before you run the 2-component case, you have to copy the files in "CompositionalFlow_fortran/Fullgrid_2c" to "CompositionalFlow_fortran/Fullgrid"; before you run the 3-component case, you have to copy the files in "CompositionalFlow_fortran/Fullgrid_3c" to "CompositionalFlow_fortran/Fullgrid".

It is also your responsibility to set the variable "TABLESIZE" in the file "RST_compositionalTwoPhaseFlow.F90" correctly. "TABLESIZE" should be equal with the size of the Sparse Grids surrogate. For example, if there are only 2 components and there are 129 points in p direction and 129 points in z1 direction, then TABLESIZE = 129*129. If there are 3 components and there are 17 points in p direction and 257 points in z1 direction and 257 points in z2 direction, then TABLESIZE = 17*(1+257)/2*257.  

If you have 2 components in the input file, you have to rename the file "RST_dummyFlashcalculation-2c.F90" in "CompositionalFlow_fortran" as "RST_dummyFlashcalculation.F90"; If you have 3 components in the input file, you have to rename the file "RST_dummyFlashcalculation-3c.F90" in "CompositionalFlow_fortran" as "RST_dummyFlashcalculation.F90". 

-----------------------------------------------------------------------------------------

Before you run the cases, you have to remove the numbers at the tail of the name of the input files at first. 

-----------------------------------------------------------------------------------------

You can see the results in the document "case*", and "*" is a number which depends on your input file. In "case*", there is a file called "matlabplot.m". You can run the file to generate the matlab figures in the document "matlabplots".

