
#!/bin/sh
#@ job_name         = CompositionalFlow_fortran_hpc
#@ output           = $(job_name).$(jobid).out
#@ error            = $(job_name).$(jobid).err
#@ job_type         = parallel
#@ environment      = COPY_ALL
#@ wall_clock_limit = 01:00:00
#@ node             = 1
#@ tasks_per_node   = 4
#@ queue

$MPIEXEC -x LD_LIBRARY_PATH -np $LOADL_TOTAL_TASKS ./CompositionalFlow_fortran_hpc