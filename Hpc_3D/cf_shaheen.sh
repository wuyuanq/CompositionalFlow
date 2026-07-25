
# @ shell=/bin/csh
# @ job_name = CompositionalFlow_fortran_hpc_3D
# @ account_no = k234
# @ error = error.$(jobid)
# @ output = output.$(jobid)
# @ environment = COPY_ALL;
# @ wall_clock_limit = 2:00:00
# @ notification = always
# @ job_type = bluegene
# @ bg_size = 64
# @ queue

/bgsys/drivers/ppcfloor/bin/mpirun -mode VN -np 80 ./CompositionalFlow_fortran_hpc_3D