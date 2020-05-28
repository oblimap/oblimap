#! /bin/csh -f
# Thomas Reerink
#
# ${1} the first  argument is the program name, e.g.: fortran_example_model_using_namelist_program
# ${2} the second argument is the config file
# ${3} the third  optional argument is some extra label
# If a last optional argument equals 'bg' the job will run in the background
#
# Run example:
#  ./long-local-run.csh fortran_mpi_win_allocate_shared_example_program config-files/config-mpi-example
#  ./long-local-run.csh fortran_mpi_win_allocate_shared_example_program config-files/config-mpi-example test bg
#

set processes_per_node = 2


if(${processes_per_node} > 1) then
 set run_parallel          = 'mpiexec -n '${processes_per_node}
else
 set run_parallel          = ''
endif


if($#argv > 1) then
 set svn_revision = `svn info | grep Revision | sed -e 's/Revision: //'`

 set background_job = 'no'
 # Check whether the last command line argument equals 'bg' ($#argv is the number of command line arguments)
 if($argv[$#argv] == 'bg') then
  set background_job = 'yes'
  echo $#argv $argv[$#argv]
 endif

 set program_name            = ${1}
 set source_path_config_file = ${2}
 set label                   = ${svn_revision}
 if($#argv > 2)   set label  = ${svn_revision}'-'${3}
 if(${3} == 'bg') set label  = ${svn_revision}

 # Check if program is available:
 if(-e ./src/${program_name}) then

  set running_program_name = ${program_name}'_r'${label}

  # Extracting the relevant label info of the input config file name:
  set label_1              = `echo ${source_path_config_file} | sed -e 's/^.*\///' -e 's/config_//' -e 's/config-//'`
  set config_name          = 'config_'${label_1}'_r'${label}
  set running_dir          = ${label_1}'-r'${label}

  mkdir -p ${running_dir}
  cp -f ./src/${program_name}      ${running_dir}/${running_program_name}
  cp -f ${source_path_config_file} ${running_dir}/${config_name}

  # COPY THE SVN DIFFERENCES RELATIVE TO THE CURRENT REVISION:
  svn diff ./src/*.f90 > ${running_dir}/svn-diff-of-source-files.txt
 #mkdir -p ${running_dir}/src-copy/
 #rsync -a ./src/*.f90 ${running_dir}/src-copy/

  cd ${running_dir}
  if(${background_job} == 'yes') then
   # Running a job in the background:
   ${run_parallel} ./${running_program_name} ${config_name} > output.txt ; if(-e gmon.out) gprof ${running_program_name} gmon.out > time_measurements.txt ; if(-e gmon.out) rm -f gmon.out ; if(-e fort.2) rm -f fort.2 &

   echo ' The results can be found in the directory  ' ${running_dir}
  else
   # Running a job in the foreground:
   ${run_parallel} ./${running_program_name} ${config_name}

   if(-e gmon.out) then
    gprof ${running_program_name} gmon.out > time_measurements.txt
    rm -f gmon.out
   endif
   if(-e fort.2) rm -f fort.2

   echo  ''
   echo ' The results can be found in the directory  ' ${running_dir}
  #echo ' The results can be found in the directory  ' ${running_dir} '  and can be viewed with:'
  #ls -1 *.nc |sed -e 's/^/  ncview '${running_dir}'\//'
  endif
  echo  ''

  echo ' This run has been executed by running the script '$0' script like:' >  log-executed-run-command.txt
  echo '  ' $0 $*                                                            >> log-executed-run-command.txt
  echo '  '                                                                  >> log-executed-run-command.txt
  echo ' And the fortran source has been runned in this directory by:'       >> log-executed-run-command.txt
  echo '  ./'${running_program_name} ${config_name}                          >> log-executed-run-command.txt

 else
  echo ' The program ' ${program_name} ' is not found. Compile it with:  make' ${program_name}
 endif

else
 echo ' Needs two arguments (an extra label can be given with the optional third argument). If the last argument equals "bg" the job will run in the background, e.g.:'
 echo
 echo ' Examples running the job in the background:'
 echo '  '$0' fortran_mpi_win_allocate_shared_example_program config-files/config-mpi-example test bg'
 echo '  '$0' fortran_mpi_win_allocate_shared_example_program config-files/config-mpi-example bg'
 echo
 echo ' Examples running the job in the foreground:'
 echo '  '$0' fortran_mpi_win_allocate_shared_example_program config-files/config-mpi-example test'
 echo '  '$0' fortran_mpi_win_allocate_shared_example_program config-files/config-mpi-example'
endif
