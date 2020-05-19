#! /bin/csh -f
# Thomas Reerink

if($#argv == 0 || $#argv == 1 || $#argv == 3) then

 set number_of_processes = 2
 if($#argv == 1) then
  set number_of_processes = $1
 endif

 set config_file_for_gcm_to_im = 'config-files/oblimap/racmo2.3-to-im-greenland/config-oblimap-racmo2.3-to-im-greenland-20x20km'
 set config_file_for_im_to_gcm = 'config-files/oblimap/im-to-racmo2.3-greenland/config-oblimap-im-to-racmo2.3-greenland-20x20km'
 if($#argv == 3) then
  set config_file_for_gcm_to_im = $2
  set config_file_for_im_to_gcm = $3
 endif


 mpirun -np ${number_of_processes} ./src/oblimap_par_gcm_to_im_program ${config_file_for_gcm_to_im}
 
 # Note that your number of cores is limited (for instance on a laptop) that using both programs at the same time with
 # more than 3 cores each causes instability. In that case better deactivate the second oblimap_par_im_to_gcm_program
 # call below:
 mpirun -np ${number_of_processes} ./src/oblimap_par_im_to_gcm_program ${config_file_for_im_to_gcm}
 

else
 echo ' This script runs without an argument, or with one or with three OPTIONAL arguments, e.g.:'
 echo ' ' $0 '2'
 echo ' Or:'
 echo ' ' $0 '2 config-files/oblimap/racmo2.3-to-im-greenland/config-oblimap-racmo2.3-to-im-greenland-20x20km config-files/oblimap/im-to-racmo2.3-greenland/config-oblimap-im-to-racmo2.3-greenland-20x20km'
endif
