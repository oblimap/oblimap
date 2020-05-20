#! /bin/csh -f
# Thomas Reerink

# Run examples:
# ./oblimap-to-and-fro-mapping.csh config-files/oblimap/racmo2.3-to-im-greenland/config-oblimap-racmo2.3-to-im-greenland-20x20km config-files/oblimap/im-to-racmo2.3-greenland/config-oblimap-im-to-racmo2.3-greenland-20x20km
# ./oblimap-to-and-fro-mapping.csh config-files/oblimap/ccsm-to-im/config-oblimap-ccsm-to-im-greenland-20x20km                   config-files/oblimap/im-to-ccsm/config-oblimap-im-to-ccsm-greenland-20x20km

if($#argv == 0 || $#argv == 2) then

 set config_file_for_gcm_to_im = 'config-files/oblimap/racmo2.3-to-im-greenland/config-oblimap-racmo2.3-to-im-greenland-20x20km'
 set config_file_for_im_to_gcm = 'config-files/oblimap/im-to-racmo2.3-greenland/config-oblimap-im-to-racmo2.3-greenland-20x20km'
 if($#argv == 2) then
  set config_file_for_gcm_to_im = $1
  set config_file_for_im_to_gcm = $2
 endif


 ./src/oblimap_gcm_to_im_program ${config_file_for_gcm_to_im}
 
 # The post profiling:
 if(-e gmon.out) then
  gprof ./src/oblimap_gcm_to_im_program gmon.out > profiling-oblimap_gcm_to_im_program.txt
  rm -f gmon.out
 endif

 ./src/oblimap_im_to_gcm_program ${config_file_for_im_to_gcm}
 
 # The post profiling:
 if(-e gmon.out) then
  gprof ./src/oblimap_im_to_gcm_program gmon.out > profiling-oblimap_im_to_gcm_program.txt
  rm -f gmon.out
 endif


else
 echo ' This script runs without a argument, or requires two OPTIONAL arguments, e.g.:'
 echo ' ' $0
 echo ' Or:'
 echo ' ' $0 'config-files/oblimap/racmo2.3-to-im-greenland/config-oblimap-racmo2.3-to-im-greenland-20x20km config-files/oblimap/im-to-racmo2.3-greenland/config-oblimap-im-to-racmo2.3-greenland-20x20km'
endif
