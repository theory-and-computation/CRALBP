#   Usage: vmd -dispdev text -e rotateProteinToPreferred.tcl -args
#                     [Files]
#

set INPDB [lindex $argv 0]
set YAW [lindex $argv 1]
set PITCH [lindex $argv 2]
set ROLL [lindex $argv 3]
set OUTPDB [lindex $argv 4]

mol new ${INPDB}

#roll - x, pitch - y, yaw - z

#YAW:
set sel [atomselect top "protein"]
set com [measure center $sel weight mass]

set matrix3 [transaxis z ${YAW}]
$sel moveby [vecscale -1.0 $com]
$sel move $matrix3
$sel moveby $com

#PITCH:
set sel [atomselect top "protein"]
set com [measure center $sel weight mass]

set matrix2 [transaxis y ${PITCH}]
$sel moveby [vecscale -1.0 $com]
$sel move $matrix2
$sel moveby $com

#ROLL:
set sel [atomselect top "protein"]
set com [measure center $sel weight mass]

set matrix1 [transaxis x ${ROLL}]
$sel moveby [vecscale -1.0 $com]
$sel move $matrix1
$sel moveby $com


set out [atomselect top "all"]
$out writepdb $OUTPDB

puts "Done."

exit

