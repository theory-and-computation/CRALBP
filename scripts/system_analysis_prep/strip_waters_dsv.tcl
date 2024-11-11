set psfFile [lindex $argv 0]
set pdbFile [lindex $argv 1]
set stkFile [lindex $argv 2]
set outPSF [lindex $argv 3]
set outPDB [lindex $argv 4]
set systemName [lindex $argv 5]
set outdir [lindex $argv 6]


############################ NOT COMPLETE


set instk [open $stkFile r]

fconfigure $instk -buffering line
gets $instk dcd


while { $dcd != "" } {

  mol load psf $psfFile pdb $pdbFile
  mol addfile $dcd waitfor all

  set start_frame 0
  set total_frame_num [molinfo top get numframes]
  set end_frame [expr $total_frame_num - 2]


  set line_len [string length $dcd]
  set start_dcd_string_step [expr $line_len - 9]
  set end_dcd_string_step [expr $line_len - 5]
  set dcd_step [string range $dcd $start_dcd_string_step $end_dcd_string_step]

  puts "Trajectory is: $dcd"
  puts "Start frame is: $start_frame. End frame is: $end_frame"

  set outDCD "${systemName}-${dcd_step}-stripped.dcd"

  set sel [atomselect top "not waters"]

  animate read dcd $dcd beg 0 end $end_frame waitfor all top
  animate write dcd $outdir/$outDCD beg 0 end $end_frame waitfor all sel $sel top

  $sel delete
  animate delete all

  gets $instk dcd

}


mol load psf $psfFile pdb $pdbFile

set sel [atomselect top "not waters"]

$sel writepsf $outPSF
$sel writepdb $outPDB


exit