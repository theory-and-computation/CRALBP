package require volutil

set inDX1 [lindex $argv 0]
set inDX2 [lindex $argv 1]
set outDX [lindex $argv 2]

volutil -add -o $outDX $inDX1 $inDX2

exit
