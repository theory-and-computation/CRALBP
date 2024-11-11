#   Usage: vmd -dispdev text -e getPBCSize.tcl -args
#                     [Files]
#

package require pbctools

set INPSF [lindex $argv 0]
set INDCD [lindex $argv 1]

mol new ${INPSF}
mol addfile ${INDCD} waitfor all


set n [molinfo top get numframes]
set output [open "RSMC_PBC_SIZE.txt" a]

for {set i 0} {$i < $n} {incr i} {

    molinfo top set frame $i
    set cell [pbc get -now]
    puts "\t \t progress: $i/$n"

    puts $output $cell
}

puts "\t \t progress: $n/$n"
puts "Done."
puts "output file: RSMC_PBC_SIZE.txt"
close $output

exit

