#   Usage: vmd -dispdev text -e getRotationalMatrices.tcl -args
#                     [Files]
#

set FILE [lindex $argv 0]
set INPSF [lindex $argv 1]
set INDCD [lindex $argv 2]
set REFPDB [lindex $argv 3]

mol new ${REFPDB}
mol new ${INPSF}
mol addfile ${INDCD} waitfor all

set n [molinfo top get numframes]
set output [open ${FILE} a+]

for {set i 0} {$i < $n} {incr i} {

    molinfo top set frame $i

    #Aligning membrane:
    set sel1 [atomselect 0 "segname MEMB"]
    set sel2 [atomselect 1 "segname MEMB"]
    set transformation_matrix [measure fit $sel2 $sel1]
    set move_sel [atomselect 1 "all"]
    $move_sel move $transformation_matrix

    #Obtaining transformation matrix required to align the proteins
    set sel1 [atomselect 0 "protein and backbone"]
    set sel2 [atomselect 1 "protein and backbone"]
    set transformation_matrix [measure fit $sel2 $sel1]
    puts $output $transformation_matrix
    puts "\t \t progress: $i/$n"
}

puts "\t \t progress: $n/$n"
puts "Done."
puts "output file: ${FILE}"
close $output

exit

