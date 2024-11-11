#   Usage: vmd -dispdev text -e getLipidProtInter-txt.tcl -args
#                     [Files]
#

set INPSF [lindex $argv 0]
set INDCD [lindex $argv 1]
set OUTTXT [lindex $argv 2]

mol new ${INPSF}
mol addfile ${INDCD} waitfor all

set n [molinfo top get numframes]
set output [open ${OUTTXT} a]

for {set i 0} {$i < $n} {incr i} {
    molinfo top set frame $i

    set sel [atomselect top "protein and same residue as (within 3 of segname MEMB)"]
    #set sel [atomselect top "protein"]
    set amino_resids [lsort -unique [$sel get resid]]

    foreach amino_resid $amino_resids {
        set out_line ""
        set amino_select [atomselect top "protein and resid $amino_resid"]

        set amino_resname [lsort -unique [$amino_select get resname]]
        set COORDS [measure center $amino_select weight mass]

        set bound_lipids [atomselect top "same residue as (within 3 of (resid $amino_resid and protein)) and segname MEMB"]
        set lipid_resnames [lsort -unique [$bound_lipids get resname]]

        set new_i [expr $i + 66525]
        #set new_i $i

        append out_line "$amino_resid,$amino_resname,$lipid_resnames,$new_i,$COORDS"
        puts $output $out_line
    }
    puts "\t \t progress: $i/$n"
}

puts "\t \t progress: $n/$n"
puts "Done."
puts "output file: ${OUTTXT}"
close $output
exit
