package require volutil

set inPSF [lindex $argv 0]
set inDCD [lindex $argv 1]
set tmpDX [lindex $argv 2]
set outDX [lindex $argv 3]

mol new $inPSF
mol addfile $inDCD waitfor all 0

set nframes [molinfo top get numframes]
volmap density [atomselect top "same residue as (segname MEMB and resname POPC and within 3 of protein) and not ions"] -minmax {{-110 -110 -110} {110 110 110}} -res 0.5 -weight mass -allframes -combine avg -mol top -o ${tmpDX}
volutil -smult ${nframes} -o ${outDX} ${tmpDX}

rm ${tmpDX}

exit
