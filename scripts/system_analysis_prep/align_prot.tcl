package require pbctools

set psfFile [lindex $argv 0]
set dcdFile [lindex $argv 1]
set refPDB [lindex $argv 2]
set outDCD [lindex $argv 3]

mol new $psfFile
mol addfile $dcdFile waitfor all 


pbc set {281 281 281} -all

#pbc wrap -all -center com -centersel "protein" -compound fragment
#pbc wrap -all -center com -centersel "resname RTA" -compound fragment

pbc wrap -all -center com -centersel "segname MEMB" - compound fragment
#pbc wrap -all center com -centersel "segname MEMB" - compound fragment
#pbc wrap -all -center com -centersel "resid 18 and fragment 43 to 48" -compound fragment
#pbc wrap -all -center com -centersel "resid 18 and fragment 50 to 55" -compound fragment
#pbc wrap -all -center com -centersel "resid 18 and fragment 57 to 62" -compound fragment
#pbc wrap -all -center com -centersel "resid 18 and fragment 64 to 69" -compound fragment
#pbc wrap -all -center com -centersel "resid 18 and fragment 71 to 76" -compound fragment
#pbc wrap -all -center com -centersel "resid 18 and fragment 42 49 56 63 70" -compound fragment

#pbc wrap -all -center com -centersel "resid 18 and protein" -compound fragment

####TESTING FOR MEMB####

mol new $refPDB
set nframes [molinfo 0 get numframes]
#set ref_sel [atomselect 1 "protein"]
set ref_sel [atomselect 1 "segname MEMB"]

for {set i 0} {$i < $nframes} {incr i} {
  #set align_sel [atomselect 0 "protein" frame $i]
  set align_sel [atomselect 0 "segname MEMB"]
  set all0 [atomselect 0 all frame $i]

  set M [measure fit $align_sel $ref_sel]
  $all0 move $M

  $align_sel delete
  $all0 delete
}

animate write dcd $outDCD waitfor all 0

exit
