# Orients the protein principal axes to the (x, y, z) axes
#   and recenters the center to the origin
#
# Usage: vmd -dispdev text -e orient.tcl -args
#          [inputPDB] [outputPDB]

set inputPDB  [lindex $argv 0]
set outputPDB [lindex $argv 1]

mol new ${inputPDB}

#
# Orient the principal axes of protein to (x, y, z)
#

lappend auto_path /Users/danielsantos/Downloads/la1.0
lappend auto_path /Users/danielsantos/Downloads/orient

package require Orient
namespace import Orient::orient

set sel [atomselect top "all"]
set I [draw principalaxes $sel]
set A [orient $sel [lindex $I 2] {0 0 1}]
$sel move $A
set I [draw principalaxes $sel]
set A [orient $sel [lindex $I 1] {0 1 0}]
$sel move $A
set I [draw principalaxes $sel]

#
# Recenter the protein to origin
#
set minmax [measure minmax $sel]

set min [lindex $minmax 0]
set max [lindex $minmax 1]

set center [vecscale 0.5 [vecadd $min $max]]
$sel moveby [vecinvert $center]

# Write out the molecule
$sel writepdb ${outputPDB}

exit
