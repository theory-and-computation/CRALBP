#!/usr/bin/python
import os

def genConf(conf_base_name, angstrom_window) -> None:
 """ Generates an LSF submission script and returns the filename """

 script_name = f"{conf_base_name}-{angstrom_window}.conf"
 output = open(script_name, 'w')
 output.write(r"""#############################################################
## JOB DESCRIPTION                                         ##
#############################################################

# CRALBP Steered Dynamics - lower pore

#############################################################
## ADJUSTABLE PARAMETERS                                   ##
#############################################################

set basename       $env(basename)
set nsteps         $env(nsteps)
set outputname	   $env(output)
set firstjob       $env(firstjob)
set temperature    $env(temp)

set lPP            $env(langPistPeriod)
set lPD            $env(langPistDecay)
set rstFreq        $env(rstFreq)

structure          ${basename}.psf
coordinates        ${basename}.pdb

if { $firstjob == "yes" } {
  temperature        ${temperature}
  bincoordinates     ${basename}.coor
}
if { $firstjob == "no" } {
  set rstname	       $env(rstname)
  bincoordinates	   ${rstname}.coor
  binvelocities	     ${rstname}.vel
  extendedSystem	   ${rstname}.xsc
}

firsttimestep      0

#############################################################
## SIMULATION PARAMETERS                                   ##
#############################################################

# Input
paraTypeCharmm	    on
parameters          ff/par_all36m_prot.prm
parameters          ff/par_all36_na.prm
parameters          ff/par_all36_carb.prm
parameters          ff/par_all36_lipid.prm
parameters          ff/par_all36_cgenff.prm
parameters          ff/rta.prm
parameters          ff/toppar_water_ions_namd.str

# Force-Field Parameters
exclude             scaled1-4
1-4scaling          1.0
cutoff              12.0    ;# ctofnb
switching           on
switchdist          10.0    ;# ctonnb
pairlistdist        15.0    ;# cutnb


# Integrator Parameters
timestep            2.0  ;# 2fs/step
rigidBonds          all  ;# needed for 2fs steps
nonbondedFreq       1
fullElectFrequency  2  
stepspercycle       10


# Constant Temperature Control
langevin            on    ;# do langevin dynamics
langevinDamping     1     ;# damping coefficient (gamma) of 1/ps
langevinTemp        $temperature
langevinHydrogen    off    ;# don't couple langevin bath to hydrogens


# Periodic Boundary Conditions
if { $firstjob == "yes"} {
  cellBasisVector1    73.0    0.0    0.0
  cellBasisVector2     0.0   73.0    0.0
  cellBasisVector3     0.0    0.0   73.0
  cellOrigin           0.0    0.0    0.0
}

wrapAll             on


# PME (for full-system periodic electrostatics)
PME                 yes
PMEGridSpacing      1.0
PMEInterpOrder      8


# Constant Pressure Control (variable volume)
useGroupPressure      yes ;# needed for rigidBonds
useFlexibleCell       no
useConstantRatio      no

langevinPiston        on
langevinPistonTarget  1.01325 ;#  in bar -> 1 atm
langevinPistonPeriod  $lPP   # 100.0
langevinPistonDecay   $lPD   # 50.0
langevinPistonTemp    $temperature
langevinPistonBarrier off

# Output
ldBalancer          hybrid
ldbPeriod           20000
outputName          $outputname
outputTiming        10000

restartfreq         $rstFreq  ;# 500steps = every 1ps
dcdfreq             1000     ;# 250steps
#forceDCDfreq        10000
xstFreq             10000     ;# 250 steps
outputEnergies      10000     ;# 100 steps
outputPressure      10000

twoAwayX yes 
twoAwayY yes

#############################################################
## EXTRA PARAMETERS                                        ##
#############################################################

colvars		on
colvarsConfig	3hy5_lowerpore-%d.cvc	


#############################################################
## EXECUTION SCRIPT                                        ##
#############################################################

run ${nsteps}
 """% angstrom_window)
 output.close()
 os.system('chmod +x %s' % script_name)

def genCVC(cvc_base_name, angstrom_window) -> None:
 """ Generates an LSF submission script and returns the filename """

 script_name = f"{cvc_base_name}-{angstrom_window}.cvc"
 output = open(script_name, 'w')
 output.write(r"""colvarsTrajFrequency    1

#Steered MD CRALBP ; Group 1 (protein residues) | Group 2 (RTA ligand)

colvar {
  name steeringDyn
  distanceZ {
    ref {
      atomsFile       restraint_atom_select/3hy5_RTA_NEW_LP_RESTRAINT.pdb
      atomsCol        B
      atomsColValue   0.10
    }
    main {
      atomsFile       restraint_atom_select/3hy5_RTA_NEW_LP_RESTRAINT.pdb

      atomsCol        B
      atomsColValue   0.20
    }
    axis {
      (1.0, 0.0, 0.0)
    }
  }
}

harmonic {
  name center_of_mass1
  colvars steeringDyn
  centers %d # initial distance
  forceConstant 4
}
""" % angstrom_window)
 output.close()
 os.system('chmod +x %s' % script_name)

#--------------#
# Main Routine #
#--------------#

if __name__ == '__main__':
 start = 1
 end  = 30
 basecvc = '3hy5_lowerpore'
 baseconf = 'umb_lp'


 for window in range(start, end + 1):
  genCVC(window, basecvc)
  genConf(window, baseconf)

