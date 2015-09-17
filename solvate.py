#!/usr/local/bin/python2.7
# encoding: utf-8
'''
solvate -- run the bevanlab Justin Lemkul's solvation free energy tutorial, and
more!

@author:     jtp4kc, Tyler P

@copyright:  2015 Shirts Group, University of Virginia. All rights reserved.

@contact:    jtp4kc@virginia.edu
@deffield    updated: Updated
'''

import sys
import os
import math

from argparse import ArgumentParser
# from argparse import RawDescriptionHelpFormatters

__all__ = []
__version__ = 0.1
__date__ = '2015-09-04'
__updated__ = '2015-09-04'

DEBUG = 0
TESTRUN = 0
PROFILE = 0

class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg
    def __str__(self):
        return self.msg
    def __unicode__(self):
        return self.msg

def job_sh():
    return """#!/bin/bash

# Set some environment variables 
FREE_ENERGY=/home/justin/Free_Energy
echo "Free energy home directory set to $FREE_ENERGY"

MDP=$FREE_ENERGY/MDP
echo ".mdp files are stored in $MDP"

LAMBDA=0

# A new directory will be created for each value of lambda and
# at each step in the workflow for maximum organization.

mkdir Lambda_$LAMBDA
cd Lambda_$LAMBDA

#################################
# ENERGY MINIMIZATION 1: STEEP  #
#################################
echo "Starting minimization for lambda = $LAMBDA..." 

mkdir EM_1 
cd EM_1

# Iterative calls to grompp and mdrun to run the simulations

grompp -f $MDP/EM/em_steep_$LAMBDA.mdp -c $FREE_ENERGY/Methane/methane_water.gro -p $FREE_ENERGY/Methane/topol.top -o min$LAMBDA.tpr

mdrun -nt 2 -deffnm min$LAMBDA


#################################
# ENERGY MINIMIZATION 2: L-BFGS #
#################################

cd ../
mkdir EM_2
cd EM_2

grompp -f $MDP/EM/em_l-bfgs_$LAMBDA.mdp -c ../EM_1/min$LAMBDA.gro -p $FREE_ENERGY/Methane/topol.top -o min$LAMBDA.tpr

# Run L-BFGS in serial (cannot be run in parallel)

mdrun -nt 1 -deffnm min$LAMBDA

echo "Minimization complete."



#####################
# NVT EQUILIBRATION #
#####################
echo "Starting constant volume equilibration..."

cd ../
mkdir NVT
cd NVT

grompp -f $MDP/NVT/nvt_$LAMBDA.mdp -c ../EM_2/min$LAMBDA.gro -p $FREE_ENERGY/Methane/topol.top -o nvt$LAMBDA.tpr

mdrun -nt 2 -deffnm nvt$LAMBDA

echo "Constant volume equilibration complete."



#####################
# NPT EQUILIBRATION #
#####################
echo "Starting constant pressure equilibration..."

cd ../
mkdir NPT
cd NPT

grompp -f $MDP/NPT/npt_$LAMBDA.mdp -c ../NVT/nvt$LAMBDA.gro -p $FREE_ENERGY/Methane/topol.top -t ../NVT/nvt$LAMBDA.cpt -o npt$LAMBDA.tpr

mdrun -nt 2 -deffnm npt$LAMBDA

echo "Constant pressure equilibration complete."


#################
# PRODUCTION MD #
#################
echo "Starting production MD simulation..."

cd ../
mkdir Production_MD
cd Production_MD

grompp -f $MDP/Production_MD/md_$LAMBDA.mdp -c ../NPT/npt$LAMBDA.gro -p $FREE_ENERGY/Methane/topol.top -t ../NPT/npt$LAMBDA.cpt -o md$LAMBDA.tpr

mdrun -nt 2 -deffnm md$LAMBDA

echo "Production MD complete."

# End
echo "Ending. Job completed for lambda = $LAMBDA"
"""

def em_steep_MDP():
    return """; Run control
integrator               = steep 
nsteps                   = 5000
; EM criteria and other stuff
emtol                    = 100
emstep                   = 0.01
niter                    = 20
nbfgscorr                = 10
; Output control
nstlog                   = 1
nstenergy                = 1
; Neighborsearching and short-range nonbonded interactions
nstlist                  = 1
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.0
; Electrostatics
coulombtype              = PME
rcoulomb                 = 1.0
; van der Waals
vdw-type                 = switch
rvdw-switch              = 0.8
rvdw                     = 0.9
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                  = EnerPres
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; EWALD/PME/PPPM parameters
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
optimize_fft             = no
; Temperature and pressure coupling are off during EM
tcoupl                   = no
pcoupl                   = no
; Free energy control stuff
free_energy              = yes
init_lambda              = 0.0
delta_lambda             = 0
foreign_lambda           = 0.05
sc-alpha                 = 0.5
sc-power                 = 1.0
sc-sigma                 = 0.3 
couple-moltype           = Methane  ; name of moleculetype to decouple
couple-lambda0           = vdw      ; only van der Waals interactions
couple-lambda1           = none     ; turn off everything, in this case only vdW
couple-intramol          = no
nstdhdl                  = 10
; Generate velocities to start
gen_vel                  = no 
; options for bonds
constraints              = h-bonds  ; we only have C-H bonds here
; Type of constraint algorithm
constraint-algorithm     = lincs
; Do not constrain the starting configuration
continuation             = no
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 12
"""

def em_lbfgs_MDP():
    return """; Run control
integrator               = l-bfgs
nsteps                   = 5000
define                   = -DFLEXIBLE
; EM criteria and other stuff
emtol                    = 100
emstep                   = 0.01
niter                    = 20
nbfgscorr                = 10
; Output control
nstlog                   = 1
nstenergy                = 1
; Neighborsearching and short-range nonbonded interactions
nstlist                  = 1
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.0
; Electrostatics
coulombtype              = PME
rcoulomb                 = 1.0
; van der Waals
vdw-type                 = switch
rvdw-switch              = 0.8
rvdw                     = 0.9
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                  = EnerPres
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; EWALD/PME/PPPM parameters
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
optimize_fft             = no
; Temperature and pressure coupling are off during EM
tcoupl                   = no
pcoupl                   = no
; Free energy control stuff
free_energy              = yes
init_lambda              = 0.0
delta_lambda             = 0
foreign_lambda           = 0.05
sc-alpha                 = 0.5
sc-power                 = 1.0
sc-sigma                 = 0.3 
couple-moltype           = Methane  ; name of moleculetype to decouple
couple-lambda0           = vdw      ; only van der Waals interactions
couple-lambda1           = none     ; turn off everything, in this case only vdW
couple-intramol          = no
nstdhdl                  = 10
; Generate velocities to start
gen_vel                  = no 
; options for bonds
constraints              = none     ; L-BFGS doesn't work with constraints 
; Type of constraint algorithm
constraint-algorithm     = lincs
; Do not constrain the starting configuration
continuation             = no
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 12
"""

def nvt_MDP():
    return """; Run control
integrator               = sd       ; Langevin dynamics
tinit                    = 0
dt                       = 0.002
nsteps                   = 50000    ; 100 ps
nstcomm                  = 100
; Output control
nstxout                  = 500
nstvout                  = 500
nstfout                  = 0
nstlog                   = 500
nstenergy                = 500
nstxtcout                = 0
xtc-precision            = 1000
; Neighborsearching and short-range nonbonded interactions
nstlist                  = 10
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.0
; Electrostatics
coulombtype              = PME
rcoulomb                 = 1.0
; van der Waals
vdw-type                 = switch
rvdw-switch              = 0.8
rvdw                     = 0.9
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                  = EnerPres
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; EWALD/PME/PPPM parameters
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
optimize_fft             = no
; Temperature coupling
; tcoupl is implicitly handled by the sd integrator
tc_grps                  = system
tau_t                    = 1.0
ref_t                    = 300
; Pressure coupling is off for NVT
Pcoupl                   = No
tau_p                    = 0.5
compressibility          = 4.5e-05
ref_p                    = 1.0 
; Free energy control stuff
free_energy              = yes
init_lambda              = 0.0
delta_lambda             = 0
foreign_lambda           = 0.05
sc-alpha                 = 0.5
sc-power                 = 1.0
sc-sigma                 = 0.3 
couple-moltype           = Methane  ; name of moleculetype to decouple
couple-lambda0           = vdw      ; only van der Waals interactions
couple-lambda1           = none     ; turn off everything, in this case only vdW
couple-intramol          = no
nstdhdl                  = 10
; Generate velocities to start
gen_vel                  = yes
gen_temp                 = 300
gen_seed                 = -1
; options for bonds
constraints              = h-bonds  ; we only have C-H bonds here
; Type of constraint algorithm
constraint-algorithm     = lincs
; Do not constrain the starting configuration
continuation             = no
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 12
"""

def npt_MDP():
    return """; Run control
integrator               = sd       ; Langevin dynamics
tinit                    = 0
dt                       = 0.002
nsteps                   = 50000    ; 100 ps
nstcomm                  = 100
; Output control
nstxout                  = 500
nstvout                  = 500
nstfout                  = 0
nstlog                   = 500
nstenergy                = 500
nstxtcout                = 0
xtc-precision            = 1000
; Neighborsearching and short-range nonbonded interactions
nstlist                  = 10
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.0
; Electrostatics
coulombtype              = PME
rcoulomb                 = 1.0
; van der Waals
vdw-type                 = switch
rvdw-switch              = 0.8
rvdw                     = 0.9
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                  = EnerPres
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; EWALD/PME/PPPM parameters
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
optimize_fft             = no
; Temperature coupling
; tcoupl is implicitly handled by the sd integrator
tc_grps                  = system
tau_t                    = 1.0
ref_t                    = 300
; Pressure coupling is on for NPT
Pcoupl                   = Parrinello-Rahman 
tau_p                    = 0.5
compressibility          = 4.5e-05
ref_p                    = 1.0 
; Free energy control stuff
free_energy              = yes
init_lambda              = 0.0
delta_lambda             = 0
foreign_lambda           = 0.05
sc-alpha                 = 0.5
sc-power                 = 1.0
sc-sigma                 = 0.3 
couple-moltype           = Methane  ; name of moleculetype to decouple
couple-lambda0           = vdw      ; only van der Waals interactions
couple-lambda1           = none     ; turn off everything, in this case only vdW
couple-intramol          = no
nstdhdl                  = 10
; Do not generate velocities
gen_vel                  = no 
; options for bonds
constraints              = h-bonds  ; we only have C-H bonds here
; Type of constraint algorithm
constraint-algorithm     = lincs
; Constrain the starting configuration
; since we are continuing from NVT
continuation             = yes 
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 12
"""

def md_MDP():
    return """; Run control
integrator               = sd       ; Langevin dynamics
tinit                    = 0
dt                       = 0.002
nsteps                   = 2500000  ; 5 ns
nstcomm                  = 100
; Output control
nstxout                  = 500
nstvout                  = 500
nstfout                  = 0
nstlog                   = 500
nstenergy                = 500
nstxtcout                = 0
xtc-precision            = 1000
; Neighborsearching and short-range nonbonded interactions
nstlist                  = 10
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.0
; Electrostatics
coulombtype              = PME
rcoulomb                 = 1.0
; van der Waals
vdw-type                 = switch
rvdw-switch              = 0.8
rvdw                     = 0.9
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                  = EnerPres
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; EWALD/PME/PPPM parameters
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
optimize_fft             = no
; Temperature coupling
; tcoupl is implicitly handled by the sd integrator
tc_grps                  = system
tau_t                    = 1.0
ref_t                    = 300
; Pressure coupling is on for NPT
Pcoupl                   = Parrinello-Rahman 
tau_p                    = 0.5
compressibility          = 4.5e-05
ref_p                    = 1.0 
; Free energy control stuff
free_energy              = yes
init_lambda              = 0.0
delta_lambda             = 0
foreign_lambda           = 0.05
sc-alpha                 = 0.5
sc-power                 = 1.0
sc-sigma                 = 0.3 
couple-moltype           = Methane  ; name of moleculetype to decouple
couple-lambda0           = vdw      ; only van der Waals interactions
couple-lambda1           = none     ; turn off everything, in this case only vdW
couple-intramol          = no
nstdhdl                  = 10
; Do not generate velocities
gen_vel                  = no 
; options for bonds
constraints              = h-bonds  ; we only have C-H bonds here
; Type of constraint algorithm
constraint-algorithm     = lincs
; Constrain the starting configuration
; since we are continuing from NPT
continuation             = yes 
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 12
"""

def em_steep_mdp(lam="0.0", fol="0.05", mol="Methane", cp1="vdw", cp2="none"):
    return """; Run control
integrator               = steep 
nsteps                   = 5000
; EM criteria and other stuff
emtol                    = 100
emstep                   = 0.01
niter                    = 20
nbfgscorr                = 10
; Output control
nstlog                   = 1
nstenergy                = 1
; Neighborsearching and short-range nonbonded interactions
nstlist                  = 1
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.0
; Electrostatics
coulombtype              = PME
rcoulomb                 = 1.0
; van der Waals
vdw-type                 = switch
rvdw-switch              = 0.8
rvdw                     = 0.9
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                  = EnerPres
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; EWALD/PME/PPPM parameters
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
optimize_fft             = no
; Temperature and pressure coupling are off during EM
tcoupl                   = no
pcoupl                   = no
; Free energy control stuff
free_energy              = yes
init_lambda              = {lam}
delta_lambda             = 0
foreign_lambda           = {fol}
sc-alpha                 = 0.5
sc-power                 = 1.0
sc-sigma                 = 0.3 
couple-moltype           = {mol}  ; name of moleculetype to decouple
couple-lambda0           = {cp1}      ; only van der Waals interactions
couple-lambda1           = {cp2}     ; turn off everything, in this case only vdW
couple-intramol          = no
nstdhdl                  = 10
; Generate velocities to start
gen_vel                  = no 
; options for bonds
constraints              = h-bonds  ; we only have C-H bonds here
; Type of constraint algorithm
constraint-algorithm     = lincs
; Do not constrain the starting configuration
continuation             = no
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 12
""".format(lam=lam, fol=fol, mol=mol, cp1=cp1, cp2=cp2)

def em_lbfgs_mdp(lam="0.0", fol="0.05", mol="Methane", cp1="vdw", cp2="none"):
    return """; Run control
integrator               = l-bfgs
nsteps                   = 5000
define                   = -DFLEXIBLE
; EM criteria and other stuff
emtol                    = 100
emstep                   = 0.01
niter                    = 20
nbfgscorr                = 10
; Output control
nstlog                   = 1
nstenergy                = 1
; Neighborsearching and short-range nonbonded interactions
nstlist                  = 1
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.0
; Electrostatics
coulombtype              = PME
rcoulomb                 = 1.0
; van der Waals
vdw-type                 = switch
rvdw-switch              = 0.8
rvdw                     = 0.9
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                  = EnerPres
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; EWALD/PME/PPPM parameters
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
optimize_fft             = no
; Temperature and pressure coupling are off during EM
tcoupl                   = no
pcoupl                   = no
; Free energy control stuff
free_energy              = yes
init_lambda              = {lam}
delta_lambda             = 0
foreign_lambda           = {fol}
sc-alpha                 = 0.5
sc-power                 = 1.0
sc-sigma                 = 0.3 
couple-moltype           = {mol}  ; name of moleculetype to decouple
couple-lambda0           = {cp1}      ; only van der Waals interactions
couple-lambda1           = {cp2}     ; turn off everything, in this case only vdW
couple-intramol          = no
nstdhdl                  = 10
; Generate velocities to start
gen_vel                  = no 
; options for bonds
constraints              = none     ; L-BFGS doesn't work with constraints 
; Type of constraint algorithm
constraint-algorithm     = lincs
; Do not constrain the starting configuration
continuation             = no
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 12
""".format(lam=lam, fol=fol, mol=mol, cp1=cp1, cp2=cp2)

def nvt_mdp(lam="0.0", fol="0.05", mol="Methane", cp1="vdw", cp2="none"):
    return """; Run control
integrator               = sd       ; Langevin dynamics
tinit                    = 0
dt                       = 0.002
nsteps                   = 50000    ; 100 ps
nstcomm                  = 100
; Output control
nstxout                  = 500
nstvout                  = 500
nstfout                  = 0
nstlog                   = 500
nstenergy                = 500
nstxtcout                = 0
xtc-precision            = 1000
; Neighborsearching and short-range nonbonded interactions
nstlist                  = 10
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.0
; Electrostatics
coulombtype              = PME
rcoulomb                 = 1.0
; van der Waals
vdw-type                 = switch
rvdw-switch              = 0.8
rvdw                     = 0.9
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                  = EnerPres
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; EWALD/PME/PPPM parameters
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
optimize_fft             = no
; Temperature coupling
; tcoupl is implicitly handled by the sd integrator
tc_grps                  = system
tau_t                    = 1.0
ref_t                    = 300
; Pressure coupling is off for NVT
Pcoupl                   = No
tau_p                    = 0.5
compressibility          = 4.5e-05
ref_p                    = 1.0 
; Free energy control stuff
free_energy              = yes
init_lambda              = {lam}
delta_lambda             = 0
foreign_lambda           = {fol}
sc-alpha                 = 0.5
sc-power                 = 1.0
sc-sigma                 = 0.3 
couple-moltype           = {mol}  ; name of moleculetype to decouple
couple-lambda0           = {cp1}      ; only van der Waals interactions
couple-lambda1           = {cp2}     ; turn off everything, in this case only vdW
couple-intramol          = no
nstdhdl                  = 10
; Generate velocities to start
gen_vel                  = yes
gen_temp                 = 300
gen_seed                 = -1
; options for bonds
constraints              = h-bonds  ; we only have C-H bonds here
; Type of constraint algorithm
constraint-algorithm     = lincs
; Do not constrain the starting configuration
continuation             = no
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 12
""".format(lam=lam, fol=fol, mol=mol, cp1=cp1, cp2=cp2)

def npt_mdp(lam="0.0", fol="0.05", mol="Methane", cp1="vdw", cp2="none"):
    return """; Run control
integrator               = sd       ; Langevin dynamics
tinit                    = 0
dt                       = 0.002
nsteps                   = 50000    ; 100 ps
nstcomm                  = 100
; Output control
nstxout                  = 500
nstvout                  = 500
nstfout                  = 0
nstlog                   = 500
nstenergy                = 500
nstxtcout                = 0
xtc-precision            = 1000
; Neighborsearching and short-range nonbonded interactions
nstlist                  = 10
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.0
; Electrostatics
coulombtype              = PME
rcoulomb                 = 1.0
; van der Waals
vdw-type                 = switch
rvdw-switch              = 0.8
rvdw                     = 0.9
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                  = EnerPres
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; EWALD/PME/PPPM parameters
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
optimize_fft             = no
; Temperature coupling
; tcoupl is implicitly handled by the sd integrator
tc_grps                  = system
tau_t                    = 1.0
ref_t                    = 300
; Pressure coupling is on for NPT
Pcoupl                   = Parrinello-Rahman 
tau_p                    = 0.5
compressibility          = 4.5e-05
ref_p                    = 1.0 
; Free energy control stuff
free_energy              = yes
init_lambda              = {lam}
delta_lambda             = 0
foreign_lambda           = {fol}
sc-alpha                 = 0.5
sc-power                 = 1.0
sc-sigma                 = 0.3 
couple-moltype           = {mol}  ; name of moleculetype to decouple
couple-lambda0           = {cp1}      ; only van der Waals interactions
couple-lambda1           = {cp2}     ; turn off everything, in this case only vdW
couple-intramol          = no
nstdhdl                  = 10
; Do not generate velocities
gen_vel                  = no 
; options for bonds
constraints              = h-bonds  ; we only have C-H bonds here
; Type of constraint algorithm
constraint-algorithm     = lincs
; Constrain the starting configuration
; since we are continuing from NVT
continuation             = yes 
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 12
""".format(lam=lam, fol=fol, mol=mol, cp1=cp1, cp2=cp2)

def md_mdp(lam="0.0", fol="0.05", mol="Methane", cp1="vdw", cp2="none"):
    return """; Run control
integrator               = sd       ; Langevin dynamics
tinit                    = 0
dt                       = 0.002
nsteps                   = 2500000  ; 5 ns
nstcomm                  = 100
; Output control
nstxout                  = 500
nstvout                  = 500
nstfout                  = 0
nstlog                   = 500
nstenergy                = 500
nstxtcout                = 0
xtc-precision            = 1000
; Neighborsearching and short-range nonbonded interactions
nstlist                  = 10
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.0
; Electrostatics
coulombtype              = PME
rcoulomb                 = 1.0
; van der Waals
vdw-type                 = switch
rvdw-switch              = 0.8
rvdw                     = 0.9
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                  = EnerPres
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; EWALD/PME/PPPM parameters
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
optimize_fft             = no
; Temperature coupling
; tcoupl is implicitly handled by the sd integrator
tc_grps                  = system
tau_t                    = 1.0
ref_t                    = 300
; Pressure coupling is on for NPT
Pcoupl                   = Parrinello-Rahman 
tau_p                    = 0.5
compressibility          = 4.5e-05
ref_p                    = 1.0 
; Free energy control stuff
free_energy              = yes
init_lambda              = {lam}
delta_lambda             = 0
foreign_lambda           = {fol}
sc-alpha                 = 0.5
sc-power                 = 1.0
sc-sigma                 = 0.3 
couple-moltype           = {mol}  ; name of moleculetype to decouple
couple-lambda0           = {cp1}      ; only van der Waals interactions
couple-lambda1           = {cp2}     ; turn off everything, in this case only vdW
couple-intramol          = no
nstdhdl                  = 10
; Do not generate velocities
gen_vel                  = no 
; options for bonds
constraints              = h-bonds  ; we only have C-H bonds here
; Type of constraint algorithm
constraint-algorithm     = lincs
; Constrain the starting configuration
; since we are continuing from NPT
continuation             = yes 
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 12
""".format(lam=lam, fol=fol, mol=mol, cp1=cp1, cp2=cp2)

class MakeMDP:

    def __init__(self):
        self.integrator = None
        self.tinit = None
        self.dt = None
        self.nsteps = None
        self.nstcomm = None
        self.define = None
        self.emtol = None
        self.emstep = None
        self.niter = None
        self.nbfgscorr = None
        self.nstxout = None
        self.nstvout = None
        self.nstfout = None
        self.nstlog = None
        self.nstenergy = None
        self.nstxtcout = None
        self.xtc_precision = None
        self.nstlist = None
        self.ns_type = None
        self.pbc = None
        self.rlist = None
        self.coulombtype = None
        self.rcoulomb = None
        self.vdw_type = None
        self.rvdw_switch = None
        self.rvdw = None
        self.DispCorr = None
        self.fourierspacing = None
        self.pme_order = None
        self.ewald_rtol = None
        self.epsilon_surface = None
        self.optimize_fft = None
        self.tcoupl = None
        self.tc_grps = None
        self.tau_t = None
        self.ref_t = None
        self.pcoupl = None
        self.tau_p = None
        self.compressibility = None
        self.ref_p = None
        self.free_energy = None
        self.init_lambda = None
        self.delta_lambda = None
        self.foreign_lambda = None
        self.sc_alpha = None
        self.sc_power = None
        self.sc_sigma = None
        self.couple_moltype = None
        self.couple_lambda0 = None
        self.couple_lambda1 = None
        self.couple_intramol = None
        self.nstdhdl = None
        self.gen_vel = None
        self.gen_temp = None
        self.gen_seed = None
        self.constraints = None
        self.constraint_algorithm = None
        self.continuation = None
        self.lincs_order = None
        # comments
        self.comment_coupling1 = None
        self.comment_coupling2 = None
        self.comment_coupling3 = None
        self.comment_velocities1 = None
        self.comment_constraints1 = None
        self.comment_constraints2 = None

    def core(self):
        # reset
        self.__init__()
        # define
        self.ns_type = "grid"
        self.pbc = "xyz"
        self.rlist = "1.0"
        self.coulombtype = "PME"
        self.rcoulomb = "1.0"
        self.vdw_type = "switch"
        self.rvdw_switch = "0.8"
        self.rvdw = "0.9"
        self.DispCorr = "EnerPres"
        self.fourierspacing = "0.12"
        self.pme_order = "6"
        self.ewald_rtol = "1e-06"
        self.epsilon_surface = "0"
        self.optimize_fft = "no"
        self.free_energy = "yes"
        self.init_lambda = "0.0"
        self.delta_lambda = "0"
        self.foreign_lambda = "0.05"
        self.sc_alpha = "0.5"
        self.sc_power = "1.0"
        self.sc_sigma = "0.3 "
        self.couple_moltype = "Methane  ; name of moleculetype to decouple"
        self.couple_lambda0 = "vdw      ; only van der Waals interactions"
        self.couple_lambda1 = "none     ; turn off everything, in this case only vdW"
        self.couple_intramol = "no"
        self.nstdhdl = "10"
        self.constraint_algorithm = "lincs"
        self.lincs_order = "12"

    def enermin_steep(self):
        self.core()
        self.integrator = "steep "
        self.nsteps = "5000"
        self.emtol = "100"
        self.emstep = "0.01"
        self.niter = "20"
        self.nbfgscorr = "10"
        self.nstlog = "1"
        self.nstenergy = "1"
        self.nstlist = "1"
        self.comment_coupling1 = ("Temperature and pressure coupling are" +
            " off during EM")
        self.tcoupl = "no"
        self.pcoupl = "no"
        self.comment_velocities1 = "Generate velocities to start"
        self.gen_vel = "no "
        self.constraints = "h-bonds  ; we only have C-H bonds here"
        self.comment_constraints1 = ("Do not constrain the starting" +
            " configuration")
        self.continuation = "no"

    def enermin_lbfgs(self):
        self.core()
        self.integrator = "l-bfgs"
        self.nsteps = "5000"
        self.define = "-DFLEXIBLE"
        self.emtol = "100"
        self.emstep = "0.01"
        self.niter = "20"
        self.nbfgscorr = "10"
        self.nstlog = "1"
        self.nstenergy = "1"
        self.nstlist = "1"
        self.comment_coupling1 = ("Temperature and pressure coupling are" +
            " off during EM")
        self.tcoupl = "no"
        self.pcoupl = "no"
        self.comment_velocities1 = "Generate velocities to start"
        self.gen_vel = "no "
        self.constraints = "none     ; L-BFGS doesn't work with constraints "
        self.comment_constraints1 = ("Do not constrain the starting" +
            " configuration")
        self.continuation = "no"

    def equilibrate_nvt(self):
        self.core()
        self.integrator = "sd       ; Langevin dynamics"
        self.tinit = "0"
        self.dt = "0.002"
        self.nsteps = "50000    ; 100 ps"
        self.nstcomm = "100"
        self.nstxout = "500"
        self.nstvout = "500"
        self.nstfout = "0"
        self.nstlog = "500"
        self.nstenergy = "500"
        self.nstxtcout = "0"
        self.xtc_precision = "1000"
        self.nstlist = "10"
        self.comment_coupling1 = "Temperature coupling"
        self.comment_coupling2 = ("tcoupl is implicitly handled by the" +
            " sd integrator")
        self.tc_grps = "system"
        self.tau_t = "1.0"
        self.ref_t = "300"
        self.comment_coupling3 = "Pressure coupling is off for NVT"
        self.pcoupl = "no"
        self.tau_p = "0.5"
        self.compressibility = "4.5e-05"
        self.ref_p = "1.0 "
        self.comment_velocities1 = "Generate velocities to start"
        self.gen_vel = "yes"
        self.gen_temp = "300"
        self.gen_seed = "-1"
        self.constraints = "h-bonds  ; we only have C-H bonds here"
        self.comment_constraints1 = ("Do not constrain the starting" +
            " configuration")
        self.continuation = "no"

    def equilibrate_npt(self):
        self.core()
        self.integrator = "sd       ; Langevin dynamics"
        self.tinit = "0"
        self.dt = "0.002"
        self.nsteps = "50000    ; 100 ps"
        self.nstcomm = "100"
        self.nstxout = "500"
        self.nstvout = "500"
        self.nstfout = "0"
        self.nstlog = "500"
        self.nstenergy = "500"
        self.nstxtcout = "0"
        self.xtc_precision = "1000"
        self.nstlist = "10"
        self.comment_coupling1 = "Temperature coupling"
        self.comment_coupling2 = ("tcoupl is implicitly handled by the" +
            " sd integrator")
        self.tc_grps = "system"
        self.tau_t = "1.0"
        self.ref_t = "300"
        self.comment_coupling3 = "Pressure coupling is on for NPT"
        self.pcoupl = "Parrinello-Rahman"
        self.tau_p = "0.5"
        self.compressibility = "4.5e-05"
        self.ref_p = "1.0 "
        self.comment_velocities1 = "Do not generate velocities"
        self.gen_vel = "no "
        self.constraints = "h-bonds  ; we only have C-H bonds here"
        self.comment_constraints1 = ("Constrain the starting configuration")
        self.comment_constraints2 = ("since we are continuing from NVT")
        self.continuation = "yes "

    def production_md(self):
        self.core()
        self.integrator = "sd       ; Langevin dynamics"
        self.tinit = "0"
        self.dt = "0.002"
        self.nsteps = "2500000  ; 5 ns"
        self.nstcomm = "100"
        self.nstxout = "500"
        self.nstvout = "500"
        self.nstfout = "0"
        self.nstlog = "500"
        self.nstenergy = "500"
        self.nstxtcout = "0"
        self.xtc_precision = "1000"
        self.nstlist = "10"
        self.comment_coupling1 = "Temperature coupling"
        self.comment_coupling2 = ("tcoupl is implicitly handled by the" +
            " sd integrator")
        self.tc_grps = "system"
        self.tau_t = "1.0"
        self.ref_t = "300"
        self.comment_coupling3 = "Pressure coupling is on for NPT"
        self.pcoupl = "Parrinello-Rahman"
        self.tau_p = "0.5"
        self.compressibility = "4.5e-05"
        self.ref_p = "1.0 "
        self.comment_velocities1 = "Do not generate velocities"
        self.gen_vel = "no "
        self.constraints = "h-bonds  ; we only have C-H bonds here"
        self.comment_constraints1 = ("Constrain the starting configuration")
        self.comment_constraints2 = ("since we are continuing from NPT")
        self.continuation = "yes "

    def compile(self):
        text = self.run_control()
        text += self.em_criteria()
        text += self.output_control()
        text += self.neighborsearching()
        text += self.electrostatics()
        text += self.vanderWaals()
        text += self.corrections()
        text += self.pme_pppm()
        text += self.coupling()
        text += self.free_energy_control()
        text += self.velocities()
        text += self.bond_constraints()
        return text

    def run_control(self):
        text = """; Run control
"""
        if self.integrator != None:
            text += """integrator               = {0}
""".format(self.integrator)
        if self.tinit != None:
            text += """tinit                    = {0}
""".format(self.tinit)
        if self.dt != None:
            text += """dt                       = {0}
""".format(self.dt)
        if self.nsteps != None:
            text += """nsteps                   = {0}
""".format(self.nsteps)
        if self.nstcomm != None:
            text += """nstcomm                  = {0}
""".format(self.nstcomm)
        if self.define != None:
            text += """define                   = {0}
""".format(self.define)
        return text

    def em_criteria(self):
        if ((self.emtol == None) and (self.emstep == None)
            and (self.niter == None) and (self.nbfgscorr == None)):
            return ""

        text = """; EM criteria and other stuff
"""
        if self.emtol != None:
            text += """emtol                    = {0}
""".format(self.emtol)
        if self.emstep != None:
            text += """emstep                   = {0}
""".format(self.emstep)
        if self.niter != None:
            text += """niter                    = {0}
""".format(self.niter)
        if self.nbfgscorr != None:
            text += """nbfgscorr                = {0}
""".format(self.nbfgscorr)
        return text

    def output_control(self):
        text = """; Output control
"""
        if self.nstxout != None:
            text += """nstxout                  = {0}
""".format(self.nstxout)
        if self.nstvout != None:
            text += """nstvout                  = {0}
""".format(self.nstvout)
        if self.nstfout != None:
            text += """nstfout                  = {0}
""".format(self.nstfout)
        if self.nstlog != None:
            text += """nstlog                   = {0}
""".format(self.nstlog)
        if self.nstenergy != None:
            text += """nstenergy                = {0}
""".format(self.nstenergy)
        if self.nstxtcout != None:
            text += """nstxtcout                = {0}
""".format(self.nstxtcout)
        if self.xtc_precision != None:
            text += """xtc-precision            = {0}
""".format(self.xtc_precision)
        return text

    def neighborsearching(self):
        text = """; Neighborsearching and short-range nonbonded interactions
"""
        if self.nstlist != None:
            text += """nstlist                  = {0}
""".format(self.nstlist)
        if self.ns_type != None:
            text += """ns_type                  = {0}
""".format(self.ns_type)
        if self.pbc != None:
            text += """pbc                      = {0}
""".format(self.pbc)
        if self.rlist != None:
            text += """rlist                    = {0}
""".format(self.rlist)
        return text

    def electrostatics(self):
        text = """; Electrostatics
"""
        if self.coulombtype != None:
            text += """coulombtype              = {0}
""".format(self.coulombtype)
        if self.rcoulomb != None:
            text += """rcoulomb                 = {0}
""".format(self.rcoulomb)
        return text

    def vanderWaals(self):
        text = """; van der Waals
"""
        if self.vdw_type != None:
            text += """vdw-type                 = {0}
""".format(self.vdw_type)
        if self.rvdw_switch != None:
            text += """rvdw-switch              = {0}
""".format(self.rvdw_switch)
        if self.rvdw != None:
            text += """rvdw                     = {0}
""".format(self.rvdw)
        return text

    def corrections(self):
        text = """; Apply long range dispersion corrections for Energy and Pressure
"""
        if self.DispCorr != None:
            text += """DispCorr                  = {0}
""".format(self.DispCorr)
        return text

    def pme_pppm(self):
        text = """; Spacing for the PME/PPPM FFT grid
"""
        if self.fourierspacing != None:
            text += """fourierspacing           = {0}
""".format(self.fourierspacing)
        text += """; EWALD/PME/PPPM parameters
"""
        if self.pme_order != None:
            text += """pme_order                = {0}
""".format(self.pme_order)
        if self.ewald_rtol != None:
            text += """ewald_rtol               = {0}
""".format(self.ewald_rtol)
        if self.epsilon_surface != None:
            text += """epsilon_surface          = {0}
""".format(self.epsilon_surface)
        if self.optimize_fft != None:
            text += """optimize_fft             = {0}
""".format(self.optimize_fft)
        return text

    def coupling(self):
        text = ""
        if self.comment_coupling1 != None:
            text += """; {0}
""".format(self.comment_coupling1)
        if self.comment_coupling2 != None:
            text += """; {0}
""".format(self.comment_coupling2)
        if self.tcoupl != None:
            text += """tcoupl                   = {0}
""".format(self.tcoupl)
        if self.tc_grps != None:
            text += """tc_grps                  = {0}
""".format(self.tc_grps)
        if self.tau_t != None:
            text += """tau_t                    = {0}
""".format(self.tau_t)
        if self.ref_t != None:
            text += """ref_t                    = {0}
""".format(self.ref_t)

        if self.comment_coupling3 != None:
            text += """; {0}
""".format(self.comment_coupling3)
        if self.pcoupl != None:
            text += """pcoupl                   = {0}
""".format(self.pcoupl)
        if self.tau_p != None:
            text += """tau_p                    = {0}
""".format(self.tau_p)
        if self.compressibility != None:
            text += """compressibility          = {0}
""".format(self.compressibility)
        if self.ref_p != None:
            text += """ref_p                    = {0}
""".format(self.ref_p)
        return text

    def free_energy_control(self):
        text = """; Free energy control stuff
"""
        if self.free_energy != None:
            text += """free_energy              = {0}
""".format(self.free_energy)
        if self.init_lambda != None:
            text += """init_lambda              = {0}
""".format(self.init_lambda)
        if self.delta_lambda != None:
            text += """delta_lambda             = {0}
""".format(self.delta_lambda)
        if self.foreign_lambda != None:
            text += """foreign_lambda           = {0}
""".format(self.foreign_lambda)
        if self.sc_alpha != None:
            text += """sc-alpha                 = {0}
""".format(self.sc_alpha)
        if self.sc_power != None:
            text += """sc-power                 = {0}
""".format(self.sc_power)
        if self.sc_sigma != None:
            text += """sc-sigma                 = {0}
""".format(self.sc_sigma)
        if self.couple_moltype != None:
            text += """couple-moltype           = {0}
""".format(self.couple_moltype)
        if self.couple_lambda0 != None:
            text += """couple-lambda0           = {0}
""".format(self.couple_lambda0)
        if self.couple_lambda1 != None:
            text += """couple-lambda1           = {0}
""".format(self.couple_lambda1)
        if self.couple_intramol != None:
            text += """couple-intramol          = {0}
""".format(self.couple_intramol)
        if self.nstdhdl != None:
            text += """nstdhdl                  = {0}
""".format(self.nstdhdl)
        return text

    def velocities(self):
        text = ""
        if self.comment_velocities1 != None:
            text += """; {0}
""".format(self.comment_velocities1)
        if self.gen_vel != None:
            text += """gen_vel                  = {0}
""".format(self.gen_vel)
        if self.gen_temp != None:
            text += """gen_temp                 = {0}
""".format(self.gen_temp)
        if self.gen_seed != None:
            text += """gen_seed                 = {0}
""".format(self.gen_seed)
        return text

    def bond_constraints(self):
        text = """; options for bonds
"""
        if self.constraints != None:
            text += """constraints              = {0}
""".format(self.constraints)
        text += """; Type of constraint algorithm
"""
        if self.constraint_algorithm != None:
            text += """constraint-algorithm     = {0}
""".format(self.constraint_algorithm)

        if self.comment_constraints1 != None:
            text += """; {0}
""".format(self.comment_constraints1)
        if self.comment_constraints2 != None:
            text += """; {0}
""".format(self.comment_constraints2)
        if self.continuation != None:
            text += """continuation             = {0}
""".format(self.continuation)
        text += """; Highest order in the expansion of the constraint coupling matrix
"""
        if self.lincs_order != None:
            text += """lincs-order              = {0}
""".format(self.lincs_order)
        return text

class MakeSLURM:

    def __init__(self, job_name, suffix, workdir, ntasks=4):
        self.double_precision = True
        self.use_mpi = False
        self.gromacs5 = False
        self.calc_wt = True
        self.fields_general = dict()
        self.fields_header = dict()
        self.fields_general['job-name'] = job_name
        self.fields_general['suffix'] = suffix
        self.fields_general['ntasks'] = ntasks
        self.fields_general['workdir'] = os.path.realpath(workdir)
        self.fields_header['outputfile'] = os.path.realpath(
            os.path.join(workdir, job_name))

    def walltime(self):
        return '{0}-{1:02}:{2:02}:{3:02}'.format(0, 12, 0, 0)

    def get_header(self):
        fields = dict()
        fields.update(self.fields_general)
        fields.update(self.fields_header)
        fields["nnodes"] = 1
        fields["partition"] = "serial"
        if self.use_mpi:
            # 20 cores physically exist on a Rivanna node
            fields["nnodes"] = int(math.ceil(fields["ntasks"] / 20.0))
        if fields["nnodes"] == 0:
            fields["nnodes"] = 1
        elif fields["nnodes"] > 1:
            fields["partition"] = "parallel"
        if self.calc_wt:
            fields["queue-time"] = self.walltime()
        self.header = """#!/bin/sh
#SBATCH --job-name={job-name}{suffix}
#SBATCH --partition={partition}
#SBATCH --nodes={nnodes}
#SBATCH --ntasks={ntasks}
#SBATCH --time={queue-time}
#SBATCH --signal=15 --comment="15 = SIGTERM"
#SBATCH --output={outputfile}.out
#SBATCH --workdir={workdir}
#SBATCH --checkpoint=06:00:00
#SBATCH --checkpoint-dir=/scratch/jtp4kc/checkpoints
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-user=jtp4kc@virginia.edu
#SBATCH --comment="#SBATCH --array-0,1,3"

MODULEPATH=/h3/n1/shirtsgroup/modules:$MODULEPATH
MODULEPATH=$HOME/modules:$MODULEPATH

module load jtp4kc
module load gromacs-jtp4kc
""".format(**fields)
        return self.header

    def get_ligandextract(self):
        text = ("python ~/git/expdens/gromod.py -n t41meth-in.gro" +
            " -o ligand.gro -i TMP -v -s -t 1methylpyrrole -m center\n")
        text += ("python ~/git/expdens/gromod.py -n t41meth-in.gro" +
            " -o solvent.gro -i SOL -v -t water\n")
        return text

    def get_boxgen(self):
        return ("genbox -cp ligand.gro -cs solvent.gro" +
            " -ci solvent.gro -p {top} -o {gro}" +
            " -nmol 6500 -maxsol 10000\n\n")

    def get_steep(self):
        return """
#################################
# ENERGY MINIMIZATION 1: STEEP  #
#################################
echo "Starting minimization for lambda = $LAMBDA..." 

# Iterative calls to grompp and mdrun to run the simulations

grompp{_d} -f em_steep.mdp -c {gro} -p {top} -o mins.tpr
count=$[ 0 ]
while [ ! -f mins.tpr ]
do
    grompp{_d} -f em_steep.mdp -c {gro} -p {top} -o mins.tpr
    count=$[ $count + 1 ]
    if [ "$count" -gt 5 ]
    then
        echo "max tries to make mins.tpr encountered"
        break
    fi
done
count=$[ 0 ]
mdrun{_d} -nt 4 -deffnm mins
while [ ! -f mins.gro ]
do
    mdrun{_d} -nt 4 -deffnm mins
    count=$[ $count + 1 ]
    if [ "$count" -gt 5 ]
    then
        echo "max tries to make mins.gro encountered"
        break
    fi
done
"""

    def get_lbfgs(self):
        return """
#################################
# ENERGY MINIMIZATION 2: L-BFGS #
#################################
grompp{_d} -f em_l-bfgs.mdp -c mins.gro -p {top} -o minl.tpr
count=$[ 0 ]
while [ ! -f minl.tpr ]
do
    grompp{_d} -f em_l-bfgs.mdp -c mins.gro -p {top} -o minl.tpr
    count=$[ $count + 1 ]
    if [ "$count" -gt 5 ]
    then
        echo "max tries to make minl.tpr encountered"
        break
    fi
done
# Run L-BFGS in serial (cannot be run in parallel)
mdrun{_d} -nt 1 -deffnm minl
count=$[ 0 ]
while [ ! -f minl.gro ]
do
    mdrun{_d} -nt 1 -deffnm minl
    count=$[ $count + 1 ]
    if [ "$count" -gt 5 ]
    then
        echo "max tries to make minl.gro encountered"
        break
    fi
done
"""

    def get_nvt(self, lbfgs=True):
        name = "mins.gro"
        if lbfgs:
            name = "minl.gro"
        return """
#####################
# NVT EQUILIBRATION #
#####################
echo "Starting constant volume equilibration..."

grompp{_d} -f nvt.mdp -c """ + name + """ -p {top} -o nvt.tpr
count=$[ 0 ]
while [ ! -f nvt.tpr ]
do
    grompp{_d} -f nvt.mdp -c """ + name + """ -p {top} -o nvt.tpr
    count=$[ $count + 1 ]
    if [ "$count" -gt 5 ]
    then
        echo "max tries to make nvt.tpr encountered"
        break
    fi
done
mdrun{_d} -nt 4 -deffnm nvt
count=$[ 0 ]
while [ ! -f nvt.gro ]
do
    mdrun{_d} -nt 4 -deffnm nvt
    count=$[ $count + 1 ]
    if [ "$count" -gt 5 ]
    then
        echo "max tries to make nvt.gro encountered"
        break
    fi
done

echo "Constant volume equilibration complete."
"""

    def get_npt(self):
        return """
#####################
# NPT EQUILIBRATION #
#####################
echo "Starting constant pressure equilibration..."

grompp{_d} -f npt.mdp -c nvt.gro -p {top} -t nvt.cpt -o npt.tpr
count=$[ 0 ]
while [ ! -f npt.tpr ]
do
    grompp{_d} -f npt.mdp -c nvt.gro -p {top} -t nvt.cpt -o npt.tpr
    count=$[ $count + 1 ]
    if [ "$count" -gt 5 ]
    then
        echo "max tries to make npt.tpr encountered"
        break
    fi
done
mdrun{_d} -nt 4 -deffnm npt
count=$[ 0 ]
while [ ! -f npt.gro ]
do
    mdrun{_d} -nt 4 -deffnm npt
    count=$[ $count + 1 ]
    if [ "$count" -gt 5 ]
    then
        echo "max tries to make npt.gro encountered"
        break
    fi
done

echo "Constant pressure equilibration complete."
"""

    def get_md(self):
        return """
#################
# PRODUCTION MD #
#################
echo "Starting production MD simulation..."

grompp{_d} -f md.mdp -c npt.gro -p {top} -t npt.cpt -o md.tpr
count=$[ 0 ]
while [ ! -f md.tpr ]
do
    grompp{_d} -f md.mdp -c npt.gro -p {top} -t npt.cpt -o md.tpr
    count=$[ $count + 1 ]
    if [ "$count" -gt 5 ]
    then
        echo "max tries to make md.tpr encountered"
        break
    fi
done
mdrun{_d} -nt 4 -deffnm md
count=$[ 0 ]
while [ ! -f md.gro ]
do
    mdrun{_d} -nt 4 -deffnm md
    count=$[ $count + 1 ]
    if [ "$count" -gt 5 ]
    then
        echo "max tries to make md.gro encountered"
        break
    fi
done

echo "Production MD complete."
"""

    def get_text(self, use_lbfgs=True, use_boxgen=False):
        text = self.get_header()
        if use_boxgen:
            text += self.get_ligandextract()
            text += self.get_boxgen()
        text += self.get_steep()
        if use_lbfgs:
            text += self.get_lbfgs()
        text += """
echo "Minimization complete."
""" + self.get_nvt(use_lbfgs) + self.get_npt() + self.get_md() + """
# End
echo "Ending. Job completed for lambda = {lam}"
"""
        return text

    def compile(self, gro, top, lam, use_lbfgs=True, use_boxgen=False):
        _d = ""
        if self.double_precision:
            _d = "_d"
        return self.get_text(use_lbfgs, use_boxgen).format(_d=_d,
                    gro=gro, top=top, lam=lam)

def format_lam(lambda_):
    lam = "{0:0.2f}".format(lambda_)
    while (lam.endswith("0")):
        lam = lam[0:-1]
        if lam.endswith("."):
            lam = lam[0:-1]
            break
    return lam

def format_fol(foreign):
    fol = ""
    for frn in foreign:
        if frn >= 0 and frn <= 1:
            fol += format_lam(frn) + " "
    return fol

def output(name, string):
    out = open(name, "w")
    out.write(string)
    out.close()

def launch_meth():
    jobname = "lambda"
    topfile = "alab.top"
    grofile = "methane_water.gro"
    dl = 0.05
    schedule = [dl * x for x in range(0, int(1.0 / dl) + 1)]
    for lambda_ in schedule:
        foreign = [lambda_ - dl, lambda_ + dl]
        lam = format_lam(lambda_)
        fol = format_fol(foreign)
        folder = jobname + "_" + lam
        if not os.path.exists(folder):
            os.mkdir(folder)
        os.chdir(folder)
        slurm = MakeSLURM(jobname, "_" + lam, ".")
        outtext = slurm.compile(os.path.join("..", grofile), os.path.join("..",
            topfile), lam)
        output("em_steep.mdp", em_steep_mdp(lam, fol))
        output("em_l-bfgs.mdp", em_lbfgs_mdp(lam, fol))
        output("nvt.mdp", nvt_mdp(lam, fol))
        output("npt.mdp", npt_mdp(lam, fol))
        output("md.mdp", md_mdp(lam, fol))
        output("job.slurm", outtext)
        os.system("sbatch job.slurm")
        os.chdir("..")

def launch_meth2():
    jobname = "lambda"
    topfile = "alab.top"
    grofile = "methane_water.gro"
    dl = 0.05
    schedule = [dl * x for x in range(0, int(1.0 / dl) + 1)]
    for lambda_ in schedule:
        foreign = [lambda_ - dl, lambda_ + dl]
        lam = format_lam(lambda_)
        fol = format_fol(foreign)
        folder = jobname + "_" + lam
        if not os.path.exists(folder):
            os.mkdir(folder)
        os.chdir(folder)
        slurm = MakeSLURM(jobname, "_" + lam, ".")
        outtext = slurm.compile(os.path.join("..", grofile), os.path.join("..",
            topfile), lam)
        mdpgen = MakeMDP()
        mdpgen.enermin_steep()
        mdpgen.init_lambda = lam
        mdpgen.foreign_lambda = fol
        output("em_steep.mdp", mdpgen.compile())
        mdpgen.enermin_lbfgs()
        mdpgen.init_lambda = lam
        mdpgen.foreign_lambda = fol
        output("em_l-bfgs.mdp", mdpgen.compile())
        mdpgen.equilibrate_nvt()
        mdpgen.init_lambda = lam
        mdpgen.foreign_lambda = fol
        output("nvt.mdp", mdpgen.compile())
        mdpgen.equilibrate_npt()
        mdpgen.init_lambda = lam
        mdpgen.foreign_lambda = fol
        output("npt.mdp", mdpgen.compile())
        mdpgen.production_md()
        mdpgen.init_lambda = lam
        mdpgen.foreign_lambda = fol
        output("md.mdp", mdpgen.compile())
        output("job.slurm", outtext)
        os.system("sbatch job.slurm")
        os.chdir("..")

def launch2():
    jobname = "1methsolv"
    topfile = "1meth.top"
    grofile = "1meth.gro"
    mol = "TMP"

    os.system("python ~/git/expdens/gromod.py -n t41meth-in.gro" +
        " -o ligand.gro -i TMP -v -s -t 1methylpyrrole -m center")
    os.system("python ~/git/expdens/gromod.py -n t41meth-in.gro" +
            " -o solvent.gro -i SOL -v -t water")
    os.system("genbox -cp ligand.gro -cs solvent.gro" +
            " -ci solvent.gro -p " + topfile + " -o " + grofile +
            " -nmol 6500 -maxsol 10000")

    coul1 = "vdw-q"
    coul2 = "vdw"
    dc = 0.25
    dl = 0.05
    coulomb = [dc * x for x in range(0, int(1.0 / dc) + 1)]
    schedule = [dl * x for x in range(0, int(1.0 / dl) + 1)]
    for lambda_ in coulomb:
        foreign = [lambda_ - dc, lambda_ + dc]
        lam = format_lam(lambda_)
        fol = format_fol(foreign)
        folder = jobname + "C" + lam
        if not os.path.exists(folder):
            os.mkdir(folder)
        os.chdir(folder)
        slurm = MakeSLURM(jobname, "C" + lam, ".")
        slurm.double_precision = True
        outtext = slurm.compile(os.path.join("..", grofile),
            os.path.join("..", topfile), lam, use_lbfgs=False)
        output("em_steep.mdp", em_steep_mdp(lam, fol, mol, coul1, coul2))
        # output("em_l-bfgs.mdp", em_lbfgs_mdp(lam, fol, mol, coul1, coul2))
        output("nvt.mdp", nvt_mdp(lam, fol, mol, coul1, coul2))
        output("npt.mdp", npt_mdp(lam, fol, mol, coul1, coul2))
        output("md.mdp", md_mdp(lam, fol, mol, coul1, coul2))
        output("job.slurm", outtext)
        os.system("sbatch job.slurm")
        os.chdir("..")
    for lambda_ in schedule:
        foreign = [lambda_ - dl, lambda_ + dl]
        lam = format_lam(lambda_)
        fol = format_fol(foreign)
        folder = jobname + "_" + lam
        if not os.path.exists(folder):
            os.mkdir(folder)
        os.chdir(folder)
        slurm = MakeSLURM(jobname, "_" + lam, ".")
        outtext = slurm.compile(os.path.join("..", grofile),
            os.path.join("..", topfile), lam, use_lbfgs=False)
        output("em_steep.mdp", em_steep_mdp(lam, fol, mol))
        # output("em_l-bfgs.mdp", em_lbfgs_mdp(lam, fol, mol))
        output("nvt.mdp", nvt_mdp(lam, fol, mol))
        output("npt.mdp", npt_mdp(lam, fol, mol))
        output("md.mdp", md_mdp(lam, fol, mol))
        output("job.slurm", outtext)
        os.system("sbatch job.slurm")
        os.chdir("..")

def do_set(mdpgen, lam, fol, mol, coul1, coul2):
    mdpgen.init_lambda = lam
    mdpgen.foreign_lambda = fol
    mdpgen.couple_moltype = mol
    mdpgen.couple_lambda0 = coul1
    mdpgen.couple_lambda1 = coul2

def array(mdpgen, lam, fol, mol, coul1="vdw", coul2="none", lbfgs=False):
    mdpgen.enermin_steep()
    do_set(mdpgen, lam, fol, mol, coul1, coul2)
    output("em_steep.mdp", mdpgen.compile())

    if lbfgs:
        mdpgen.enermin_lbfgs()
        do_set(mdpgen, lam, fol, mol, coul1, coul2)
        output("em_l-bfgs.mdp", mdpgen.compile())

    mdpgen.equilibrate_nvt()
    do_set(mdpgen, lam, fol, mol, coul1, coul2)
    output("nvt.mdp", mdpgen.compile())

    mdpgen.equilibrate_npt()
    do_set(mdpgen, lam, fol, mol, coul1, coul2)
    output("npt.mdp", mdpgen.compile())

    mdpgen.production_md()
    do_set(mdpgen, lam, fol, mol, coul1, coul2)
    output("md.mdp", mdpgen.compile())

def launch():
    jobname = "1methsolv"
    topfile = "1meth.top"
    grofile = "1meth.gro"
    mol = "TMP"

    mdpgen = MakeMDP()

    os.system("python ~/git/expdens/gromod.py -n t41meth-in.gro" +
        " -o ligand.gro -i TMP -v -s -t 1methylpyrrole -m center")
    os.system("python ~/git/expdens/gromod.py -n t41meth-in.gro" +
            " -o solvent.gro -i SOL -v -t water")
    os.system("genbox -cp ligand.gro -cs solvent.gro" +
            " -ci solvent.gro -p " + topfile + " -o " + grofile +
            " -nmol 6500 -maxsol 10000")

    coul1 = "vdw-q"
    coul2 = "vdw"
    dc = 0.25
    dl = 0.05
    coulomb = [dc * x for x in range(0, int(1.0 / dc) + 1)]
    schedule = [dl * x for x in range(0, int(1.0 / dl) + 1)]
    for lambda_ in coulomb:
        foreign = [lambda_ - dc, lambda_ + dc]
        lam = format_lam(lambda_)
        fol = format_fol(foreign)
        folder = jobname + "C" + lam
        if not os.path.exists(folder):
            os.mkdir(folder)
        os.chdir(folder)
        slurm = MakeSLURM(jobname, "C" + lam, ".")
        slurm.double_precision = True
        outtext = slurm.compile(os.path.join("..", grofile),
            os.path.join("..", topfile), lam, use_lbfgs=False)
        array(mdpgen, lam, fol, mol, coul1, coul2)
        output("job.slurm", outtext)
        os.system("sbatch job.slurm")
        os.chdir("..")
    for lambda_ in schedule:
        foreign = [lambda_ - dl, lambda_ + dl]
        lam = format_lam(lambda_)
        fol = format_fol(foreign)
        folder = jobname + "_" + lam
        if not os.path.exists(folder):
            os.mkdir(folder)
        os.chdir(folder)
        slurm = MakeSLURM(jobname, "_" + lam, ".")
        outtext = slurm.compile(os.path.join("..", grofile),
            os.path.join("..", topfile), lam, use_lbfgs=False)
        array(mdpgen, lam, fol, mol)
        output("job.slurm", outtext)
        os.system("sbatch job.slurm")
        os.chdir("..")

def testmdp():
    output("0steep.mdp", em_steep_MDP())
    output("0lbfgs.mdp", em_lbfgs_MDP())
    output("0nvt.mdp", nvt_MDP())
    output("0npt.mdp", npt_MDP())
    output("0md.mdp", md_MDP())

    output("1steep.mdp", em_steep_mdp())
    output("1lbfgs.mdp", em_lbfgs_mdp())
    output("1nvt.mdp", nvt_mdp())
    output("1npt.mdp", npt_mdp())
    output("1md.mdp", md_mdp())

    maker = MakeMDP()
    maker.enermin_steep()
    output("2steep.mdp", maker.compile())
    maker.enermin_lbfgs()
    output("2lbfgs.mdp", maker.compile())
    maker.equilibrate_nvt()
    output("2nvt.mdp", maker.compile())
    maker.equilibrate_npt()
    output("2npt.mdp", maker.compile())
    maker.production_md()
    output("2md.mdp", maker.compile())

    os.system("diff 1steep.mdp 2steep.mdp &> diff_steep.txt")
    os.system("diff 1lbfgs.mdp 2lbfgs.mdp &> diff_lbfgs.txt")
    os.system("diff 1nvt.mdp 2nvt.mdp &> diff_nvt.txt")
    os.system("diff 1npt.mdp 2npt.mdp &> diff_npt.txt")
    os.system("diff 1md.mdp 2md.mdp &> diff_md.txt")

def testplain():
    folder_old = "5382947-old"
    folder_new = "5382947-new"
    os.system("rm -r " + folder_old)
    os.mkdir(folder_old)
    os.chdir(folder_old)
    launch_meth()
    os.chdir("..")
    os.system("rm -r " + folder_new)
    os.mkdir(folder_new)
    os.chdir(folder_new)
    launch_meth2()
    os.chdir("..")
    # launch2()
    os.system("diff " + folder_old + "/lambda_0.3/em_steep.mdp" +
        " " + folder_new + "/lambda_0.3/em_steep.mdp &> diff_steep.txt")
    os.system("diff " + folder_old + "/lambda_0.3/em_l-bfgs.mdp" +
        " " + folder_new + "/lambda_0.3/em_l-bfgs.mdp &> diff_lbfgs.txt")
    os.system("diff " + folder_old + "/lambda_0.3/nvt.mdp " + folder_new +
        "/lambda_0.3/nvt.mdp &> diff_nvt.txt")
    os.system("diff " + folder_old + "/lambda_0.3/npt.mdp " + folder_new +
        "/lambda_0.3/npt.mdp &> diff_npt.txt")
    os.system("diff " + folder_old + "/lambda_0.3/md.mdp " + folder_new +
        "/lambda_0.3/md.mdp &> diff_md.txt")

def main(argv=None):  # IGNORE:C0111
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

  Created by user_name on %s.
  Copyright 2015 organization_name. All rights reserved.

  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license)  # , formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-r", "--recursive", dest="recurse", action="store_true", help="recurse into subfolders [default: %(default)s]")
        parser.add_argument("-v", "--verbose", dest="verbose", action="count", help="set verbosity level [default: %(default)s]")
        parser.add_argument("-i", "--include", dest="include", help="only include paths matching this regex pattern. Note: exclude is given preference over include. [default: %(default)s]", metavar="RE")
        parser.add_argument("-e", "--exclude", dest="exclude", help="exclude paths matching this regex pattern. [default: %(default)s]", metavar="RE")
        parser.add_argument('-V', '--version', action='version', version=program_version_message)
        parser.add_argument(dest="paths", help="paths to folder(s) with source file(s) [default: %(default)s]", metavar="path", nargs='*')

        # Process arguments
        args = parser.parse_args()

        paths = args.paths
        verbose = args.verbose
        recurse = args.recurse
        inpat = args.include
        expat = args.exclude

        if verbose > 0:
            print("Verbose mode on")
            if recurse:
                print("Recursive mode on")
            else:
                print("Recursive mode off")

        if inpat and expat and inpat == expat:
            raise CLIError("include and exclude pattern are equal! Nothing will be processed.")

        for inpath in paths:
            ### do something with inpath ###
            print(inpath)

        launch()

        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
#     except Exception, e:
#         if DEBUG or TESTRUN:
#             raise(e)
#         indent = len(program_name) * " "
#         sys.stderr.write(program_name + ": " + repr(e) + "\n")
#         sys.stderr.write(indent + "  for help use --help")
#         return 2

if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-h")
        sys.argv.append("-v")
        sys.argv.append("-r")
    if TESTRUN:
        import doctest
        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'solvate_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())
