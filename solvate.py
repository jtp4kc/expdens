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

import mdp_template
from argparse import ArgumentParser
from mdp_template import MDP
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

class MakeMDP(mdp_template.MDP):

    def core(self):
        self.reset()
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
        self.nstxout = "0"  # "500"
        self.nstvout = "0"  # "500"
        self.nstfout = "0"
        self.nstlog = "500"
        self.nstenergy = "500"
        self.nstxtcout = "500"  # "0"
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

    def modify_1meth(self, double=True, pcouple=False):
        if double:
            tau_t = 0.5
            tau_p = 20.0
            shake_tol = 1e-12
        else:
            tau_t = 1.0
            tau_p = 5.0
            shake_tol = 5e-6

        self.include = "-I/nv/blue/jtp4kc/gromacs/itp/"
        self.integrator = "md-vv"
        self.comm_mode = "Linear"
        self.nstcomm = "1"
        self.nsttcouple = "1"  #
        self.nstpcouple = "1"  #
        # neighbors
        self.rlist = "1.2"
        # output control
        self.nstcalcenergy = "1"  #
        # interactions
        self.cutoff_scheme = "group"  #
        self.coulombtype = "PME"  #
        self.coulomb_modifier = "Potential-Switch"  #
        self.rcoulomb_switch = "0.88"  #
        self.rcoulomb = "0.9"
        self.vdw_type = "Cut-off"
        self.vdw_modifier = "Potential-switch"  #
        self.rvdw_switch = "0.85"
        self.rvdw = "0.9"
        self.DispCorr = "AllEnerPres"
        self.fourierspacing = "0.12"
        self.fourier_nx = "0"  #
        self.fourier_ny = "0"  #
        self.fourier_nz = "0"  #
        self.pme_order = "4"
        self.ewald_rtol = "1e-05"  #
        self.ewald_geometry = "3d"  #
        # bonds
        self.constraint_algorithm = "shake"
        self.shake_tol = str(shake_tol)  #
        # coupling
        self.tc_grps = "System"
        self.tcoupl = "Nose-Hoover"  # don't use MTTK with v-rescale
        self.tau_t = str(tau_t)
        self.ref_t = "300.0"
        if pcouple:
            self.pcoupl = "MTTK"
            self.pcoupltype = "isotropic"  #
            self.tau_p = str(tau_p)
            self.compressibility = "4.5e-5 ; 1/bar"
            self.ref_p = "1.0 ; bar"
        else:
            self.pcoupl = "no"
        # free energy
        self.dhdl_print_energy = "yes"

class MakeSLURM:

    def __init__(self, job_name, suffix, workdir, ntasks=1):
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
        return '{0}-{1:02}:{2:02}:{3:02}'.format(3, 12, 0, 0)

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
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jtp4kc@virginia.edu
#SBATCH --comment="#SBATCH --array-0,1,3"

MODULEPATH=/h3/n1/shirtsgroup/modules:$MODULEPATH
MODULEPATH=$HOME/modules:$MODULEPATH

module load jtp4kc
module load gromacs-jtp4kc

export GMX_SUPPRESS_DUMP=1 #prevent step file output

sleep 1
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
        n = str(self.fields_general["ntasks"])
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
mdrun{_d} -nt """ + n + """ -deffnm mins
while [ ! -f mins.gro ]
do
    mdrun{_d} -nt """ + n + """ -deffnm mins
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
        n = str(self.fields_general["ntasks"])
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
mdrun{_d} -nt """ + n + """ -deffnm nvt
count=$[ 0 ]
while [ ! -f nvt.gro ]
do
    mdrun{_d} -nt """ + n + """ -deffnm nvt
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
        n = str(self.fields_general["ntasks"])
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
mdrun{_d} -nt """ + n + """ -deffnm npt
count=$[ 0 ]
while [ ! -f npt.gro ]
do
    mdrun{_d} -nt """ + n + """ -deffnm npt
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
        n = str(self.fields_general["ntasks"])
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
mdrun{_d} -nt """ + n + """ -deffnm md
count=$[ 0 ]
while [ ! -f md.gro ]
do
    mdrun{_d} -nt """ + n + """ -deffnm md
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

def genscript(double=False, gromacs5=False, volume_control=False,
        fixed_state=False, initial_weights=False, use_gibbs=True,
        use_metro=True):
    mdp = MakeMDP()

    if double:
        tau_t = 0.5
        tau_p = 20.0
        shake_tol = 1e-12
    else:
        tau_t = 1.0
        tau_p = 5.0
        shake_tol = 5e-6

    # params
    mdp.integrator = "md-vv"
    mdp.tinit = "0"
    mdp.dt = "0.002"
    mdp.nsteps = "2500000   ; 5 ns"
    mdp.comm_mode = "Linear"
    mdp.nstcomm = "1"
    mdp.nsttcouple = "1"  #
    mdp.nstpcouple = "1"  #
    # neighbors
    mdp.nstlist = "10"
    mdp.ns_type = "grid"
    mdp.pbc = "xyz"
    mdp.rlist = "1.2"
    # output control
    mdp.nstxout = "0"
    mdp.nstvout = "0"
    mdp.nstfout = "0"
    mdp.nstlog = "500"
    mdp.nstcalcenergy = "1"  #
    mdp.nstenergy = "500"
    if gromacs5:
        mdp.nstxout_compressed = "500"  #
        mdp.compressed_x_precision = "1000"  #
    else:
        mdp.nstxtcout = "500"
        mdp.xtc_precision = "1000"
    # interactions
    mdp.cutoff_scheme = "group"  #
    mdp.coulombtype = "PME"  #
    mdp.coulomb_modifier = "Potential-Switch"  #
    mdp.rcoulomb_switch = "0.88"  #
    mdp.rcoulomb = "0.9"
    mdp.vdw_type = "Cut-off"
    mdp.vdw_modifier = "Potential-switch"  #
    mdp.rvdw_switch = "0.85"
    mdp.rvdw = "0.9"
    mdp.DispCorr = "AllEnerPres"
    mdp.fourierspacing = "0.12"
    mdp.fourier_nx = "0"  #
    mdp.fourier_ny = "0"  #
    mdp.fourier_nz = "0"  #
    mdp.pme_order = "4"
    mdp.ewald_rtol = "1e-05"  #
    mdp.ewald_geometry = "3d"  #
    # bonds
    mdp.constraints = "hbonds ; constrain bonds to hydrogen"
    mdp.constraint_algorithm = "shake"
    mdp.shake_tol = str(shake_tol)  #
    # velocities
    mdp.gen_vel = "yes"
    mdp.gen_temp = "300.0 ; K"
    mdp.gen_seed = "10200"
    # coupling
    mdp.tc_grps = "System"
    mdp.tcoupl = "v-rescale"
    mdp.tau_t = str(tau_t)
    mdp.ref_t = "300.0"
    if volume_control:
        mdp.pcoupl = "no"
    else:
        mdp.pcoupl = "MTTK"
        mdp.pcoupltype = "isotropic"  #
        mdp.tau_p = str(tau_p)
        mdp.compressibility = "4.5e-5 ; 1/bar"
        mdp.ref_p = "1.0 ; bar"
    # free energy
    mdp.sc_alpha = "0.5"
    mdp.couple_moltype = "TMP"
    mdp.couple_lambda0 = "vdw-q"
    mdp.couple_lambda1 = "none"
    mdp.couple_intramol = "no"
    mdp.fep_lambdas = ("[0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0" +
        " 0.0 0.0 0.0 0.0 0.0]")  #
    mdp.coul_lambdas = ("[0.0 0.3 0.5 1.0 1.0 1.0 1.0 1.0 1.0" +
        " 1.0 1.0 1.0 1.0 1.0]")  #
    mdp.vdw_lambdas = ("[0.0 0.0 0.0 0.0 0.1 0.2 0.3 0.4 0.5" +
        " 0.6 0.7 0.8 0.9 1.0]")  #
    mdp.symmetrized_transition_matrix = "yes"  #
    mdp.nst_transition_matrix = "100000"  #
    mdp.nstdhdl = "500"
    if gromacs5:
        mdp.dhdl_print_energy = "total"  #
    else:
        mdp.dhdl_print_energy = "yes"  #
    if fixed_state:
        mdp.free_energy = "yes"
        mdp.init_lambda = "0"
    else:
        mdp.free_energy = "expanded"
        mdp.init_lambda_state = "0"  #
    # expanded ensemble
    mdp.nstexpanded = "20"  #
    mdp.lmc_stats = "wang-landau"  #
    mdp.lmc_weights_equil = "wl-delta"  #
    mdp.lmc_seed = "10200"  #
    if use_gibbs and use_metro:
        mdp.lmc_move = "metropolized-gibbs"  #
        mdp.lmc_gibbsdelta = "-1"  #
    elif use_gibbs:
        mdp.lmc_move = "gibbs"
        mdp.lmc_gibbsdelta = "-1"
    elif use_metro:
        mdp.lmc_move = "metropolis"
    else:
        mdp.lmc_move = "no"
    if initial_weights:
        mdp.init_lambda_weights = "[0.0000]"  #
    else:
        mdp.weight_equil_wl_delta = "0.001"  #
    # monte carlo
    mdp.wl_scale = "0.8"
    mdp.wl_ratio = "0.8"
    mdp.init_wl_delta = "0.50"
    mdp.wl_oneovert = "yes"

    return mdp

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

def dos_set(mdpgen, lam, mol, coup1, coup2, fep, vdw, coul):
    mdpgen.delta_lambda = None
    mdpgen.foreign_lambda = None
    mdpgen.init_lambda = None
    mdpgen.init_lambda_state = lam
    mdpgen.calc_lambda_neighbors = "-1"
    mdpgen.couple_moltype = mol
    mdpgen.couple_lambda0 = coup1
    mdpgen.couple_lambda1 = coup2
    mdpgen.fep_lambdas = "["
    mdpgen.vdw_lambdas = "["
    mdpgen.coul_lambdas = "["
    for (f, v, c) in zip(fep, vdw, coul):
        mdpgen.fep_lambdas += "{0:0.2f} ".format(f)
        mdpgen.vdw_lambdas += "{0:0.2f} ".format(v)
        mdpgen.coul_lambdas += "{0:0.2f} ".format(c)
    mdpgen.fep_lambdas += "]"
    mdpgen.vdw_lambdas += "]"
    mdpgen.coul_lambdas += "]"

def array(mdpgen, lam, fol, mol, coul1="vdw", coul2="none", lbfgs=False,
        methylpyrrole=False):
    mdpgen.enermin_steep()
    do_set(mdpgen, lam, fol, mol, coul1, coul2)
    output("em_steep.mdp", mdpgen.compile())

    if lbfgs:
        mdpgen.enermin_lbfgs()
        do_set(mdpgen, lam, fol, mol, coul1, coul2)
        output("em_l-bfgs.mdp", mdpgen.compile())

    mdpgen.equilibrate_nvt()
    if methylpyrrole:
        mdpgen.modify_1meth()
    do_set(mdpgen, lam, fol, mol, coul1, coul2)
    output("nvt.mdp", mdpgen.compile())

    mdpgen.equilibrate_npt()
    if methylpyrrole:
        mdpgen.modify_1meth(pcouple=True)
    do_set(mdpgen, lam, fol, mol, coul1, coul2)
    output("npt.mdp", mdpgen.compile())

    mdpgen.production_md()
    if methylpyrrole:
        mdpgen.modify_1meth(pcouple=True)
    do_set(mdpgen, lam, fol, mol, coul1, coul2)
    output("md.mdp", mdpgen.compile())

def array2(mdpgen, lam, mol, fep, coul, vdw, coup1="vdw-q",
        coup2="none", lbfgs=False, methylpyrrole=False):

    mdpgen.enermin_steep()
    dos_set(mdpgen, lam, mol, coup1, coup2, fep, vdw, coul)
    output("em_steep.mdp", mdpgen.compile())

    if lbfgs:
        mdpgen.enermin_lbfgs()
        dos_set(mdpgen, lam, mol, coup1, coup2, fep, vdw, coul)
        output("em_l-bfgs.mdp", mdpgen.compile())

    mdpgen.equilibrate_nvt()
    if methylpyrrole:
        mdpgen.modify_1meth()
    dos_set(mdpgen, lam, mol, coup1, coup2, fep, vdw, coul)
    output("nvt.mdp", mdpgen.compile())

    mdpgen.equilibrate_npt()
    if methylpyrrole:
        mdpgen.modify_1meth(pcouple=True)
    dos_set(mdpgen, lam, mol, coup1, coup2, fep, vdw, coul)
    output("npt.mdp", mdpgen.compile())

    mdpgen.production_md()
    if methylpyrrole:
        mdpgen.modify_1meth(pcouple=True)
    dos_set(mdpgen, lam, mol, coup1, coup2, fep, vdw, coul)
    output("md.mdp", mdpgen.compile())

def launch(skip=False):
    jobname = "1methsolv"
    topfile = "1meth.top"
    grofile = "1meth.gro"
    mol = "TMP"

    mdpgen = MakeMDP()

    if not skip:
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
        array(mdpgen, lam, fol, mol, coul1, coul2, methylpyrrole=True)
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
        array(mdpgen, lam, fol, mol, methylpyrrole=True)
        output("job.slurm", outtext)
        os.system("sbatch job.slurm")
        os.chdir("..")

def launch_fep(skip=False, name="1meth", lname="1methylpyrrole"):
    jobname = name + "solv"
    topfile = name + ".top"
    grofile = name + ".gro"
    mol = "TMP"

    mdpgen = MakeMDP()

    if not skip:
        os.system("python -u ~/git/expdens/gromod.py -n t4" + name +
            "-in.gro -o ligand.gro -i TMP -v -s -t " + lname + " -m center")
        os.system("python -u ~/git/expdens/gromod.py -n t4" + name +
                "-in.gro -o solvent.gro -i SOL -v -t water")
        os.system("genbox -cp ligand.gro -cs solvent.gro" +
                " -ci solvent.gro -p " + topfile + " -o " + grofile +
                " -nmol 6500 -maxsol 10000")

    coup1 = "vdw-q"
    coup2 = "none"
    coul = [0.0, 0.3, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0]
    vdw = [0.0, 0.0, 0.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
        0.9, 1.0]
    n = len(coul)
    fep = [0.0] * n
    for l in range(n):
        lam = "{0:02}".format(l)
        folder = jobname + "-" + lam
        if not os.path.exists(folder):
            os.mkdir(folder)
        os.chdir(folder)
        slurm = MakeSLURM(jobname, "-" + lam, ".")
        slurm.double_precision = True
        outtext = slurm.compile(os.path.join("..", grofile),
            os.path.join("..", topfile), lam, use_lbfgs=False)
        array2(mdpgen, l, mol, fep, coul, vdw, coup1, coup2,
            methylpyrrole=True)
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
        parser.add_argument("-v", "--verbose", dest="verbose", action="count", help="set verbosity level [default: %(default)s]")
        parser.add_argument('-V', '--version', action='version', version=program_version_message)
        parser.add_argument('-s', '--skip', action='store_true',
            help="skip gromod")
        parser.add_argument('-f', '--fep', action='store_true',
            help="use fep explicitly")
        parser.add_argument('-n', '--name', help="name of job", default=None)
        parser.add_argument('-l', '--lname', help="title for grofile",
                            default=None)
        parser.add_argument(dest="paths", help="paths to folder(s) with source file(s) [default: %(default)s]", metavar="path", nargs='*')

        # Process arguments
        args = parser.parse_args()

        verbose = args.verbose

        if verbose > 0:
            print("Verbose mode on")

        if args.fep:
            if (args.name != None) and (args.lname != None):
                launch_fep(skip=args.skip, name=args.name, lname=args.lname)
            elif (args.name != None):
                launch_fep(skip=args.skip, name=args.name)
            elif (args.lname != None):
                launch_fep(skip=args.skip, lname=args.lname)
            else:
                launch_fep(skip=args.skip)
        else:
            launch(skip=args.skip)

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
