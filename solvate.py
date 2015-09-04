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
from argparse import RawDescriptionHelpFormatter
from openpyxl.reader.strings import get_text

__all__ = []
__version__ = 0.1
__date__ = '2015-09-04'
__updated__ = '2015-09-04'

DEBUG = 0
TESTRUN = 1
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

def em_steep_mdp(lam="0.0", fol="0.05"):
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
""".format(lam=lam, fol=fol)

def em_lbfgs_mdp(lam="0.0", fol="0.05"):
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
""".format(lam=lam, fol=fol)

def nvt_mdp(lam="0.0", fol="0.05"):
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
""".format(lam=lam, fol=fol)

def npt_mdp(lam="0.0", fol="0.05"):
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
""".format(lam=lam, fol=fol)

def md_mdp(lam="0.0", fol="0.05"):
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
""".format(lam=lam, fol=fol)

class MakeMDP:

    def __init__(self):
        self.integrator = "steep"
        self.nsteps = "5000"
        self.emtol = "100"
        self.emstep = "0.01"
        self.niter = "20"
        self.nbfgscorr = "10"
        self.nstlog = "1"
        self.nstenergy = "1"
        self.nstlist = "1"
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
        self.tcoupl = "no"
        self.pcoupl = "no"
        self.free_energy = "yes"
        self.init_lambda = "0.0"
        self.delta_lambda = "0"
        self.foreign_lambda = "0.05"
        self.sc_alpha = "0.5"
        self.sc_power = "1.0"
        self.sc_sigma = "0.3"
        self.couple_moltype = "Methane  ; name of moleculetype to decouple"
        self.couple_lambda0 = "vdw      ; only van der Waals interactions"
        self.couple_lambda1 = "none     ; turn off everything, in this case only vdW"
        self.couple_intramol = "no"
        self.nstdhdl = "10"
        self.gen_vel = "no "
        self.constraints = "h-bonds  ; we only have C-H bonds here"
        self.constraint_algorithm = "lincs"
        self.continuation = "no"
        self.lincs_order = "12"

    def enermin_steep(self):
        self.integrator = "steep"

    def enermin_lbfgs(self):
        self.integrator = "l-bfgs"
        self.define = "-DFLEXIBLE"
        self.constraints = "none     ; L-BFGS doesn't work with constraints "

    def nvt(self):
        self.nsteps = "50000    ; 100 ps"

    def npt(self):
        pass

    def md(self):
        pass

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
        text += self.free_energy()
        text += self.velocities()
        text += self.bond_constraints()
        return text

    def run_control(self):
        text = """; Run control
integrator               = steep 
nsteps                   = 5000
"""
        return text

    def em_criteria(self):
        text = """; EM criteria and other stuff
emtol                    = 100
emstep                   = 0.01
niter                    = 20
nbfgscorr                = 10
"""
        return text

    def output_control(self):
        text = """; Output control
nstlog                   = 1
nstenergy                = 1
"""
        return text

    def neighborsearching(self):
        text = """; Neighborsearching and short-range nonbonded interactions
nstlist                  = 1
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.0
"""
        return text

    def electrostatics(self):
        text = """; Electrostatics
coulombtype              = PME
rcoulomb                 = 1.0
"""
        return text

    def vanderWaals(self):
        text = """; van der Waals
vdw-type                 = switch
rvdw-switch              = 0.8
rvdw                     = 0.9
"""
        return text

    def corrections(self):
        text = """; Apply long range dispersion corrections for Energy and Pressure
DispCorr                  = EnerPres
"""
        return text

    def pme_pppm(self):
        text = """; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; EWALD/PME/PPPM parameters
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
optimize_fft             = no
"""
        return text

    def coupling(self):
        text = """; Temperature and pressure coupling are off during EM
tcoupl                   = no
pcoupl                   = no
"""
        return text

    def free_energy(self):
        text = """; Free energy control stuff
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
"""
        return text

    def velocities(self):
        text = """; Generate velocities to start
gen_vel                  = no 
"""
        return text

    def bond_constraints(self):
        text = """; options for bonds
"""
        if self.constraints:
            text += """constraints              = h-bonds  ; we only have C-H bonds here
""".format(self.constraints)
        text += """; Type of constraint algorithm
"""
        if self.contstraint_algorithm:
            text += """constraint-algorithm     = {0}
""".format(self.constraint_algorithm)
        text += """; Do not constrain the starting configuration
"""
        if self.continuation:
            text += """continuation             = {0}
""".format(self.continuation)
        text += """; Highest order in the expansion of the constraint coupling matrix
"""
        if self.lincs_order:
            text += """lincs-order              = {0}
""".format(self.lincs_order)
        return text

class MakeSLURM:

    def __init__(self, job_name, suffix, workdir, ntasks=2):
        self.input_gro = "meth.gro"
        self.double_precision = True
        self.use_mpi = False
        self.gromacs5 = False
        self.calc_wt = True
        self.fields_general = dict()
        self.fields_header = dict()
        self.fields_general['job-name'] = job_name
        self.fields_general['suffix'] = suffix
        self.fields_general['ntasks'] = ntasks
        self.fields_general['workdir'] = workdir
        self.fields_header['outputfile'] = os.path.join(workdir, job_name)

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
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtp4kc@virginia.edu
#SBATCH --comment="#SBATCH --array-0,1,3"
""".format(**fields)
        return self.header

    def get_text(self):
        return self.get_header() + """
module load jtp4kc
module load gromacs-jtp4kc
        
#################################
# ENERGY MINIMIZATION 1: STEEP  #
#################################
echo "Starting minimization for lambda = $LAMBDA..." 

# Iterative calls to grompp and mdrun to run the simulations

grompp{_d} -f em_steep.mdp -c {gro} -p {top} -o mins.tpr

mdrun{_d} -nt 2 -deffnm min


#################################
# ENERGY MINIMIZATION 2: L-BFGS #
#################################
grompp{_d} -f em_l-bfgs.mdp -c mins.gro -p {top} -o minl.tpr
# Run L-BFGS in serial (cannot be run in parallel)
mdrun{_d} -nt 1 -deffnm minl

echo "Minimization complete."


#####################
# NVT EQUILIBRATION #
#####################
echo "Starting constant volume equilibration..."

grompp{_d} -f nvt.mdp -c minl.gro -p {top} -o nvt.tpr
mdrun{_d} -nt 2 -deffnm nvt

echo "Constant volume equilibration complete."


#####################
# NPT EQUILIBRATION #
#####################
echo "Starting constant pressure equilibration..."

grompp{_d} -f npt.mdp -c nvt.gro -p {top} -t nvt.cpt -o npt.tpr
mdrun{_d} -nt 2 -deffnm npt

echo "Constant pressure equilibration complete."


#################
# PRODUCTION MD #
#################
echo "Starting production MD simulation..."

grompp{_d} -f md.mdp -c npt.gro -p {top} -t npt.cpt -o md.tpr
mdrun{_d} -nt 2 -deffnm md

echo "Production MD complete."

# End
echo "Ending. Job completed for lambda = {lam}"
"""

    def compile(self, gro, top, lam):
        _d = ""
        if self.double_precision:
            _d = "_d"
        return self.get_text().format(_d=_d, gro=gro, top=top, lam=lam)

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

def launch():
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
        slurm = MakeSLURM(jobname, "_" + lam, "")
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
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
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
