#!/usr/bin/env python
'''
Created on Jun 19, 2015

@author: Tyler P
'''
import sys, os
from optparse import OptionParser
import json
import datetime
import math
import backup
from parameters import Parameters
from parameters import Keys

verbose = 0
SUBS = dict()  # available subcommands, as a dictionary
def _subcomment():
    return "One of: " + ", ".join(SUBS.keys())

class MyKeys(Keys):
    """Container class for the options' key strings
    @see binding_ensemble_parameters()
    """

    def __init__(self):
        Keys.__init__(self)
        # Mostly values that can be specified by the user in the options file
        ################################################
        # General
        self.base_name = 'base-name'
        self.job_name = 'job-name'
        self._job_name = '_job-name'  # private
        self._params = '_params'  # private
        self._params_out = '_params-out'  # private
        self._calling_dir = '_calling-dir'  # private
        self.work_dir = 'sim-dir'
        self.script_dir = 'scripts-dir'
        self.ligand = 'ligand-res-name'

        section = "General"
        self.add_key(self.base_name, section, 'Name to prefix all files')
        self.add_keys(section, self.job_name, self.work_dir, self.script_dir,
                      self.ligand)
        ################################################
        # Sim Options
        self.partition = 'partition-name'
        self.sim_time = 'sim-time-ns'
        self.sim_weights = 'sim-initial-weights'
        self.sim_fixed_weights = 'sim-fixed-weights'
        self.sim_weight_values = 'sim-value-of-weights'
        self.sim_fep_values = 'sim-values-of-fep'
        self.sim_vdw_values = 'sim-values-of-vdw'
        self.sim_coul_values = 'sim-values-of-coul'
        self.sim_genxcoupled = 'sim-generate-x-coupled-states'
        self.sim_genxuncupld = 'sim-generate-x-uncoupled-states'
        self.sim_wgtxcoupled = 'sim-weight-x-coupled-states'
        self.sim_wgtxuncupld = 'sim-weight-x-uncoupled-states'
        self.sim_incrementor = 'sim-weight-incrementor'
        self.sim_init_lambda = 'sim-init-lambda'
        self.sim_fixed_lambda = 'sim-init-lambda-only_no-expdens'
        self.sim_use_gibbs = 'sim-use-gibbs-sampling'
        self.sim_gibbs_delta = 'sim-gibbs-max-step'
        self.sim_nstout = 'sim-nstout'
        self.sim_nst_mc = 'sim-nst-mc'
        self.sim_pressure = 'sim-pressure-coupling'

        section = "Simulation"
        self.add_key(self.partition, section, ('Known queues: economy,' +
                                               ' serial, parallel'))
        self.add_key(self.sim_genxcoupled, section, ("Number of states to" +
            " generate by adding entries (0.0's) at the beginning of the" +
            " state index list"))
        self.add_key(self.sim_genxuncupld, section, ("Number of states to" +
            " generate by adding entries (1.0's) at the end of the state" +
            " index list"))
        self.add_key(self.sim_wgtxcoupled, section, ("Emulate these number of" +
            " states by adding ln(x) to all but the beginning state"))
        self.add_key(self.sim_wgtxuncupld, section, ("Emulate these number of" +
            " states by adding ln(x) to just the end state"))
        self.add_keys(section, self.sim_time,
            self.sim_weights, self.sim_fixed_weights,
            self.sim_weight_values, self.sim_fep_values,
            self.sim_vdw_values, self.sim_coul_values,
            self.sim_incrementor, self.sim_init_lambda,
            self.sim_fixed_lambda, self.sim_use_gibbs,
            self.sim_gibbs_delta, self.sim_nstout,
            self.sim_nst_mc, self.sim_pressure)
        ################################################
        # MDP Array
        self.mdr_count = 'mdr-number-of-simulations'
        self.mdr_threads = 'mdr-number-of-threads'
        self.mdr_queue_time = 'mdr-queue-time'
        self.mdr_genseed = 'mdr-generate-seeds'

        section = "mdrun Array"
        self.add_key(self.mdr_threads, section, 'Max 20 threads (for now)')
        self.add_keys(self.mdr_count, self.mdr_queue_time, self.mdr_genseed)
        ################################################
        # Automation
        self.auto_time_equil = 'auto-equil-time-ns'
        self.auto_time_rand = 'auto-rand-time-ns'
        self.auto_time_array = 'auto-array-time-ns'
        self.auto_calc_wt = 'auto-calculate-walltime'

        section = "Automation"
        self.add_keys(section, self.auto_time_equil, self.auto_time_rand,
                      self.auto_time_array, self.auto_calc_wt)
        ################################################
        # Process Control
        self.verbosity = 'verbose-level'
        self.subcommand = 'subcommand'
        self.run_mdp = 'create-mdp'
        self.run_array = 'setup-array-of-jobs'
        self._submit = 'submit-jobs'  # private
        self._dryrun = 'dry-run'  # private
        self._chain_all = 'GEN_ALL'  # private

        section = "Process Control"
        self.add_key(self.subcommand, updater=_subcomment)
        self.add_keys(section, self.verbosity, self.run_mdp, self.run_array)

KEYS = MyKeys()

param = Parameters(KEYS)

def option_defaults():
    """Define the default parameters for param file.
    
    The inspiration for this procedure is taken from the GROMACS simulation
    package, and the file type .mdp

    @organization: Shirts Group
    @author: Tyler P
    
    @return: dict containing all possible options, and their default values
    """

    options = dict()

    # To help avoid spelling errors, the keys are kept in a container object,
    #    and each use of that key points back to the string listed there.
    ################################################
    # General
    options[KEYS.base_name] = None
    options[KEYS.job_name] = None
    options[KEYS.work_dir] = '.'
    options[KEYS.script_dir] = '.'
    options[KEYS.ligand] = 'TMP'
    ################################################
    # Process Control
    options[KEYS.verbosity] = 1  # Change output frequency and detail
    options[KEYS.subcommand] = 'exit'
    options[KEYS.run_mdp] = False
    options[KEYS.run_array] = False
    ################################################
    # Sim Array
    options[KEYS.partition] = 'serial'
    options[KEYS.sim_time] = 0.2  # ns
    options[KEYS.sim_weights] = False
    options[KEYS.sim_fixed_weights] = False
    options[KEYS.sim_weight_values] = [0.0, 1.57031, 2.43164, 3.62305, 3.87891,
        4.06836, 4.26562, 4.4043, 4.51367, 4.60156, 4.62695, 4.6582, 4.77344,
        4.93164]
    options[KEYS.sim_fep_values] = [0.0] * 14
    options[KEYS.sim_vdw_values] = [0.0, 0.0, 0.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5,
        0.6, 0.7, 0.8, 0.9, 1.0]
    options[KEYS.sim_coul_values] = [0.0, 0.3, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    options[KEYS.sim_genxcoupled] = 0
    options[KEYS.sim_genxuncupld] = 0
    options[KEYS.sim_wgtxcoupled] = 0
    options[KEYS.sim_wgtxuncupld] = 0
    options[KEYS.sim_incrementor] = 1
    options[KEYS.sim_init_lambda] = -1  # last index
    options[KEYS.sim_fixed_lambda] = False
    options[KEYS.sim_use_gibbs] = False
    options[KEYS.sim_gibbs_delta] = -1
    options[KEYS.sim_nstout] = 500
    options[KEYS.sim_nst_mc] = 50
    options[KEYS.sim_pressure] = True
    ################################################
    # MDP Array
    options[KEYS.mdr_count] = 1
    options[KEYS.mdr_threads] = 10
    options[KEYS.mdr_queue_time] = '1-00:00:00'
    options[KEYS.mdr_genseed] = False
    ################################################
    # Automation
    options[KEYS.auto_time_equil] = 0.2  # ns
    options[KEYS.auto_time_rand] = 0.2  # ns
    options[KEYS.auto_time_array] = 0.2  # ns
    options[KEYS.auto_calc_wt] = True

    return options

# provide custom default library
param.option_defaults = option_defaults

#################################################################################
#################################################################################

def place_simarray_vars(fields):
    fields['outputfile'] = os.path.join(fields['path'], fields['job-name'] +
        fields['suffix'], fields['job-name'])
    fields['workdir-name'] = fields['job-name'] + fields['suffix']
    fields['workdir'] = os.path.join(fields['path'], fields['workdir-name'])
    fields['infolder'] = os.path.join('..', '')
    fields['jobstatus'] = os.path.join('..', 'jobstatus.txt')

    fields['move-cmds'] = ''
    if 'move-par' in fields and fields['move-par']:
        fields['move-cmds'] += "mv {infolder}{param-out} {param-out} \n"
    if 'move-slurm' in fields and fields['move-slurm']:
        fields['move-cmds'] += "mv {infolder}{slurm-out} {slurm-out} \n"
    if 'move-gro' in fields and fields['move-gro']:
        fields['move-cmds'] += "mv {infolder}{gro-in}-in.gro {gro-in}-in.gro \n"

    extins = 'extra-instructions'
    if not extins in fields:
        fields[extins] = ''
    else:
        fields[extins] = fields[extins].format(**fields)

    if fields['calc-walltime']:
        rate = 7.5 / 20  # 7.5 ns/day per 20 cores
        speed = rate * fields['ntasks']  # ns/day
        est = fields['time-ns'] / speed  # days
        time_ = est * 1.20  # add 20% buffer
        if time_ > 7:
            time_ = 7  # seven day limit on serial queue
        frac, days = math.modf(time_)
        frac, hours = math.modf(frac * 24.0)
        frac, mins = math.modf(frac * 60.0)
        frac, secs = math.modf(frac * 60.0)

        days = int(days)
        hours = int(hours)
        mins = int(mins)
        secs = int(secs)
        fields['queue-time'] = '{0}-{1:02}:{2:02}:{3:02}'.format(days, hours,
            mins, secs)

    text_ = """#!/bin/sh
#SBATCH --job-name={job-name}{suffix}
#SBATCH --partition={partition}
#SBATCH --nodes=1
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

# submit a single SLURM job for a free energy calc
module load gromacs/5.0.2-sse
export THREADINFO="-nt {ntasks} "
export GMX_SUPPRESS_DUMP=1 #prevent step file output when parameters don't converge

#export OMP_NUM_THREADS={ntasks}

sleep 1

#Change the job status to 'SUBMITTED'
echo "SUBMITTED {job-name}{suffix} `date`" >> {jobstatus}
echo "Job Started at `date`"

#NVT PRODUCTION
grompp_d -c {infolder}{gro-in}-in.gro -p {infolder}{base}.top -n {infolder}{base}.ndx -f {mdp-in}.mdp -o {job-name}.tpr -maxwarn 15
if [ -f {job-name}.tpr ];
then
    echo "MD Simulation launching..."
    mdrun_d ${{THREADINFO}} -deffnm {job-name}
else
    echo "Error generating {job-name}.tpr"
fi
{move-cmds}
{extra-instructions}
echo "FINISHED {job-name}{suffix} `date`" >> {jobstatus}
# print end time
echo
echo "Job Ended at `date`"
echo "###################################################################"
""".format(**fields)

    return text_

def place_simarray_mdp_vars(fields):
    fields['dt'] = 0.002;  # ps
    fields['nstrun'] = fields['time-ns'] * (1000 / fields['dt'])

    fields['fep-lambdas'] = ''
    fields['coul-lambdas'] = ''
    fields['vdw-lambdas'] = ''
    for f, c, v in zip(fields['fep'], fields['coul'], fields['vdw']):
        fields['fep-lambdas'] += '{0:0.1f} '.format(f)
        fields['coul-lambdas'] += '{0:0.1f} '.format(c)
        fields['vdw-lambdas'] += '{0:0.1f} '.format(v)

    fields['wl-weights'] = ''

    if fields['initial-weights']:
        for w in fields['weights']:
            fields['wl-weights'] += '{0:0.5f} '.format(w)
        if not len(fields['fep']) == len(fields['weights']):
            print('Number of weights not equal to number of states, please' +
                ' address')
            return False
        fields['ilw'] = ''
    else:
        fields['ilw'] = '; '

    if fields['fixed-weights'] and fields['initial-weights']:
        fields['weights-equil'] = 'yes'
        fields['wed'] = '; '
    else:
        fields['weights-equil'] = 'wl-delta'
        fields['wed'] = ''

    fields['init-lambda'] = (fields['init-state-index'] + 1) / len(fields['fep'])
    if fields['init-lambda'] < 0:
        fields['init-lambda'] = 0
    elif fields['init-lambda'] > 1:
        fields['init-lambda'] = 1

    if fields['fixed-lambda']:
        fields['free-energy-flag'] = 'yes'
        fields['ilv'] = ''
        fields['ils'] = '; '
    else:
        fields['free-energy-flag'] = 'expanded'
        fields['ilv'] = '; '
        fields['ils'] = ''

    if fields['use-pressure']:
        fields['npt'] = ''
        fields['pcoupl'] = 'MTTK'
    else:
        fields['npt'] = '; '
        fields['pcoupl'] = 'no'


    text_ = """; RUN CONTROL PARAMETERS = 
integrator               = md-vv
; start time and timestep in ps = 
tinit                    = 0
dt                       = {dt:0.3f}
nsteps                   = {nstrun:<12.0f}; {time-ns:0.2f} ns
; mode for center of mass motion removal = 
comm-mode                = Linear
; number of steps for center of mass motion removal = 
nstcomm                  = 1
; group(s) for center of mass motion removal = 
;comm-grps                = 
nsttcouple               = 1
nstpcouple               = 1
; NEIGHBORSEARCHING PARAMETERS = 
; nblist update frequency = 
nstlist                  = 10
; ns algorithm (simple or grid) = 
ns_type                  = grid
; Periodic boundary conditions: xyz or no = 
pbc                      = xyz
; nblist cut-off         = 
rlist                    = 1.2

; OUTPUT CONTROL OPTIONS
; Output frequency for coords (x), velocities (v) and forces (f)
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
; Output frequency for energies to log file and energy file
nstlog                   = {nstout:0.0f}  ; changing this allows you to see the frequency the weights are computed.
nstcalcenergy            = 1
nstenergy                = {nstout:0.0f}
; Output frequency and precision for .xtc file
nstxout-compressed       = {nstout:0.0f}  ; change this to control how frequently the structures are printed out.
compressed-x-precision   = 1000

; OPTIONS FOR ELECTROSTATICS AND VDW = 
; Method for doing electrostatics = 
cutoff-scheme            = group
coulombtype              = PME
coulomb-modifier         = Potential-Switch
rcoulomb-switch          = 0.88
rcoulomb                 = 0.9

; Method for doing Van der Waals = 
vdw-type                 = Cut-off
vdw-modifier             = Potential-switch
; cut-off lengths        = 
rvdw-switch              = 0.85
rvdw                     = 0.9
; Apply long range dispersion corrections for Energy and Pressure = 
DispCorr                 = AllEnerPres
; Spacing for the PME/PPPM FFT grid = 
fourierspacing           = 0.12
; FFT grid size, when a value is 0 fourierspacing will be used = 
fourier_nx               = 0
fourier_ny               = 0
fourier_nz               = 0
; EWALD/PME/PPPM parameters = 
pme_order                = 4
ewald_rtol               = 1e-05
ewald_geometry           = 3d

; OPTIONS FOR BONDS     = 
constraints              = hbonds ; constrain bonds to hydrogen
; Type of constraint algorithm = 
constraint-algorithm     = shake
; Highest order in the expansion of the constraint coupling matrix = 
shake-tol                = 1e-12 ; 5e-6 ~TP

; GENERATE VELOCITIES FOR STARTUP RUN = 
gen_vel                  = yes
gen_temp                 = 300.0 ; K
gen_seed                 = {gen-seed}


; OPTIONS FOR WEAK COUPLING ALGORITHMS = 

; Groups to couple separately = 
tc-grps                  = System
; Time constant (ps) and reference temperature (K) = 
tcoupl                   = nose-hoover
tau_t                    = 0.5; 1.0 ~TP
ref_t                    = 300.0 

; Pressure coupling      = 
Pcoupl                   = {pcoupl}
{npt}Pcoupltype               = isotropic
{npt}; Time constant (ps), compressibility (1/bar) and reference P (bar) = 
{npt}tau_p                    = 20.0 ; 5.0 ; ps ~TP
{npt}compressibility          = 4.5e-5 ; 1/bar
{npt}ref_p                    = 1.0 ; bar

; OPTIONS FOR EXPANDED ENSEMBLE SIMULATIONS
; Free energy control stuff  
free-energy               = {free-energy-flag}
sc-alpha                  = 0.5
couple-moltype            = {ligand}
couple-lambda0            = vdw-q
couple-lambda1            = none
couple-intramol           = no
fep-lambdas               = {fep-lambdas}
coul-lambdas              = {coul-lambdas}
vdw-lambdas               = {vdw-lambdas}
{ilv}init-lambda               = {init-lambda}
{ils}init-lambda-state         = {init-state-index}
symmetrized-transition-matrix = yes
nst-transition-matrix     = 100000
nstdhdl                   = {nstout:0.0f}
dhdl-print-energy         = total

; expanded ensemble
nstexpanded              = {nst-mc:0.0f}
lmc-stats                = wang-landau
lmc-move                 = {metropolis}
lmc-weights-equil        = {weights-equil}
{wed}weight-equil-wl-delta    = 0.001
lmc-seed                 = {lmc-seed}
{ilw}init-lambda-weights      = {wl-weights}
lmc-gibbsdelta           = {gibbs-delta}

; Seed for Monte Carlo in lambda space
wl-scale                 = 0.5
wl-ratio                 = 0.8
init-wl-delta            = {incrementor:0.3f}
wl-oneovert              = yes
""".format(**fields)

    return text_

def generate(opts):
    randseed = opts[KEYS.mdr_genseed]
    submit = opts[KEYS._submit]
    job_name = opts[KEYS.job_name]
    if not job_name:
        print('Job name not specified')
        return
    base_name = opts[KEYS.base_name]
    if not base_name:
        base_name = job_name

    cur_dir = os.getcwd()
    dir_ = backup.expandpath(opts[KEYS.script_dir])
    if not os.path.exists(dir_):
        os.mkdir(dir_)
    os.chdir(dir_)
    path = backup.expandpath(opts[KEYS.work_dir])
    for i in range(opts[KEYS.mdr_count]):
        suffix = '_{0:0>2}'.format(i)
        if opts[KEYS.mdr_count] < 2:
            suffix = ''
        file_name = job_name + suffix + '.slurm'

        fields = dict()
        fields['base'] = base_name
        fields['job-name'] = job_name
        fields['suffix'] = suffix
        fields['path'] = path
        fields['partition'] = opts[KEYS.partition]
        fields['ntasks'] = opts[KEYS.mdr_threads]
        fields['queue-time'] = opts[KEYS.mdr_queue_time]
        fields['calc-walltime'] = opts[KEYS.auto_calc_wt]
        fields['gro-in'] = base_name + suffix
        fields['mdp-in'] = job_name + suffix
        fields['pathsep'] = os.path.join('A', '').replace('A', '')
        fields['move-par'] = True
        fields['move-slurm'] = True
        fields['slurm-out'] = file_name
        if opts[KEYS._chain_all] and (opts[KEYS.subcommand] == 'equil' or
            opts[KEYS.subcommand] == 'rand' or
            opts[KEYS.subcommand] == 'array'):

            fields['extra-instructions'] = """
# genscript {this-cmd}"""
            if opts[KEYS.subcommand] == 'array':
                fields['extra-instructions'] += """
cd {workdir}
"""
            else:
                fields['extra-instructions'] += """
cd {calling-dir}
python {path-to-here} {next-cmd} {special-flag} --submit {params-option}
cd {workdir}
"""
            fields['this-cmd'] = opts[KEYS.subcommand]
            if opts[KEYS.subcommand] == 'equil':
                fields['next-cmd'] = 'rand'
            if opts[KEYS.subcommand] == 'rand':
                fields['next-cmd'] = 'array'
            fields['script-dir'] = dir_
            fields['path-to-here'] = os.path.realpath(__file__)
            fields['special-flag'] = KEYS._chain_all
            fields['calling-dir'] = opts[KEYS._calling_dir]
            fields['param-out'] = opts[KEYS._params_out]
            fields['params-option'] = ""
            if opts[KEYS._params]:
                fields['params-option'] = "--par " + opts[KEYS._params]

        cont = False
        if randseed:
            fields['gro-in'] = base_name
            try:
                import random
                # random.seed(0)
                seed = random.randint(0, 65536)
                fields['gen-seed'] = seed + i
                fields['lmc-seed'] = seed + i
                cont = parameters(opts, name=fields['mdp-in'], fields=fields)
            except IOError as ioex:
                print(ioex)
        else:
            try:
                cont = parameters(opts, name=job_name + suffix, fields=fields)
            except IOError as ioex:
                print(ioex)

        if not cont:
            return

        file_ = open(file_name, 'w')
        file_.write(place_simarray_vars(fields))
        file_.close()
        if submit:
            dir_name = os.path.join(path, job_name + suffix)
            if not os.path.exists(dir_name):
                os.mkdir(dir_name)
            # yes, this seems to be unnecessary, but it helps in debugging if
            #    directory doesn't actually have to exist to run without submit
            os.system("mv " + fields['mdp-in'] + ".mdp " +
                os.path.join(dir_name, fields['mdp-in']) + ".mdp")
            if not opts[KEYS._dryrun]:
                os.system("sbatch " + file_name)
                print("Sbatch'd Job " + job_name + suffix)
            else:
                print("DRYRUN: Would sbatch job " + job_name + suffix)

    os.chdir(cur_dir)

def parameters(opts, dir_='.', name=None, fields=None):
    job_name = opts[KEYS.job_name]
    if not job_name:
        print('Job name not specified')
        return

    if not name:
        name = job_name

    cur_dir = os.getcwd()
    dir_ = os.path.expandvars(os.path.expanduser(dir_))
    if not os.path.exists(dir_):
        os.mkdir(dir_)
    os.chdir(dir_)

    if not fields:
        fields = dict()
    fields['ligand'] = opts[KEYS.ligand]
    fields['time-ns'] = opts[KEYS.sim_time]
    fields['initial-weights'] = opts[KEYS.sim_weights]
    fields['fixed-weights'] = opts[KEYS.sim_fixed_weights]
    fields['fixed-lambda'] = opts[KEYS.sim_fixed_lambda]

    if not 'gen-seed' in fields:
        fields['gen-seed'] = 10200
    if not 'lmc-seed' in fields:
        fields['lmc-seed'] = 10200

    fep = opts[KEYS.sim_fep_values]
    coul = opts[KEYS.sim_coul_values]
    vdw = opts[KEYS.sim_vdw_values]
    weights = opts[KEYS.sim_weight_values]
    if verbose > 2:
        print('Debug, genstates')

    if verbose > 2:
        print('Value coup {0}'.format(opts[KEYS.sim_wgtxcoupled]))
    if opts[KEYS.sim_wgtxcoupled] > 0:
        if verbose > 2:
            print('making more states, coupled')
        for i in range(1, len(weights)):
            weights[i] = weights[i] - math.log(opts[KEYS.sim_wgtxcoupled])
    if verbose > 2:
        print('Value uncoup {0}'.format(opts[KEYS.sim_wgtxuncupld]))
    if opts[KEYS.sim_wgtxuncupld] > 0:
        if verbose > 2:
            print('making more states, uncoupled')
        if len(weights) > 0:
            weights[-1] = weights[-1] + math.log(opts[KEYS.sim_wgtxuncupld])
    if opts[KEYS.sim_genxcoupled] > 0:
        fep = [0.0] * int(opts[KEYS.sim_genxcoupled]) + fep
        vdw = [0.0] * int(opts[KEYS.sim_genxcoupled]) + vdw
        coul = [0.0] * int(opts[KEYS.sim_genxcoupled]) + coul
        weights = [weights[0]] * int(opts[KEYS.sim_genxcoupled]) + weights
    if opts[KEYS.sim_genxuncupld] > 0:
        fep = fep + [0.0] * int(opts[KEYS.sim_genxuncupld])
        vdw = vdw + [1.0] * int(opts[KEYS.sim_genxuncupld])
        coul = coul + [1.0] * int(opts[KEYS.sim_genxuncupld])
        weights = weights + [weights[-1]] * int(opts[KEYS.sim_genxuncupld])
    if verbose > 2:
        print(weights)

    fields['fep'] = fep
    fields['coul'] = coul
    fields['vdw'] = vdw
    if opts[KEYS.sim_init_lambda] >= 0:
        fields['init-state-index'] = opts[KEYS.sim_init_lambda]
    else:
        fields['init-state-index'] = len(fields['fep']) + opts[
            KEYS.sim_init_lambda]
    fields['weights'] = weights
    fields['incrementor'] = opts[KEYS.sim_incrementor]

    fields['metropolis'] = 'metropolis'
    if opts[KEYS.sim_use_gibbs]:
        fields['metropolis'] = 'metropolized-gibbs'
    fields['gibbs-delta'] = int(opts[KEYS.sim_gibbs_delta])

    fields['nstout'] = opts[KEYS.sim_nstout]
    fields['nst-mc'] = opts[KEYS.sim_nst_mc]

    fields['use-pressure'] = opts[KEYS.sim_pressure]

    file_name = name + '.mdp'
    file_ = open(file_name, 'w')
    value = place_simarray_mdp_vars(fields)
    if value:
        file_.write(value)
    file_.close()

    os.chdir(cur_dir)
    return value

def gen_exit(_):
    sys.exit(0)

def gen_all(opts):
    opts[KEYS._chain_all] = True

    opts2 = dict(opts)  # backup
    gen_equil(opts)

    if opts[KEYS._dryrun]:
        opts = dict(opts2)  # restore backup
        gen_rand(opts)

        opts = dict(opts2)  # restore backup
        gen_array(opts)

def gen_equil(opts):
    opts[KEYS.subcommand] = 'equil'
    opts[KEYS.mdr_count] = 1
    opts[KEYS.sim_fixed_weights] = False
    opts[KEYS.sim_weights] = False
    opts[KEYS.sim_fixed_lambda] = False
    opts[KEYS.sim_time] = opts[KEYS.auto_time_equil]
    if opts[KEYS._job_name]:
        opts[KEYS.job_name] = opts[KEYS._job_name] + '-equil'
    generate(opts)

def gen_rand(opts):
    if opts[KEYS._chain_all]:
        # attempt to parse output of the equil step
        try:
            cur_dir = os.getcwd()
            work_dir = backup.expandpath(opts[KEYS.work_dir])
            os.chdir(os.path.join(work_dir, opts[KEYS.job_name] + '-equil'))
            file_name = '{0}-equil.log'.format(opts[KEYS.job_name])
            if not os.path.exists(file_name):
                raise Exception('Simulation file not found ' + file_name)
            os.system('tail -300 {0} &>temp.log'.format(file_name))
            line_count = -1
            weights = []
            for line in open('temp.log'):
                if 'N' in line and 'Count' in line and 'G(in kT)' in line:
                    weights = []
                    line_count = 1
                if line_count > 0:
                    line_count -= 1
                elif line_count == 0:
                    if line.isspace():
                        line_count = -1
                    else:
                        splt = line.split()
                        weight = float(splt[5])
                        weights.append(weight)
            opts[KEYS.sim_weight_values] = weights
            os.remove('temp.log')

            os.chdir(work_dir)
            param.write_options('rand.save', {'weights': weights},
                None, None)

            os.chdir(cur_dir)
        except Exception as ex:
            print(ex)
            os.chdir(cur_dir)
            if not opts[KEYS._dryrun]:
                return

    opts[KEYS.subcommand] = 'rand'
    opts[KEYS.mdr_count] = 1
    opts[KEYS.sim_fixed_weights] = False
    opts[KEYS.sim_weights] = False
    opts[KEYS.sim_fixed_lambda] = True
    opts[KEYS.sim_time] = opts[KEYS.auto_time_rand]
    if opts[KEYS._job_name]:
        opts[KEYS.job_name] = opts[KEYS._job_name] + '-rand'
    generate(opts)

def gen_array(opts):
    if opts[KEYS._chain_all]:
        # attempt to parse output of the rand step
        try:
            cur_dir = os.getcwd()
            work_dir = backup.expandpath(opts[KEYS.work_dir])
            os.chdir(work_dir)
            file_name = 'rand.save'
            if not os.path.exists(file_name):
                raise Exception('Save file not found ' + file_name)
            save_data = param.parse_options(file_name, dict())
            opts[KEYS.sim_weight_values] = save_data['weights']

            os.chdir(os.path.join(work_dir, opts[KEYS.job_name] + '-rand'))
            file_name = '{0}-rand.log'.format(opts[KEYS.job_name])
            if not os.path.exists(file_name):
                raise Exception('Simulation file not found ' + file_name)
            os.system('tail -300 {0} &>temp.log'.format(file_name))
            line_count = -1
            num_of_steps = 0
            for line in open('temp.log'):
                if 'Step' in line and 'Time' in line and 'Lambda' in line:
                    line_count = 1
                if line_count == 0:
                    splt = line.split()
                    num_of_steps = int(splt[0])
                if line_count >= 0:
                    line_count -= 1
            num_frames = math.floor(num_of_steps / 500.0)  # hardcoded 500
            # 500 is taken from nstxout-compressed variable in mdp
            os.remove('temp.log')

            base_name = opts[KEYS.base_name]
            if not base_name:
                base_name = opts[KEYS.job_name]
            xtc_name = '{0}-rand.xtc'.format(opts[KEYS.job_name])
            tpr_name = '{0}-rand.tpr'.format(opts[KEYS.job_name])
            if not os.path.exists(xtc_name):
                raise Exception('Simulation files not found ' + xtc_name +
                    ' and/or ' + tpr_name)
            cmnd = ('echo "System" | trjconv_d -f {xtc} -s {tpr} -o {gro}' +
                ' -b {frame} -e {frame} &>../trjconv{index}.log')
            segments = (opts[KEYS.mdr_count] - 1)
            if segments == 0:
                spacing = 1  # will only be multiplied by 0
            else:
                spacing = num_frames / segments

            if ('GMXBIN' in os.environ and os.environ['GMXBIN'] and
                os.environ['GMXBIN'] not in os.environ['PATH']):
                print('Adding GMXBIN to PATH...')
                os.environ['PATH'] = os.environ['PATH'] + ':' + os.environ['GMXBIN']

            for i in range(opts[KEYS.mdr_count]):
                suffix = '_{0:0>2}'.format(i)
                frame_num = math.floor(i * spacing)  # evenly spaced frames
                gro_name = '../' + base_name + suffix + '-in.gro'
                fmt = cmnd.format(xtc=xtc_name, tpr=tpr_name, gro=gro_name,
                    frame=frame_num, index=suffix)
                os.system(fmt)
                if not os.path.exists(gro_name):
                    raise Exception('Extraction of starting frames' +
                        ' unsuccessful - Err: ' + gro_name)

            os.chdir(cur_dir)
        except Exception as ex:
            print(ex)
            os.chdir(cur_dir)
            if not opts[KEYS._dryrun]:
                return

    opts[KEYS.subcommand] = 'array'
    opts[KEYS.sim_fixed_weights] = True
    opts[KEYS.sim_weights] = True
    opts[KEYS.sim_fixed_lambda] = False
    opts[KEYS.sim_time] = opts[KEYS.auto_time_array]
    if opts[KEYS._job_name]:
        opts[KEYS.job_name] = opts[KEYS._job_name]
    generate(opts)

def gen_opt(opts):
    if opts[KEYS.run_mdp]:
        parameters(opts)
    if opts[KEYS.run_array]:
        generate(opts)

def main(argv=None):
    if argv == None:
        argv = sys.argv[1:]

    cur_dir = os.getcwd()
    name_ = os.path.basename(__file__).replace('.py', '')
    print('Invocation of %s:\n\t%s ' % (__file__, name_) + " ".join(argv) + '\n')

    opts = option_defaults()

    parser = OptionParser()
    parser.set_description("Generates Rivanna SLURM scripts for" +
        " running Gromacs simulations")
    parser.add_option("-v", "--verbose", help="Increase output frequency" +
        " and detail. Stacks three times.", action="count")
    parser.add_option("--par", help="Optional configuration file to specify" +
        " command line parameters and more.",
        default=None, metavar="file.par")
    parser.add_option("--slurm", help="""Prepare scripts to run on Rivanna""",
        default=None, action='store_true')
    parser.add_option("--mdp", help="Output the Molecular Dynamics Parameters" +
        " (.mdp) configuration file",
        default=None, action='store_true')
    parser.add_option("-n", help="""Job name for output""", metavar="sim")
    parser.add_option("-b", help="""Base name for job input""", metavar="sim")
    parser.add_option("--submit", action='store_true')
    parser.add_option("--dryrun", action='store_true')

    parser.set_usage(name_ + '.py' + " subcommand [options]\n"
        "\tThe following subcommands may be used:\n" +
        "\t  {0:<12s} - run an entire array production\n".format('all') +
        "\t  {0:<12s} - run the Wang-Landau equilibriation\n".format('equil') +
        "\t  {0:<12s} - run a position randomizing simulation\n".format('rand') +
        "\t  {0:<12s} - run a batch of simulations\n".format('array') +
        "\t  {0:<12s} - output a formatted mdp file\n".format('mdp') +
        "\t  {0:<12s} - perform tasks based on a param file\n".format('opt') +
        "\t  {0:<12s} - print usage and exit".format('exit'))

    if len(argv) < 1:
        parser.print_help()
        sys.exit(0)

    (options, args) = parser.parse_args(argv)

    # handle options, reading from file if requested
    if options.par:
        print('Reading parameters from ' + options.par)
        opts = param.parse_options(options.par, opts)
        if not opts:
            print('Error reading from parameters file.')
            sys.exit(1)

    global verbose  # ensure that we are talking about the same verbose here
    if options.verbose:
        verbose = options.verbose
    else:
        verbose = opts[KEYS.verbosity]
    opts[KEYS.verbosity] = verbose

    if args:
        subcommand = args[0]
        if subcommand not in SUBS:
            subcommand = 'exit'
    else:
        subcommand = opts[KEYS.subcommand]

    if subcommand == 'exit':
        parser.print_help()

    if options.n:
        opts[KEYS.job_name] = options.n
    if options.b:
        opts[KEYS.base_name] = options.b
    if options.mdp:
        opts[KEYS.run_mdp] = options.mdp
    if options.slurm:
        opts[KEYS.run_array] = options.slurm
    opts[KEYS.subcommand] = subcommand

    # output run options
    if opts[KEYS.job_name]:
        param_name = opts[KEYS.job_name] + param.file_ext
        param_out = backup.backup_file('./', param_name, verbose=verbose)
    else:
        param_name = param.file_ext
        param_out = backup.backup_file('', param_name, verbose=verbose)
    param.write_options(param_out, opts)

    opts[KEYS._calling_dir] = cur_dir
    opts[KEYS._params] = options.par
    opts[KEYS._params_out] = param_name
    opts[KEYS._chain_all] = False
    if args and len(args) > 1:
        flag = args[1]
        if flag == KEYS._chain_all:
            opts[KEYS._chain_all] = True

    opts[KEYS._job_name] = opts[KEYS.job_name]  # backup the original name
    opts[KEYS._submit] = options.submit or options.dryrun
    opts[KEYS._dryrun] = options.dryrun

    # perform requested file output and job submissions
    SUBS[subcommand](opts)

SUBS.update({'exit': gen_exit, 'all': gen_all, 'equil': gen_equil,
    'rand': gen_rand, 'array': gen_array, 'mdp': parameters,
    'opt': gen_opt})
KEYS.update()
if __name__ == '__main__':
    # initialize the subcommand list
    main(sys.argv[1:])















