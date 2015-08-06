#!/usr/bin/env python
'''
Created on Jun 19, 2015

@author: Tyler P
'''
import sys, os
from optparse import OptionParser
import math
import backup
from param_versions.version_2_0 import Parameters
from param_versions.version_2_0 import Keys

verbose = 0
SUBS = dict()  # available subcommands, as a dictionary
def _subcomment():
    return "One of: " + ", ".join(SUBS.keys())

class MyKeys(Keys):
    """Container class for the options' key strings
    @see option_defaults()
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
        self.sim_use_mpi = "sim-use-mpi-parallelization"
        self.sim_time = 'sim-time-ns'
        self.sim_temperature = "sim-temperature"
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
        self.sim_use_gibbs = 'sim-use-gibbs-state-sampling'
        self.sim_use_metro = 'sim-use-metropolis-rejection'
        self.sim_gibbs_delta = 'sim-gibbs-max-step'
        self.sim_nstout = 'sim-nstout'
        self.sim_nst_mc = 'sim-nst-mc'
        self.sim_temp_alg = 'sim-temperature-algorithm'
        self.sim_pressure = 'sim-pressure-coupling'
        self.sim_precision = 'sim-use-double-precision'
        self.sim_gromacs5 = 'sim-gromacs5+'

        section = "Simulation"
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
        self.add_keys(section, self.sim_use_mpi, self.sim_time,
            self.sim_temperature, self.sim_weights, self.sim_fixed_weights,
            self.sim_weight_values, self.sim_fep_values,
            self.sim_vdw_values, self.sim_coul_values,
            self.sim_incrementor, self.sim_init_lambda,
            self.sim_fixed_lambda, self.sim_use_gibbs, self.sim_use_metro,
            self.sim_gibbs_delta, self.sim_nstout,
            self.sim_nst_mc, self.sim_pressure, self.sim_precision)
        ################################################
        # MDP Array
        self.mdr_count = 'mdr-number-of-simulations'
        self.mdr_threads = 'mdr-number-of-threads'
        self.mdr_queue_time = 'mdr-queue-time'
        self.mdr_genseed = 'mdr-generate-seeds'
        self._expected_frames = "mdr-expected-number-of-frames"

        section = "mdrun Array"
        self.add_key(self.mdr_threads, section=section, 'Max 20 threads (for now)')
        self.add_keys(section, self.mdr_count, self.mdr_queue_time,
            self.mdr_genseed)
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
        self.add_key(self.subcommand, section=section, updater=_subcomment)
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
    options[KEYS.sim_gromacs5] = False
    options[KEYS.sim_use_mpi] = False
    options[KEYS.sim_time] = 0.2  # ns
    options[KEYS.sim_temperature] = 300.0  # K
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
    options[KEYS.sim_use_metro] = True
    options[KEYS.sim_gibbs_delta] = -1
    options[KEYS.sim_nstout] = 500
    options[KEYS.sim_nst_mc] = 50
    options[KEYS.sim_temp_alg] = "v-rescale"  # "nose-hoover"
    options[KEYS.sim_pressure] = True
    options[KEYS.sim_precision] = True
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

class FileScan:

    TMPNAME = "#FilescanTempFile.out#"

    def __init__(self, path):
        self.filepath = path
        self.oldscan = []
        self.newscan = []
        self.weights = None
        self.num_of_steps = None

    def scan(self):
        if not os.path.exists(self.filepath):
            raise Exception('Simulation file not found ' + self.filepath)
        os.system('tail -500 {0} &>{1}'.format(self.filepath, FileScan.TMPNAME))
        if not os.path.exists(FileScan.TMPNAME):
            raise Exception('Tail command appears to have failed.')

        for line in open(FileScan.TMPNAME):
            if 'Step' in line and 'Time' in line and 'Lambda' in line:
                self.oldscan = self.newscan
                self.newscan = []
            self.newscan.append(line)
        os.remove(FileScan.TMPNAME)

        # assume a short scan indicates the file was cut off
        #     which could get fooled by dumping the state matrix, but it's only
        #     one frame off
        if len(self.oldscan) > len(self.newscan):
            self.newscan = self.oldscan

        count = 0
        capture_weights = False
        for line in self.scan:
            count += 1
            if count == 2:
                splt = line.split()
                self.num_of_steps = int(splt[0])
            if line.isspace():
                capture_weights = False
            if capture_weights:
                splt = line.split()
                weight = float(splt[5])
                self.weights.append(weight)
            if 'N' in line and 'Count' in line and 'G(in kT)' in line:
                capture_weights = True


    def get_wanglandau_weights(self):
        return self.weights

    def get_frame_number(self):
        return self.num_of_steps

class SlurmGen:

    def __init__(self):
        self.script = ""
        self.header = ""
        self.setup = ""
        self.checkin = ""
        self.main = ""
        self.cleanup = ""
        self.checkout = ""
        self.fields_general = dict()
        self.fields_header = dict()
        self.fields_checkin = dict()
        self.fields_setup = dict()
        self.fields_main = dict()
        self.fields_cleanup = dict()
        self.fields_checkout = dict()
        self.double_precision = False
        self.use_mpi = False
        self.gromacs5 = False
        self.calc_wt = True

    def walltime(self, ns, ntasks, nnodes):
        rate = 7.5 / 20  # 7.5 ns/day per 20 cores for double precicion
        if not self.double_precision:
            rate *= 1.8  # assume 80% speedup for single precision
        speed = rate * ntasks  # ns/day
        est = ns / speed  # days
        time_ = est * 1.20  # add 20% buffer
        time_limit = 7  # seven day limit on serial queue
        if nnodes > 1:
            time_limit = 2  # two day limit on parallel queue
        if time_ > time_limit:
            time_ = time_limit
        frac, days = math.modf(time_)
        frac, hours = math.modf(frac * 24.0)
        frac, mins = math.modf(frac * 60.0)
        frac, secs = math.modf(frac * 60.0)

        days = int(days)
        hours = int(hours)
        mins = int(mins)
        secs = int(secs)
        return '{0}-{1:02}:{2:02}:{3:02}'.format(days, hours, mins, secs)

    def compile(self):
        self.script = ""
        self.script += self.get_header()
        self.script += "\n"
        self.script += self.get_setup()
        self.script += "\n"
        self.script += self.get_checkin()
        self.script += "\n"
        self.script += self.get_main()
        self.script += "\n"
        self.script += self.get_cleanup()
        self.script += "\n"
        self.script += self.get_checkout()
        return self.script

    def get_header(self):
        fields = dict()
        fields.update(self.fields_general)
        fields.update(self.fields_header)
        fields["nnodes"] = 1
        fields["partition"] = "serial"
        if self.use_mpi:
            # 20 cores physically exist on a Rivanna node
            fields["nnodes"] = math.ceil(fields["ntasks"] / 20)
        if fields["nnodes"] > 1:
            fields["partition"] = "parallel"
        if self.calc_wt:
            fields["queue-time"] = self.walltime(fields["time-ns"],
                fields["ntasks"], fields["nnodes"])
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

    def get_setup(self):
        fields = dict()
        fields.update(self.fields_general)
        fields.update(self.fields_setup)
        self.setup = """MODULEPATH=/h3/n1/shirtsgroup/modules:$MODULEPATH
MODULEPATH=$HOME/modules:$MODULEPATH
module load jtp4kc
module load python
module load openmpi/1.8.1
"""
        if self.gromacs5:
            self.setup += "module load gromacs/5.0.2-sse"
        else:
            self.setup += "module load gromacs-shirtsgroup/4.6.7"

        self.setup += """

export THREADINFO="-nt {ntasks} "
export GMX_SUPPRESS_DUMP=1 #prevent step file output

sleep 1
diagnostics
"""
        self.setup = self.setup.format(**fields)
        return self.setup

    def get_checkin(self):
        fields = dict()
        fields.update(self.fields_general)
        fields.update(self.fields_checkin)
        self.checkin = """#Change the job status to 'SUBMITTED'
echo "SUBMITTED {job-name}{suffix} `date`" >> {jobstatus}
echo "Job Started at `date`"
""".format(**fields)
        return self.checkin

    def get_main(self):
        fields = dict()
        fields.update(self.fields_general)
        fields.update(self.fields_main)
        fields['dm'] = ""
        if self.double_precision:
            fields['dm'] = "_d"
        self.main = """#PRODUCTION
grompp{dm} -c {infolder}{gro-in}-in.gro -p {infolder}{base}.top -n {infolder}{base}.ndx -f {mdp-in}.mdp -o {job-name}.tpr -maxwarn 15
if [ -f {job-name}.tpr ];
then
    echo "MD Simulation launching..."
    mdrun{dm} ${{THREADINFO}} -deffnm {job-name}
else
    echo "Error generating {job-name}.tpr"
fi
""".format(**fields)
        return self.main

    def get_cleanup(self):
        fields = dict()
        fields.update(self.fields_general)
        fields.update(self.fields_cleanup)
        self.cleanup = ""
        if self.fields_cleanup["move-cmds"]:
            self.cleanup += "{move-cmds}\n"
        if self.fields_cleanup["extra-instructions"]:
            self.cleanup += "{extra-instructions}\n"
        self.cleanup = self.cleanup.format(**fields)
        return self.cleanup

    def get_checkout(self):
        fields = dict()
        fields.update(self.fields_general)
        fields.update(self.fields_checkout)
        self.checkout = """
echo "FINISHED {job-name}{suffix} `date`" >> {jobstatus}
# print end time
echo
echo "Job Ended at `date`"
echo "###################################################################"
""".format(**fields)
        return self.checkout

class MDPGen:

    class PkgPrecision:
        def __init__(self):
            self.tau_t = 1.0
            self.tau_p = 5.0  # ps
            self.shake_tol = 5e-6

    def __init__(self):
        self.script = ""
        self.params = ""
        self.neighbors = ""
        self.output = ""
        self.interactions = ""
        self.bonds = ""
        self.velocities = ""
        self.coupling = ""
        self.expdens = ""
        self.fields_general = dict()
        self.fields_params = dict()
        self.fields_neighbors = dict()
        self.fields_output = dict()
        self.fields_interactions = dict()
        self.fields_bonds = dict()
        self.fields_velocities = dict()
        self.fields_coupling = dict()
        self.fields_expdens = dict()
        self.ensemble = "NVT"
        self.volume_control = False
        self.fixed_state = False
        self.initial_weights = False
        self.use_gibbs = False
        self.use_metro = True
        self.pkg_single = MDPGen.PkgPrecision()
        self.pkg_single.tau_t = 1.0
        self.pkg_single.tau_p = 5.0
        self.pkg_single.shake_tol = 5e-6
        self.pkg_double = MDPGen.PkgPrecision()
        self.pkg_double.tau_t = 0.5
        self.pkg_double.tau_p = 20.0
        self.pkg_double.shake_tol = 1e-12
        self.precision = self.pkg_double
        self.gromacs5 = False

    def compile(self):
        self.script = ""
        self.script += self.get_params()
        self.script += "\n"
        self.script += self.get_neighbors()
        self.script += "\n"
        self.script += self.get_output()
        self.script += "\n"
        self.script += self.get_interactions()
        self.script += "\n"
        self.script += self.get_bonds()
        self.script += "\n"
        self.script += self.get_velocities()
        self.script += "\n"
        self.script += self.get_coupling()
        self.script += "\n"
        self.script += self.get_expdens()
        return self.script

    def set_ensemble(self, code="NVT"):
        self.ensemble = code
        if code == "NPT":
            self.volume_control = False
        else:
            self.volume_control = True

    def set_precision(self, code="double"):
        if code == "double":
            self.precision = self.pkg_double
        else:
            self.precision = self.pkg_single

    def get_params(self):
        fields = dict()
        fields.update(self.fields_general)
        fields.update(self.fields_params)
        fields['nstrun'] = fields['time-ns'] * (1000 / fields['dt'])
        self.params = """; RUN CONTROL PARAMETERS = 
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
""".format(**fields)
        return self.params

    def get_neighbors(self):
        fields = dict()
        fields.update(self.fields_general)
        fields.update(self.fields_neighbors)
        self.neighbors = """; NEIGHBORSEARCHING PARAMETERS = 
; nblist update frequency = 
nstlist                  = 10
; ns algorithm (simple or grid) = 
ns_type                  = grid
; Periodic boundary conditions: xyz or no = 
pbc                      = xyz
; nblist cut-off         = 
rlist                    = 1.2
""".format(**fields)
        return self.neighbors

    def get_output(self):
        fields = dict()
        fields.update(self.fields_general)
        fields.update(self.fields_output)
        self.output = """; OUTPUT CONTROL OPTIONS
; Output frequency for coords (x), velocities (v) and forces (f)
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
; Output frequency for energies to log file and energy file
nstlog                   = {nstout:0.0f}  ; changing this allows you to see the frequency the weights are computed.
nstcalcenergy            = 1
nstenergy                = {nstout:0.0f}
; Output frequency and precision for .xtc file"""
        if self.gromacs5:
            self.output += """
nstxout-compressed       = {nstout:0.0f}  ; trajectory output frequency
compressed-x-precision   = 1000
"""
        else:
            self.output += """
nstxtcout                = {nstout:0.0f}  ; trajectory output frequency
xtc-precision            = 1000
"""
            self.output = self.output.format(**fields)
        return self.output

    def get_interactions(self):
        fields = dict()
        fields.update(self.fields_general)
        fields.update(self.fields_interactions)
        self.interactions = """; OPTIONS FOR ELECTROSTATICS AND VDW = 
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
""".format(**fields)
        return self.interactions

    def get_bonds(self):
        fields = dict()
        fields.update(self.fields_general)
        fields.update(self.fields_bonds)
        fields.update(self.precision.__dict__)
        self.bonds = """; OPTIONS FOR BONDS     = 
constraints              = hbonds ; constrain bonds to hydrogen
; Type of constraint algorithm = 
constraint-algorithm     = shake
; Highest order in the expansion of the constraint coupling matrix = 
shake-tol                = {shake_tol} ; ~TP
""".format(**fields)
        return self.bonds

    def get_velocities(self):
        fields = dict()
        fields.update(self.fields_general)
        fields.update(self.fields_velocities)
        self.velocities = """; GENERATE VELOCITIES FOR STARTUP RUN = 
gen_vel                  = yes
gen_temp                 = {temp:0.1f} ; K
gen_seed                 = {gen-seed}
""".format(**fields)
        return self.velocities

    def get_coupling(self):
        fields = dict()
        fields.update(self.fields_general)
        fields.update(self.fields_coupling)
        fields.update(self.precision.__dict__)
        self.coupling = """; OPTIONS FOR WEAK COUPLING ALGORITHMS = 
; Groups to couple separately = 
tc-grps                  = System
; Time constant (ps) and reference temperature (K) = 
tcoupl                   = {temp_alg}
tau_t                    = {tau_t} ; ~TP
ref_t                    = {temp:0.1f}
"""
        if self.volume_control:
            self.coupling += """
; Volume control         = 
Pcoupl                   = no
"""
        else:
            self.coupling += """
; Pressure coupling      = 
Pcoupl                   = MTTK
Pcoupltype               = isotropic
; Time constant (ps), compressibility (1/bar) and reference P (bar) = 
tau_p                    = {tau_p} ; ps ~TP
compressibility          = 4.5e-5 ; 1/bar
ref_p                    = 1.0 ; bar
"""
        self.coupling = self.coupling.format(**fields)
        return self.coupling

    def get_expdens(self):
        fields = dict()
        fields.update(self.fields_general)
        fields.update(self.fields_expdens)
        self.expdens = "; OPTIONS FOR EXPANDED ENSEMBLE SIMULATIONS"
        self.expdens += self.get_free_energy()
        self.expdens += "\n"
        self.expdens += self.get_expanded_ensemble()
        self.expdens += "\n"
        self.expdens += self.get_monte_carlo()
        return self.expdens

    def get_free_energy(self):
        fields = dict()
        fields.update(self.fields_general)
        fields.update(self.fields_expdens)
        fields['energy-out'] = "yes"
        if self.gromacs5:
            fields['energy-out'] = "total"
        part = """; Free energy control stuff  
sc-alpha                  = 0.5
couple-moltype            = {ligand}
couple-lambda0            = vdw-q
couple-lambda1            = none
couple-intramol           = no
fep-lambdas               = {fep-lambdas}
coul-lambdas              = {coul-lambdas}
vdw-lambdas               = {vdw-lambdas}
symmetrized-transition-matrix = yes
nst-transition-matrix     = 100000
nstdhdl                   = {nstout:0.0f}
dhdl-print-energy         = {energy-out}
"""
        if self.fixed_state:
            part += "free-energy               = yes\n"
            part += "init-lambda               = {init-lambda}\n"
        else:
            part += "free-energy               = expanded\n"
            part += "init-lambda-state         = {init-state-index}\n"
        part = part.format(**fields)
        return part

    def get_expanded_ensemble(self):
        fields = dict()
        fields.update(self.fields_general)
        fields.update(self.fields_expdens)
        part = """; expanded ensemble
nstexpanded              = {nst-mc:0.0f}
lmc-stats                = wang-landau
lmc-weights-equil        = {weights-equil}
lmc-seed                 = {lmc-seed}
"""
        if self.use_gibbs and self.use_metro:
            part += "lmc-move                 = metropolized-gibbs\n"
            part += "lmc-gibbsdelta           = {gibbs-delta}\n"
        elif self.use_gibbs:
            part += "lmc-move                 = gibbs\n"
            part += "lmc-gibbsdelta           = {gibbs-delta}\n"
        elif self.use_metro:
            part += "lmc-move                 = metropolis\n"
        else:
            part += "lmc-move                 = no\n"

        if self.initial_weights:
            part += "init-lambda-weights      = {wl-weights}\n"
        else:
            part += "weight-equil-wl-delta    = 0.001\n"
        part = part.format(**fields)
        return part

    def get_monte_carlo(self):
        fields = dict()
        fields.update(self.fields_general)
        fields.update(self.fields_expdens)
        part = """; Seed for Monte Carlo in lambda space
wl-scale                 = 0.5
wl-ratio                 = 0.8
init-wl-delta            = {incrementor:0.3f}
wl-oneovert              = yes
""".format(**fields)
        return part

class BEPGen:

    class MyBEPKeys(Keys):
        def __init__(self):
            Keys.__init__(self)
            # Values that can be specified by the user in the options file
            ################################################
            # Some parameters and definitions.
            self.test_directory = 'test_directory'
            self.base_func_dir = 'base_func_dir'
            self.test_job_name = 'test_job_name'
            self.number_of_states = 'number_of_states'
            self.manual_num_states = 'manual_num_states'
            self.coupling_nums = 'coupling_nums'
            self.ligplot_src = 'ligplot_src'
            self.hbdir = 'hbdir'
            self.iterations_info = 'iterations_info'
            self.last_states_coupled = 'last_states_coupled'
            self.temperature = 'temperature'

            section = "Parameters"
            self.add_keys(section, self.test_directory, self.base_func_dir,
                self.test_job_name, self.manual_num_states, self.ligplot_src,
                self.iterations_info, self.last_states_coupled,
                self.temperature, self.coupling_nums, self.hbdir,
                self.number_of_states)
            ################################################
            # Process Control
            self.verbosity = 'verbose-level'

            section = "Process Control"
            self.add_keys(section, self.verbosity)
            ################################################
            # Gromacs Analysis
            self.ligand_res = 'ligand-residue-name'
            self.gro_files = 'gro_initial-frames'
            self.xtc_files = 'gro_trajectories'
            self.tpr_files = 'gro_run-input-files'
            self.edr_files = 'gro_energy-files'
            self.dhdl_files = 'gro_dhdl-files'
            self.single_state = 'single-state-simulation_state-index'

            section = "Gromacs Analysis"
            self.add_keys(section, self.gro_files, self.xtc_files,
                self.tpr_files, self.edr_files, self.dhdl_files,
                self.single_state, self.ligand_res)
            ################################################
            # Preliminary_analysis
            self.check_average_energies = 'check_average_energies'

            section = "Preliminary_analysis"
            self.add_keys(section, self.check_average_energies)

    KEYS = MyBEPKeys()
    params = Parameters(KEYS)

    def __init__(self):
        self.fields_output = dict()

    def write(self, file_name):
        BEPGen.params.write_options(file_name, self.fields_output)

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
    dir_ = os.path.join(backup.expandpath(opts[KEYS.script_dir]), "")
    if not os.path.exists(dir_):
        os.mkdir(dir_)
    os.chdir(dir_)
    path = os.path.join(backup.expandpath(opts[KEYS.work_dir]), "")

    tpr_files = []
    xtc_files = []
    gro_files = []
    edr_files = []
    xvg_files = []

    for i in range(opts[KEYS.mdr_count]):
        suffix = '_{0:0>2}'.format(i)
        if opts[KEYS.mdr_count] < 2:
            suffix = ''

        file_name = job_name + suffix + '.slurm'
        timens = opts[KEYS.sim_time]
        workdir = os.path.join(path, job_name + suffix, "")
        callingdir = opts[KEYS._calling_dir]
        indir = os.path.join('..', '')
        gro_in = base_name + suffix
        mdp_in = job_name + suffix
        move_par = True
        move_slurm = True
        move_gro = False

        cont = False
        if randseed:
            gro_in = base_name
            try:
                import random
                # random.seed(0)
                seed = random.randint(0, 65536)
                num = seed + i
                cont = make_mdp(opts, name=mdp_in, genseed=num, lmcseed=num)
            except IOError as ioex:
                print(ioex)
        else:
            try:
                cont = make_mdp(opts, name=job_name + suffix)
            except IOError as ioex:
                print(ioex)

        if not cont:
            return

        move = ""
        if move_par:
            move += ("mv " + dir_ + opts[KEYS._params_out] + " " +
                opts[KEYS._params_out] + " \n")
        if move_slurm:
            move += "mv " + dir_ + file_name + " " + file_name + "\n"
        if move_gro:
            move += "mv " + indir + gro_in + "-in.gro " + gro_in + "-in.gro\n"

        extra = ""
        subcom = opts[KEYS.subcommand]
        if opts[KEYS._chain_all] and (subcom == 'equil' or
            subcom == 'rand' or subcom == 'array'):

            if subcom == 'equil':
                next_cmd = 'rand'
            elif subcom == 'rand':
                next_cmd = 'array'
            elif subcom == 'array':
                move_gro = True
            param_opt = ""
            if opts[KEYS._params]:
                param_opt = "--par " + opts[KEYS._params]

            extra = "# genscript " + subcom + "\n"
            if subcom == 'array':
                extra += "cd " + workdir + "\n"
            else:
                extra += ("cd " + callingdir + "\n" +
                    "python " + os.path.realpath(__file__) + " " + next_cmd +
                    " " + KEYS._chain_all + " --submit " + param_opt + "\n" +
                    "cd " + workdir + "\n")

        builder = SlurmGen()
        builder.double_precision = opts[KEYS.sim_precision]
        builder.use_mpi = opts[KEYS.sim_use_mpi]
        builder.calc_wt = opts[KEYS.auto_calc_wt]
        builder.gromacs5 = opts[KEYS.sim_gromacs5]
        builder.fields_general['job-name'] = job_name
        builder.fields_general['suffix'] = suffix
        builder.fields_general['ntasks'] = opts[KEYS.mdr_threads]
        builder.fields_general['workdir'] = workdir
        builder.fields_general['jobstatus'] = os.path.join('..',
            'jobstatus.txt')
        builder.fields_header['queue-time'] = opts[KEYS.mdr_queue_time]
        builder.fields_header['outputfile'] = os.path.join(path,
            job_name + suffix, job_name)
        builder.fields_main['base'] = base_name
        builder.fields_main['gro-in'] = gro_in
        builder.fields_main['mdp-in'] = mdp_in
        builder.fields_main['infolder'] = indir
        builder.fields_cleanup['move-cmds'] = move
        builder.fields_cleanup['extra-instructions'] = extra

        file_ = open(file_name, 'w')
        file_.write(builder.compile())
        file_.close()
        if submit:
            folder = job_name + suffix
            dir_name = os.path.join(path, folder)
            if not os.path.exists(dir_name):
                os.mkdir(dir_name)
            # yes, this seems to be unnecessary, but it helps in debugging if
            #    directory doesn't actually have to exist to run without submit
            os.system("mv " + mdp_in + ".mdp " +
                os.path.join(dir_name, mdp_in) + ".mdp")
            tpr_files.append(os.path.join(folder, job_name + ".tpr"))
            xtc_files.append(os.path.join(folder, job_name + ".xtc"))
            if move_gro:
                gro_files.append(os.path.join(folder, gro_in + ".gro"))
            else:
                gro_files.append(os.path.join(gro_in + ".gro"))
            edr_files.append(os.path.join(folder, job_name + ".edr"))
            xvg_files.append(os.path.join(folder, job_name + ".xvg"))
            if not opts[KEYS._dryrun]:
                os.system("sbatch " + file_name)
                print("Sbatch'd Job " + job_name + suffix)
            else:
                print("DRYRUN: Would sbatch job " + job_name + suffix)

    if submit:
        num_coup = opts[KEYS.sim_genxcoupled]
        num_uncp = opts[KEYS.sim_genxuncupld]
        num_states = len(opts[KEYS.sim_fep_values]) + num_coup + num_uncp
        num_coup += 1  # assuming fep has one coupled
        num_uncp += 1  # assuming fep has one uncoupled
        num_frames = int(opts[KEYS._expected_frames] + 1)
        if opts[KEYS.sim_init_lambda] >= 0:
            state_index = opts[KEYS.sim_init_lambda]
        else:
            state_index = num_states + opts[KEYS.sim_init_lambda]

        analysis = BEPGen()
        analysis.fields_output[BEPGen.KEYS.test_directory] = "."
        analysis.fields_output[BEPGen.KEYS.base_func_dir] = ("~/Documents/" +
            "git/binding-ensemble-analysis")
        analysis.fields_output[BEPGen.KEYS.test_job_name] = job_name
        analysis.fields_output[BEPGen.KEYS.number_of_states] = num_states
        analysis.fields_output[BEPGen.KEYS.manual_num_states] = True
        analysis.fields_output[BEPGen.KEYS.coupling_nums] = [num_coup, num_uncp]
        analysis.fields_output[BEPGen.KEYS.ligplot_src] = ("~/Documents/" +
            "LIGPLOT/source/ligplot.scr")
        analysis.fields_output[BEPGen.KEYS.hbdir] = "~/Documents/LIGPLOT"
        analysis.fields_output[BEPGen.KEYS.iterations_info] = [num_frames, 0]
        analysis.fields_output[BEPGen.KEYS.last_states_coupled] = False
        analysis.fields_output[BEPGen.KEYS.temperature] = opts[KEYS.
            sim_temperature]
        analysis.fields_output[BEPGen.KEYS.verbosity] = 2
        analysis.fields_output[BEPGen.KEYS.ligand_res] = opts[KEYS.ligand]
        analysis.fields_output[BEPGen.KEYS.gro_files] = gro_files
        analysis.fields_output[BEPGen.KEYS.xtc_files] = xtc_files
        analysis.fields_output[BEPGen.KEYS.tpr_files] = tpr_files
        analysis.fields_output[BEPGen.KEYS.edr_files] = edr_files
        analysis.fields_output[BEPGen.KEYS.dhdl_files] = xvg_files
        analysis.fields_output[BEPGen.KEYS.check_average_energies] = False
        if opts[KEYS.sim_fixed_lambda]:
            analysis.fields_output[BEPGen.KEYS.single_state] = state_index
        analysis.write("analysis_" + opts[KEYS.job_name] + ".bep")

    os.chdir(cur_dir)

def make_mdp(opts, dir_='.', name=None, genseed=10200, lmcseed=10200):
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

    state_index = 0
    if opts[KEYS.sim_init_lambda] >= 0:
        state_index = opts[KEYS.sim_init_lambda]
    else:
        state_index = len(fep) + opts[KEYS.sim_init_lambda]

    fep_lambdas = ''
    coul_lambdas = ''
    vdw_lambdas = ''
    for f, c, v in zip(fep, coul, vdw):
        fep_lambdas += '{0:0.1f} '.format(f)
        coul_lambdas += '{0:0.1f} '.format(c)
        vdw_lambdas += '{0:0.1f} '.format(v)

    wl_weights = ''
    if opts[KEYS.sim_weights]:
        for w in weights:
            wl_weights += '{0:0.5f} '.format(w)
        if not len(fep) == len(weights):
            print('Number of weights not equal to number of states, please' +
                ' address')
            return False

    equil = ""
    if opts[KEYS.sim_fixed_weights] and opts[KEYS.sim_weights]:
        equil = 'yes'
    else:
        equil = 'wl-delta'

    init_lambda = (state_index + 1) / len(fep)
    if init_lambda < 0:
        init_lambda = 0
    elif init_lambda > 1:
        init_lambda = 1

    steps = opts[KEYS.sim_time] * (1000 / 0.002)
    opts[KEYS._expected_frames] = math.ceil(steps / opts[KEYS.sim_nstout])

    builder = MDPGen()
    builder.fixed_state = opts[KEYS.sim_fixed_lambda]
    builder.initial_weights = opts[KEYS.sim_weights]
    builder.gromacs5 = opts[KEYS.sim_gromacs5]
    builder.use_metro = opts[KEYS.sim_use_metro]
    builder.use_gibbs = opts[KEYS.sim_use_gibbs]
    if opts[KEYS.sim_pressure]:
        builder.set_ensemble("NPT")
    else:
        builder.set_ensemble("NVT")
    if opts[KEYS.sim_precision]:
        builder.set_precision("double")
    else:
        builder.set_precision("single")
    builder.fields_general['nstout'] = opts[KEYS.sim_nstout]
    builder.fields_general['temp'] = opts[KEYS.sim_temperature]
    builder.fields_params['dt'] = 0.002;  # ps
    builder.fields_params['time-ns'] = opts[KEYS.sim_time]
    builder.fields_velocities['gen-seed'] = genseed
    builder.fields_expdens['ligand'] = opts[KEYS.ligand]
    builder.fields_expdens['fep-lambdas'] = fep_lambdas
    builder.fields_expdens['coul-lambdas'] = coul_lambdas
    builder.fields_expdens['vdw-lambdas'] = vdw_lambdas
    builder.fields_expdens['init-lambda'] = init_lambda
    builder.fields_expdens['init-state-index'] = state_index
    builder.fields_expdens['nst-mc'] = opts[KEYS.sim_nst_mc]
    builder.fields_expdens['weights-equil'] = equil
    builder.fields_expdens['lmc-seed'] = lmcseed
    builder.fields_expdens['gibbs-delta'] = int(opts[KEYS.sim_gibbs_delta])
    builder.fields_expdens['wl-weights'] = wl_weights
    builder.fields_expdens['incrementor'] = opts[KEYS.sim_incrementor]

    file_name = name + '.mdp'
    file_ = open(file_name, 'w')
    file_.write(builder.compile())
    file_.close()

    os.chdir(cur_dir)
    return True

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

            logscan = FileScan(file_name)
            logscan.scan()
            weights = logscan.get_wanglandau_weights()
            opts[KEYS.sim_weight_values] = weights

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

            logscan = FileScan(file_name)
            logscan.scan()
            num_of_steps = logscan.get_frame_number()
            num_frames = math.floor(num_of_steps / opts[KEYS.sim_nstout])

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
        make_mdp(opts)
    if opts[KEYS.run_array]:
        generate(opts)

def setup(options, args, opts, parser, cur_dir):
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

    opts[KEYS._calling_dir] = os.path.join(cur_dir, '')
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
        return 0

    (options, args) = parser.parse_args(argv)

    # handle options, reading from file if requested
    if options.par:
        print('Reading parameters from ' + options.par)
        opt_list = param.parse_options(options.par, opts)
        if not opt_list:
            print('Error reading from parameters file.')
            return 1
    else:
        opt_list = opts

    if isinstance(opt_list, list):
        for opts in opt_list:
            setup(options, args, opts, parser, cur_dir)
    else:
        setup(options, args, opt_list, parser, cur_dir)

SUBS.update({'exit': gen_exit, 'all': gen_all, 'equil': gen_equil,
    'rand': gen_rand, 'array': gen_array, 'mdp': make_mdp,
    'opt': gen_opt})
KEYS.update()
if __name__ == '__main__':
    # initialize the subcommand list
    main(sys.argv[1:])



