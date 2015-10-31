#!/usr/bin/env python
'''
Created on Jun 19, 2015

@author: Tyler P
'''
import sys, os
import traceback
from optparse import OptionParser
import math
import backup
from param_versions.version_2_0 import Parameters
from param_versions.version_2_0 import Keys
import job_utils
import mdp_template
import slurm_template

verbose = 0
POST_COMMANDS = []  #
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
        self.input_dir = 'input-file-dir'
        self.ligand = 'ligand-res-name'

        section = "General"
        self.add_key(self.base_name, section, "Input file prefix")
        self.add_key(self.job_name, section, "Name to prefix all files")
        self.add_key(self.work_dir, section, "Directory to operate within," +
            " and create save files. Should generally be the same as " +
            self.script_dir)
        self.add_key(self.script_dir, section, "Directory to create scripts" +
            " and script folders within")
        self.add_key(self.input_dir, section, "Location of input files (.top," +
            " .ndx, -in.gro)")
        self.add_keys(section, self.ligand)
        ################################################
        # Sim Options
        self.sim_use_mpi = "sim-use-mpi-parallelization"
        self.sim_time = 'sim-time-ns'
        self.sim_dt = 'sim-time-step'
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
        self.sim_wl_scale = 'sim-weight-scale'
        self.sim_wl_ratio = 'sim-weight-ratio'
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
            " states by subtracting ln(x) from all but the beginning state"))
        self.add_key(self.sim_wgtxuncupld, section, ("Emulate these number of" +
            " states by adding ln(x) to just the end state"))
        self.add_keys(section, self.sim_use_mpi, self.sim_time, self.sim_dt,
            self.sim_temperature, self.sim_weights, self.sim_fixed_weights,
            self.sim_weight_values, self.sim_fep_values,
            self.sim_vdw_values, self.sim_coul_values,
            self.sim_incrementor, self.sim_wl_scale,
            self.sim_wl_ratio, self.sim_init_lambda,
            self.sim_fixed_lambda, self.sim_use_gibbs, self.sim_use_metro,
            self.sim_gibbs_delta, self.sim_nstout,
            self.sim_nst_mc, self.sim_temp_alg, self.sim_pressure,
            self.sim_precision, self.sim_gromacs5)
        ################################################
        # MDP Array
        self.mdr_count = 'mdr-number-of-simulations'
        self.mdr_threads = 'mdr-number-of-threads'
        self.mdr_queue_time = 'mdr-queue-time'
        self.mdr_genseed = 'mdr-generate-seeds'
        self.mdr_seedrand = 'mdr-seed-for-random'
        self._expected_frames = "mdr-expected-number-of-frames"

        section = "mdrun Array"
        self.add_key(self.mdr_threads, section, 'Max 20 threads (for now)')
        self.add_keys(section, self.mdr_count, self.mdr_queue_time,
            self.mdr_genseed, self.mdr_seedrand)
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
    options[KEYS.sim_dt] = 0.002  # ps
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
    options[KEYS.sim_wl_ratio] = 0.8
    options[KEYS.sim_wl_scale] = 0.8
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
    options[KEYS.mdr_seedrand] = None
    ################################################
    # Automation
    options[KEYS.auto_time_equil] = 0.2  # ns
    options[KEYS.auto_time_rand] = 0.2  # ns
    options[KEYS.auto_time_array] = 0.2  # ns
    options[KEYS.auto_calc_wt] = True

    return options

# provide custom default library
param.option_defaults = option_defaults

saver = job_utils.SaveJobs()

#################################################################################
#################################################################################

class FilesystemImpactRegister:

    def __init__(self):
        self.cwd = '.'
        self.files_created = []
        self.filename_bases = []
        self.files_removed = []

    def merge(self, other_FIR):
        fir = FilesystemImpactRegister()
        fir.files_created.extend(self.files_created)
        fir.files_created.extend(other_FIR.files_created)
        fir.filename_bases.extend(self.filename_bases)
        fir.filename_bases.extend(other_FIR.filename_bases)
        fir.files_removed.extend(self.files_removed)
        fir.files_removed.extend(other_FIR.files_removed)
        return fir

class SlurmGen(slurm_template.Slurm):

    def __init__(self):
        slurm_template.Slurm.__init__(self)
        self.double_precision = False
        self.use_mpi = False
        self.gromacs5 = False
        self.calc_wt = True
        self.extra = None

        self.opt_jobname = ""
        self.opt_suffix = ""
        self.opt_ntasks = 1
        self.opt_workdir = None
        self.opt_queuetime = "2-00:00:00"
        self.opt_timens = 1

        self.file_out = ""
        self.file_mdp = ""
        self.file_gro = ""
        self.file_top = ""
        self.file_ndx = ""
        self.file_tpr = ""

    def walltime(self, ns, ntasks, nnodes):
        rate = 7.5 / 20  # 7.5 ns/day per 20 cores for double precicion
        if not self.double_precision:
            rate *= 1.7  # assume 70% speedup for single precision
        speed = rate * ntasks  # ns/day
        est = ns / speed  # days
        time_ = est * 1.20  # add 20% buffer
        time_limit = 7  # seven day limit on serial queue
        if nnodes > 1:
            time_limit = 2  # two day limit on parallel queue
        if time_ > time_limit:
            time_ = time_limit
        if time_ < 0.5:  # Since it doesn't hurt, and algorithm isn't perfect
            time_ = 0.5  # make this a minimum run time, to get everything
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
        self.set_header()
        self.set_setup()
        self.set_commands()
        return slurm_template.Slurm.compile(self)

    def set_header(self):
        nnodes = 1
        partition = "serial"
        queuetime = self.opt_queuetime

        if self.use_mpi:
            # 20 cores physically exist on a Rivanna node
            nnodes = int(math.ceil(self.opt_ntasks / 20.0))
        if nnodes == 0:
            nnodes = 1
        elif nnodes > 1:
            partition = "parallel"
        if self.calc_wt:
            queuetime = self.walltime(self.opt_timens, self.opt_ntasks, nnodes)

        self.job_name = self.opt_jobname + self.opt_suffix
        self.partition = partition
        self.nodes = nnodes
        self.ntasks = self.opt_ntasks
        self.time = queuetime
        self.signal = 15  # 15 = SIGTERM
        self.output = self.file_out
        self.workdir = self.opt_workdir
        self.mail_types.append("REQUEUE")
        self.mail_types.append("END")
        self.mail_types.append("FAIL")
        self.mail_user = "jtp4kc@virginia.edu"

    def set_setup(self):
        self.status_filename = os.path.join('..', 'jobstatus.txt')

        self.modulepath_prepends.append("/h3/n1/shirtsgroup/modules")
        self.modulepath_prepends.append("$HOME/modules")
        self.modules.append("jtp4kc")
        self.modules.append("anaconda-jtp")
        if self.gromacs5:
            self.modules.append("module load gromacs/5.0.2-sse")
        else:
            self.modules.append("gromacs-shirtsgroup/4.6.7")

        self.exports.append('THREADINFO="-nt ' + str(self.opt_ntasks) + ' "')
        self.exports.append('GMX_SUPPRESS_DUMP=1 #prevent step file output')

    def set_commands(self):
        dm = ""
        if self.double_precision:
            dm = "_d"

        grompp = ["grompp" + dm,
                  "-c", self.file_gro,
                  "-p", self.file_top,
                  "-n", self.file_ndx,
                  "-f", self.file_mdp,
                  "-o", self.file_tpr,
                  "-maxwarn", "15"]

        self.commands.append("rm " + self.file_tpr + " # in case files changed")
        self.commands.append(" ".join(grompp))
        self.commands.append("if [ -f " + self.file_tpr + " ];")
        self.commands.append("then")
        self.commands.append('    echo "MD Simulation launching..."')
        self.commands.append("    mdrun" + dm + " ${THREADINFO} -deffnm " +
                             self.opt_jobname)
        self.commands.append("else")
        self.commands.append('    echo "Error creating ' + self.file_tpr + '"')
        self.commands.append("fi")

        if self.extra:
            self.commands.append("")
            self.commands.append(self.extra)

class MDPGen(mdp_template.MDP):

    class PkgPrecision:
        def __init__(self):
            self.tau_t = 1.0
            self.tau_p = 5.0  # ps
            self.shake_tol = 5e-6

    def __init__(self):
        mdp_template.MDP.__init__(self)
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

        self.opt_nstout = 0
        self.opt_temp = 0
        self.opt_dt = 0
        self.opt_timens = 0
        self.opt_genseed = 0
        self.opt_tempalg = ""
        self.opt_ligand = ""
        self.opt_feplambdas = ""
        self.opt_coullambdas = ""
        self.opt_vdwlambdas = ""
        self.opt_initlambda = 0
        self.opt_initlambdaindex = 0
        self.opt_nstmc = ""
        self.opt_weightsequil = ""
        self.opt_lmcseed = 0
        self.opt_gibbsdelta = ""
        self.opt_wlweights = ""
        self.opt_wlincr = 0
        self.opt_wlscale = 0
        self.opt_wlratio = 0

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

    def compile(self):
        self.set_params()
        self.set_neighbors()
        self.set_output()
        self.set_interactions()
        self.set_bonds()
        self.set_velocities()
        self.set_coupling()
        self.set_free_energy()
        self.set_expanded_ensemble()
        self.set_monte_carlo()
        return mdp_template.MDP.compile(self)

    def set_params(self):
        nstrun = self.opt_timens * (1000 / self.opt_dt)
        self.integrator = "md-vv"
        self.tinit = "0"
        self.dt = "{0:0.3f}".format(self.opt_dt)
        self.nsteps = "{0:<12.0f}; {1:0.2f} ns".format(nstrun, self.opt_timens)
        self.comm_mode = "Linear"
        self.nstcomm = "1"
        self.nsttcouple = "1"
        self.nstpcouple = "1"
        self.include = "-I/nv/blue/jtp4kc/gromacs/itp"

    def set_neighbors(self):
        self.nstlist = "10"
        self.ns_type = "grid"
        self.pbc = "xyz"
        self.rlist = "1.2"

    def set_output(self):
        self.nstxout = "0"
        self.nstvout = "0"
        self.nstfout = "0"
        self.nstlog = "{0:0.0f}".format(self.opt_nstout)
        self.nstcalcenergy = "1"
        self.nstenergy = self.nstlog
        if self.gromacs5:
            self.nstxout_compressed = self.nstlog
            self.compressed_x_precision = "1000"
        else:
            self.nstxtcout = self.nstlog
            self.xtc_precision = "1000"

    def set_interactions(self):
        self.cutoff_scheme = "group"
        self.coulombtype = "PME"
        self.coulomb_modifier = "Potential-Switch"
        self.rcoulomb_switch = "0.88"
        self.rcoulomb = "0.9"
        self.vdw_type = "Cut-off"
        self.vdw_modifier = "Potential-Switch"
        self.rvdw_switch = "0.85"
        self.rvdw = "0.9"
        self.DispCorr = "AllEnerPres"
        self.fourierspacing = "0.12"
        self.fourier_nx = "0"
        self.fourier_ny = "0"
        self.fourier_nz = "0"
        self.pme_order = "4"
        self.ewald_rtol = "1e-05"
        self.ewald_geometry = "3d"

    def set_bonds(self):
        self.constraints = "hbonds ; constrain bonds to hydrogen"
        self.constraint_algorithm = "shake"
        self.shake_tol = str(self.precision.shake_tol)

    def set_velocities(self):
        self.comment_velocities1 = "GENERATE VELOCITIES FOR STARTUP RUN = "
        self.gen_vel = "yes"
        self.gen_temp = "{0:0.1f} ; K".format(self.opt_temp)
        self.gen_seed = str(self.opt_genseed)

    def set_coupling(self):
        self.comment_coupling1 = "OPTIONS FOR WEAK COUPLING ALGORITMS"
        self.tc_grps = "System"
        self.tcoupl = self.opt_tempalg
        self.tau_t = str(self.precision.tau_t)
        self.ref_t = str(self.opt_temp)
        if self.volume_control:
            self.comment_coupling3 = "Volume control"
            self.pcoupl = "no"
        else:
            self.comment_coupling3 = "Pressure coupling"
            self.pcoupl = "MTTK"
            self.pcoupltype = "isotropic"
            self.tau_p = str(self.precision.tau_p) + "; ps"
            self.compressibility = "4.5e-5 ; 1/bar"
            self.ref_p = "1.0 ; bar"

    def set_free_energy(self):
        self.sc_alpha = "0.5"
        self.couple_moltype = self.opt_ligand
        self.couple_lambda0 = "vdw-q"
        self.couple_lambda1 = "none"
        self.couple_intramol = "no"
        self.fep_lambdas = self.opt_feplambdas
        self.coul_lambdas = self.opt_coullambdas
        self.vdw_lambdas = self.opt_vdwlambdas
        self.symmetrized_transition_matrix = "yes"
        self.nst_transition_matrix = "100000"
        self.nstdhdl = str(self.opt_nstout)
        if self.gromacs5:
            self.dhdl_print_energy = "total"
        else:
            self.dhdl_print_energy = "yes"
        if self.fixed_state:
            self.free_energy = "yes"
            self.init_lambda = str(self.opt_initlambda)
        else:
            self.free_energy = "expanded"
            self.init_lambda_state = str(self.opt_initlambdaindex)

    def set_expanded_ensemble(self):
        self.nstexpanded = "{0:0.0f}".format(self.opt_nstmc)
        self.lmc_stats = "wang-landau"
        self.lmc_weights_equil = str(self.opt_weightsequil)
        self.lmc_seed = str(self.opt_lmcseed)
        if self.use_gibbs and self.use_metro:
            self.lmc_move = "metropolized-gibbs"
            self.lmc_gibbsdelta = str(self.opt_gibbsdelta)
        elif self.use_gibbs:
            self.lmc_move = "gibbs"
            self.lmc_gibbsdelta = str(self.opt_gibbsdelta)
        elif self.use_metro:
            self.lmc_move = "metropolis"
        else:
            self.lmc_move = "no"
        if self.initial_weights:
            self.init_lambda_weights = str(self.opt_wlweights)
        else:
            self.weight_equil_wl_delta = "0.001"

    def set_monte_carlo(self):
        self.wl_scale = "{0:0.3f}".format(self.opt_wlscale)
        self.wl_ratio = "{0:0.3f}".format(self.opt_wlratio)
        self.init_wl_delta = "{0:0.3f}".format(self.opt_wlincr)
        self.wl_oneovert = "yes"

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
    fir = FilesystemImpactRegister()
    fir.cwd = os.getcwd()

    randseed = opts[KEYS.mdr_genseed]
    submit = opts[KEYS._submit]
    job_name = opts[KEYS.job_name]
    if not job_name:
        print('Job name not specified')
        return
    base_name = opts[KEYS.base_name]
    if not base_name:
        base_name = job_name

    here = os.path.join(os.curdir, "")
    dir_ = os.path.join(backup.expandrelpath(opts[KEYS.work_dir]), "")
    if dir_.startswith(here):
        dir_ = dir_.replace(here, "")
    if dir_:
        if not os.path.exists(dir_):
            os.mkdir(dir_)
        os.chdir(dir_)
    path = os.path.join(backup.expandrelpath(opts[KEYS.script_dir]), "")
    if path.startswith(here):
        path = path.replace(here, "")

    tpr_files = []
    xtc_files = []
    gro_files = []
    edr_files = []
    xvg_files = []

    for i in range(opts[KEYS.mdr_count]):
        suffix = '-{0:0>2}'.format(i)
        if opts[KEYS.mdr_count] < 2:
            suffix = ''

        folder = job_name + suffix
        dir_name = os.path.join(path, folder)
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)
        os.chdir(dir_name)
        file_name = job_name + suffix + '.slurm'
        workdir = os.path.join(path, job_name + suffix, "")
        callingdir = opts[KEYS._calling_dir]
        indir = os.path.join('..', '')
        gro_in = base_name + suffix
        mdp_in = job_name + suffix
        move_par = True
        move_gro = False

        cont = False
        if randseed:
            gro_in = base_name
            try:
                import random
                if opts[KEYS.seedrand] != None:
                    random.seed(opts[KEYS.seedrand])
                seed = random.randint(0, 65536)
                num = seed + i
                cont = make_mdp(opts, name=mdp_in, genseed=num, lmcseed=num)
            except IOError as ioex:
                print(ioex)
                traceback.print_exc()
        else:
            try:
                cont = make_mdp(opts, name=job_name + suffix)
            except IOError as ioex:
                print(ioex)
                traceback.print_exc()

        if not cont:
            return

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

        move = ""
        if move_par:
            move += ("mv " + dir_ + opts[KEYS._params_out] + " " +
                opts[KEYS._params_out] + " \n")
        if move_gro:
            move += "mv " + indir + gro_in + "-in.gro " + gro_in + "-in.gro\n"

        builder = SlurmGen()
        builder.double_precision = opts[KEYS.sim_precision]
        builder.use_mpi = opts[KEYS.sim_use_mpi]
        builder.calc_wt = opts[KEYS.auto_calc_wt]
        builder.gromacs5 = opts[KEYS.sim_gromacs5]

        builder.opt_jobname = job_name
        builder.opt_suffix = suffix
        builder.opt_ntasks = opts[KEYS.mdr_threads]
        builder.opt_workdir = workdir
        builder.opt_timens = opts[KEYS.sim_time]
        builder.opt_queuetime = opts[KEYS.mdr_queue_time]
        builder.file_out = os.path.join(path, job_name + suffix, job_name +
                                        ".out")
        builder.file_mdp = mdp_in + ".mdp"
        builder.file_gro = os.path.join("..", gro_in + "-in.gro")
        builder.file_top = os.path.join("..", base_name + ".top")
        builder.file_ndx = os.path.join("..", base_name + ".ndx")
        builder.file_tpr = job_name + ".tpr"
        builder.extra = extra

        file_ = open(file_name, 'w')
        file_.write(builder.compile())
        file_.close()

        if submit:
            tpr_files.append(os.path.join(folder, job_name + ".tpr"))
            xtc_files.append(os.path.join(folder, job_name + ".xtc"))
            if move_gro:
                gro_files.append(os.path.join(folder, gro_in + "-in.gro"))
            else:
                gro_files.append(os.path.join(gro_in + "-in.gro"))
            edr_files.append(os.path.join(folder, job_name + ".edr"))
            xvg_files.append(os.path.join(folder, job_name + ".xvg"))
            if not opts[KEYS._dryrun]:
                num = job_utils.submit_slurm(file_name, job_name + suffix)
#                 global SAVE_LIBRARY
#                 SAVE_LIBRARY[save_keys.jobs].append((job_name + suffix, num))
            else:
                print("DRYRUN: Would sbatch job " + job_name + suffix)

#         global SAVE_LIBRARY
#         jname = job_name + suffix
#         SAVE_LIBRARY[save_keys.files].append(dir_ + file_name)
#         SAVE_LIBRARY[save_keys.files].append(dir_ + opts[KEYS._params_out])
#         if submit:
#             SAVE_LIBRARY[save_keys.folders].append((dir_ + folder, jname, job_name))

    if submit:
        os.chdir(dir_)
        # Create analysis options file
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

        filename = "analysis_" + opts[KEYS.job_name] + ".bep"
        filepath = os.path.realpath(filename)
        analysis.write(filepath)

#         SAVE_LIBRARY[save_keys.files].append(filepath)

    os.chdir(fir.cwd)
    return fir

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

    fep, coul, vdw, weights = [], [], [], []
    fep.extend(opts[KEYS.sim_fep_values])
    coul.extend(opts[KEYS.sim_coul_values])
    vdw.extend(opts[KEYS.sim_vdw_values])
    weights.extend(opts[KEYS.sim_weight_values])
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
        fep_lambdas += '{0:0.2f} '.format(f)
        coul_lambdas += '{0:0.2f} '.format(c)
        vdw_lambdas += '{0:0.2f} '.format(v)

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

    dt = opts[KEYS.sim_dt]
    steps = opts[KEYS.sim_time] * (1000 / dt)
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
    builder.opt_nstout = opts[KEYS.sim_nstout]
    builder.opt_temp = opts[KEYS.sim_temperature]
    builder.opt_dt = dt;  # ps
    builder.opt_timens = opts[KEYS.sim_time]
    builder.opt_genseed = genseed
    builder.opt_tempalg = opts[KEYS.sim_temp_alg]
    builder.opt_ligand = opts[KEYS.ligand]
    builder.opt_feplambdas = fep_lambdas
    builder.opt_coullambdas = coul_lambdas
    builder.opt_vdwlambdas = vdw_lambdas
    builder.opt_initlambda = init_lambda
    builder.opt_initlambdaindex = state_index
    builder.opt_nstmc = opts[KEYS.sim_nst_mc]
    builder.opt_weightsequil = equil
    builder.opt_lmcseed = lmcseed
    builder.opt_gibbsdelta = int(opts[KEYS.sim_gibbs_delta])
    builder.opt_wlweights = wl_weights
    builder.opt_wlincr = opts[KEYS.sim_incrementor]
    builder.opt_wlscale = opts[KEYS.sim_wl_scale]
    builder.opt_wlratio = opts[KEYS.sim_wl_ratio]

    file_name = os.path.realpath(name + '.mdp')
    file_ = open(file_name, 'w')
    file_.write(builder.compile())
    file_.close()
#     SAVE_LIBRARY[save_keys.files].append(file_name)

    os.chdir(cur_dir)
    return True

def gen_exit(_):
    pass

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

            logscan = job_utils.LogScan(file_name)
            logscan.scan()
            weights = logscan.get_wanglandau_weights()
            opts[KEYS.sim_weight_values] = weights

            os.chdir(work_dir)
            param.write_options('rand#47#.save', {'weights': weights},
                None, None)

            os.chdir(cur_dir)
        except Exception as ex:
            print(ex)
            traceback.print_exc()
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
            file_name = 'rand#47#.save'
            if not os.path.exists(file_name):
                raise Exception('Save file not found ' + file_name)
            save_data = param.parse_options(file_name, dict())
            if isinstance(save_data, list):  # v2.0+
                save_data = save_data[0]
            opts[KEYS.sim_weight_values] = save_data['weights']

            os.chdir(os.path.join(work_dir, opts[KEYS.job_name] + '-rand'))
            file_name = '{0}-rand.log'.format(opts[KEYS.job_name])

            logscan = job_utils.LogScan(file_name)
            logscan.scan()
            num_of_steps = logscan.get_step_number()
            # print("ns:" + str(num_of_steps))
            # num_frames = math.floor(num_of_steps / opts[KEYS.sim_nstout])

            base_name = opts[KEYS.base_name]
            if not base_name:
                base_name = opts[KEYS.job_name]
            xtc_name = '{0}-rand.xtc'.format(opts[KEYS.job_name])
            tpr_name = '{0}-rand.tpr'.format(opts[KEYS.job_name])
            if not os.path.exists(xtc_name):
                raise Exception('Simulation files not found ' + xtc_name +
                    ' and/or ' + tpr_name)
            cmnd = ('echo "System" | trjconv -f {xtc} -s {tpr} -o {gro}' +
                ' -b {timeb} -e {timee} &>../trjconv{index}.log')
            segments = (opts[KEYS.mdr_count] - 1)
            if segments == 0:
                spacing = 1  # will only be multiplied by 0
            else:
                spacing = num_of_steps / segments
                # print("spac:" + str(spacing))

            if ('GMXBIN' in os.environ and os.environ['GMXBIN'] and
                os.environ['GMXBIN'] not in os.environ['PATH']):
                print('Adding GMXBIN to PATH...')
                os.environ['PATH'] = os.environ['PATH'] + ':' + os.environ['GMXBIN']

            for i in range(segments + 1):
                suffix = '-{0:0>2}'.format(i)
                step_num = math.floor(i * spacing)  # evenly spaced frames
                timeps = step_num * opts[KEYS.sim_dt]
                # print("ps:" + str(timeps))
                modif = 0.5 * opts[KEYS.sim_dt] * opts[KEYS.sim_nstout]
                tb = timeps - modif
                te = timeps + modif
                gro_name = '../' + base_name + suffix + '-in.gro'
                fmt = cmnd.format(xtc=xtc_name, tpr=tpr_name, gro=gro_name,
                    timeb=tb, timee=te, index=suffix)
                os.system(fmt)
                if not os.path.exists(gro_name):
                    raise Exception('Extraction of starting frames' +
                        ' unsuccessful - Err: ' + gro_name)

            os.chdir(cur_dir)
        except Exception as ex:
            print(ex)
            traceback.print_exc()
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

def sim_status(save_lib):
    for folder, job, name in save_lib[save_keys.folders]:
        logfile = os.path.join(folder, name + ".log")
        status = "STATUS UNKNOWN"
        extra = None
        warn = None
        printshake = False
        shakecount = 0
        step = 0
        if os.path.exists(logfile):
            scan = job_utils.LogScan(logfile)
            try:
                scan.scan()
                step = scan.get_step_number()
                status = "IN PROGRESS"
                outfile = os.path.join(folder, name + ".out")
                if os.path.exists(outfile):
                    for line in open(outfile, "r"):
                        if 'slurmstepd:' in line:
                            extra = line.replace("slurmstepd:", "")
                        if 'Shake did not converge' in line and printshake:
                            shakecount += 1
                            warn = 'SHAKE experienced issues x{0}'.format(
                                shakecount)
                        if ' fault' in line:
                            status = "FAULT"
                            extra = line
                if scan.finish_detected:
                    status = "FINISHED"
                if scan.cancel_detected:
                    status = "CANCELED"
                if scan.failure_detected:
                    status = "FAILED"
                    extra = scan.fail_statement
            except Exception as ex:
                print(ex)
                traceback.print_exc()
        else:
            status = "NOT STARTED"
        print("{0:<16s} (Step {1:>10}) < ".format(status, step) + job)
        if warn:
            print ("\t" + warn)
        if extra:
            extras = extra.split("\n")
            for ex in extras:
                print("\t" + ex)

def sim_submit(save_lib):
    submitted = False
    for filename in save_lib[save_keys.files]:
        filename = str(filename)
        if os.path.exists(filename) and filename.endswith('.slurm'):
            submitted = True
            job = os.path.basename(filename)
            num = job_utils.submit_slurm(filename, job)
            save_lib[save_keys.jobs].append((job, num))
    if submitted:
        saver.write_options(save_lib[save_keys.name], save_lib)

def sim_cancel(save_lib):
    for name, num in save_lib[save_keys.jobs]:
        if verbose > 0:
            print('Cancelling ' + name + ', job #' + str(num))
        os.system("scancel " + str(num))

def sim_clean(save_lib):
    for filename in save_lib[save_keys.files]:
        if os.path.exists(filename):
            if verbose > 0:
                print('Removing ' + filename)
            os.system("rm " + filename)
    for folder, _, _ in save_lib[save_keys.folders]:
        if os.path.exists(folder):
            if verbose > 0:
                print('Removing folder ' + folder)
            os.system("rm -r " + folder)

def setup(options, args, opts, parser, cur_dir, save_name):
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

    if opts[KEYS.mdr_seedrand] == None:
        import random
        opts[KEYS.mdr_seedrand] = random.randint(0, 65536)

    if subcommand in POST_COMMANDS:
        if 'status' in subcommand and save_name == None:
            here = os.getcwd()
            for dir_item in os.listdir(here):
                item = os.path.join(here, dir_item);
                if os.path.isfile(item) and item.endswith(".save"):
                    save_name = dir_item
                    break
        if save_name == None:
            print(subcommand + " requires a save file to be specified")
            return False  # don't print save files
        else:
            save_path = os.path.realpath(save_name)
            print('Reading previous launch files from ' + save_path)
            save_lib = saver.parse_options(save_path)
            if isinstance(save_lib, list):
                save_lib = save_lib[0]
            save_lib[save_keys.name] = save_path
            SUBS[subcommand](save_lib)
            return False  # don't print save files

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
    return True  # print save files

def main(argv=None):
    if argv == None:
        argv = sys.argv[1:]

    cur_dir = os.getcwd()
    name_ = os.path.basename(__file__).replace('.pyc', '').replace('.py', '')
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
    parser.add_option("--save", help="File generated by a previous run of the" +
        " code, needed to check status or other similar tasks",
        default=None, metavar=("file.save"))
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
        "\t  {0:<12s} - cancel a series of tasks run by {1}\n".format('cancel',
            name_) +
        "\t  {0:<12s} - delete a series of tasks run by {1}\n".format('clean',
            name_) +
        "\t  {0:<12s} - check a series of tasks run by {1}\n".format('status',
            name_) +
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

    save_name = None
    if options.save:
        save_name = options.save

#     global SAVE_LIBRARY
#     SAVE_LIBRARY = job_utils.save_defaults()

    do_output = False
    if isinstance(opt_list, list):
        if options.verbose:
            print("List received, parameters are v2+, {0} entries".format(
                len(opt_list)))
        for opts in opt_list:
            do_output = (setup(options, args, opts, parser, cur_dir, save_name)
                or do_output)
        opts = opt_list[0]
    else:
        if options.verbose:
            print("Dictionary received, parameters are v1.x")
        do_output = setup(options, args, opt_list, parser, cur_dir, save_name)
        opts = opt_list

    if not do_output:
        return 0

    if options.par:
        par_base = os.path.basename(options.par)
        save_name = par_base.replace(param._file_ext(), "") + ".save"
        save_out = backup.backup_file('./', save_name, verbose=verbose)
    else:
        save_name = name_ + ".save"
        save_out = backup.backup_file('', save_name, verbose=verbose)
    saver.write_options(save_out, SAVE_LIBRARY)

POST_COMMANDS.extend(['status', 'submit', 'cancel', 'clean'])
SUBS.update({'exit': gen_exit, 'all': gen_all, 'equil': gen_equil,
    'rand': gen_rand, 'array': gen_array, 'mdp': make_mdp,
    'opt': gen_opt, 'status': sim_status, 'submit': sim_submit,
    'cancel': sim_cancel, 'clean': sim_clean})
KEYS.update()
if __name__ == '__main__':
    # initialize the subcommand list
    main(sys.argv[1:])


