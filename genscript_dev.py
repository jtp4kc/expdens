#!/usr/bin/env python
'''
Created on Jun 19, 2015

@author: Tyler P
'''
import sys, os
import traceback
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
import math
import random
import backup
from param_versions.version_2_0 import Parameters
from param_versions.version_2_0 import Keys
import mdp_template
import slurm_template
from abc import abstractmethod, ABCMeta
import job_utils
import job_daemon

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
        self.work_dir = 'work-dir'
        self.script_dir = 'sim-dir'
        self.input_dir = 'input-file-dir'
        self.ligand = 'ligand-res-name'

        section = "General"
        self.add_key(self.base_name, section, "Input file prefix")
        self.add_key(self.job_name, section, "Name to prefix all files")
        self.add_key(self.work_dir, section, "Directory to operate within," +
            " and create save files. Relative to current directory.")
        self.add_key(self.script_dir, section, "Directory to create script" +
            " folders within, where the jobs will run (relative to " +
            self.work_dir + ")")
        self.add_key(self.input_dir, section, "Location of input files" +
            " (.top, .ndx, -in.gro), relative to " + self.work_dir)
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
        self.mdr_randsrc = 'mdr-randomized-xtc-for-initial-frames'
        self._expected_frames = "mdr-expected-number-of-frames"

        section = "mdrun Array"
        self.add_key(self.mdr_threads, section, 'Max 20 threads (for now)')
        self.add_key(self.mdr_randsrc, section, "Path to xtc to use for" +
            " initial configurations. There must also be a .tpr by the same" +
            " name (relative to here, to parfile, or to workdir, in order)")
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
        self._submit = '_submit-jobs'  # private
        self._dryrun = '_dry-run'  # private

        section = "Process Control"
        self.add_key(self.subcommand, section=section, updater=_subcomment)
        self.add_keys(section, self.verbosity, self.run_mdp, self.run_array)

KEYS = MyKeys()
param = Parameters(KEYS)
ATTR = job_daemon.Attr()

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
    options[KEYS.input_dir] = '.'
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
    options[KEYS.mdr_genseed] = True
    options[KEYS.mdr_seedrand] = None
    options[KEYS.mdr_randsrc] = None
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

class FilesystemImpactRegister:

    class Entry():
        __metaclass__ = ABCMeta

        def __init__(self, filename):
            if filename != None:
                fname = str(filename)
            else:
                raise Exception("Filename cannot be 'None'")
            self.filename = fname

        @abstractmethod
        def make(self, filename):
            """ Implement custom code to write using filename
            """
            pass

    class StreamEntry(Entry):
        __metaclass__ = ABCMeta

        def __init__(self, filename):
            FilesystemImpactRegister.Entry.__init__(self, filename)

        def make(self, filename):
            output = open(filename, "w")
            self.write(output)
            output.close()

        @abstractmethod
        def write(self, stream):
            """ Implement custom code to write to stream
            """
            pass

    def __init__(self, directory=None):
        self.dirs_to_create = []
        self.dirs_to_remove = []
        self.files_to_create = []
        self.files_to_remove = []
        self.directory = directory

    def add_create(self, entry_to_create):
        good = False
        try:
            if entry_to_create.make and entry_to_create.filename:
                good = True
        except:
            print("Please use the Entry inner class with this method")

        if good:
            self.files_to_create.append(entry_to_create)

    def add_remove(self, filename_to_remove):
        val = self._add_check(filename_to_remove)
        if val != None:
            self.files_to_remove.append(val)

    def add_mkdir(self, path_to_make):
        val = self._add_check(path_to_make)
        if val != None:
            self.dirs_to_create.append(val)

    def add_rmdir(self, path_to_del):
        val = self._add_check(path_to_del)
        if val != None:
            self.dirs_to_remove.append(val)

    def _add_check(self, name):
        if name != None:
            return str(name)
        else:
            print("Path cannot be 'None'")
            return None

    def merge(self, other_FIR):
        fir = FilesystemImpactRegister()
        fir.directory = self.directory
        fir.files_to_create.extend(self.files_to_create)
        fir.files_to_create.extend(other_FIR.files_to_create)
        fir.files_to_remove.extend(self.files_to_remove)
        fir.files_to_remove.extend(other_FIR.files_to_remove)
        return fir

    def describe(self, prefix="", print_method=None):
        if print_method == None:
            def alias(*args):
                for arg in args:
                    print arg,
                print
            print_method = alias

        dir_ = self.directory
        if dir_ is None:
            dir_ = "cwd"
        print_method(prefix + "From " + dir_ + ", there are instructions to" +
                     " ...")

        for path in self.dirs_to_remove:
            print_method(prefix + "remove directory " + path)

        for path in self.dirs_to_create:
            print_method(prefix + "create directory " + path)

        for entry in self.files_to_create:
            print_method(prefix + "create file " + entry.filename)

        for fname in self.files_to_remove:
            print_method(prefix + "remove file " + fname)

    def execute(self, verbose=True, quiet=False):
        """ Perform all the file creation and removal assigned to this register
        """
        v = False
        q = False
        if verbose:  # whatever the current rule for boolean evaluation is
            v = True
        if quiet:
            q = True

        here = os.getcwd()
        dir_ = self.directory
        if dir_ is None:
            dir_ = here
        if not os.path.isdir(dir_):
            if v:
                print("Creating directory " + dir_)
            os.mkdir(dir_)
        if v:
            print("Changing to directory " + dir_)
        os.chdir(dir_)

        for path in self.dirs_to_remove:
            if os.path.exists(path) and os.path.isdir(path):
                files = os.listdir(path)
                if v:
                    print("Removing directory and all its contents: " + path)
                for fname in files:
                    print("\t" + fname)
                os.system("rm -r " + path)
            elif os.path.exists(path):
                if not q:
                    print("Warning: is not a directory " + path)
            elif not q:
                print("Warning: directory does not exist " + path)

        for path in self.dirs_to_create:
            if os.path.exists(path) and os.path.isdir(path):
                if not q:
                    print("Warning: directory already exists " + path)
                continue
            if v:
                print("Creating directory " + path)
            os.mkdir(path)

        to_remove = list(self.files_to_remove)
        for entry in self.files_to_create:
            if entry.filename in to_remove:
                to_remove.remove(entry.filename)
            elif os.path.exists(entry.filename) and not q:
                print("Warning: overwriting file " + entry.filename)
            if v:
                print("Creating file " + entry.filename)
            entry.make(entry.filename)

        for fname in to_remove:
            if os.path.exists(fname):
                if v:
                    print("Removing file " + fname)
                os.remove(fname)
            elif not q:
                print("Cannot remove " + fname + ", file does not exist.")

        if v:
            print("Returning to directory " + here)
        os.chdir(here)

class SlurmGen(slurm_template.Slurm):

    def __init__(self):
        slurm_template.Slurm.__init__(self)
        self.double_precision = False
        self.use_mpi = False
        self.gromacs5 = False
        self.calc_wt = True
        self.extra = None

        self.opt_jobname = ""
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
        rate = 4.1 / 20  # 4.1 ns/day per 20 cores for double precicion
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

        self.job_name = self.opt_jobname
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
            self.modules.append("gromacs/5.0.2-sse")
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
        self.cmt_weights1 = None
        self.cmt_weights2 = None

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
            self.comment_weights1 = self.cmt_weights1
            self.comment_weights2 = self.cmt_weights2
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

class WriterEntry(FilesystemImpactRegister.Entry):

    def __init__(self, filename, writer_item):
        FilesystemImpactRegister.Entry.__init__(self, filename)
        self.item = writer_item

    def make(self, filename):
        self.item.write(filename)

class CompiledEntry(FilesystemImpactRegister.StreamEntry):

    def __init__(self, filename, compiled_item):
        FilesystemImpactRegister.Entry.__init__(self, filename)
        self.item = compiled_item

    def write(self, stream):
        text = self.item.compile()
        stream.write(text)

def handle_job(opts):
    job = opts[KEYS.job_name]
    wrkdir = os.path.join(backup.expandrelpath(opts[KEYS.work_dir]), "")
    savename = os.path.join(wrkdir, job + ".save")
    jobsave = job_utils.SaveJobs()
    jobsave.attr["name"] = job

    make_job(opts, jobsave)
    if not opts[KEYS._dryrun]:
        jobsave.save(savename)
    if opts[KEYS._submit]:
        job_daemon.reschedule_self(job, savename, time="7-00:00:00", live=True)

def make_job(opts, jobsave=None):
    job_name = opts[KEYS.job_name]
    if not job_name:
        print('Job name not specified')
        return
    base_name = opts[KEYS.base_name]
    if not base_name:
        base_name = job_name

    wrkdir = os.path.join(backup.expandrelpath(opts[KEYS.work_dir]), "")
    fir = FilesystemImpactRegister(directory=wrkdir)
    path = os.path.join(backup.expandrelpath(opts[KEYS.script_dir]), "")
    indir = os.path.join(backup.expandrelpath(opts[KEYS.input_dir]), "")

    if not os.path.isdir(path):
        fir.add_mkdir(path)
    bep_dir = path  # in case it decides to move

    tpr_files = []
    xtc_files = []
    gro_files = []
    xvg_files = []
    to_submit = []

    if opts[KEYS.mdr_seedrand] != None:
        random.seed(opts[KEYS.mdr_seedrand])

    n_sim = opts[KEYS.mdr_count]

    gro_gen = None
    frames = None
    if opts[KEYS.mdr_randsrc] is not None:
        randname = opts[KEYS.mdr_randsrc]
        randxtc = backup.expandrelpath(randname)
        if not os.path.exists(randxtc):
            randxtc = backup.expandrelpath(randname,
                                           os.path.dirname(opts[KEYS._params]))
        if not os.path.exists(randxtc):
            randxtc = backup.expandrelpath(randname, wrkdir)
        if verbose > 0:
            print("Trajectory to use to make random input frames: " + randxtc)
        gro_gen = job_utils.ExtractFrames(randxtc)
        try:
            frames = gro_gen.get_gro_frames(n_sim, framenums=True)
        except Exception as e:
            print(e)

    for i in range(n_sim):
        suffix = '-{0:0>2}'.format(i)
        if opts[KEYS.mdr_count] < 2:
            suffix = ''

        dir_name = job_name + suffix
        jobdir = os.path.join(path, dir_name, "")
        if not os.path.exists(os.path.join(wrkdir, jobdir)):
            fir.add_mkdir(jobdir)

        # filenames
        slurmname = os.path.join(jobdir, job_name + suffix + '.slurm')
        mdpname = os.path.join(jobdir, job_name + suffix + '.mdp')
        tprname = os.path.join(jobdir, job_name + suffix + '.tpr')
        xtcname = os.path.join(jobdir, job_name + suffix + '.xtc')
        xvgname = os.path.join(jobdir, job_name + suffix + '.xvg')
        logname = os.path.join(jobdir, job_name + suffix + '.log')
        outname = os.path.join(jobdir, job_name + suffix + '.out')
        ndxname = os.path.join(indir, base_name + '.ndx')
        topname = os.path.join(indir, base_name + '.top')
        groname = os.path.join(indir, base_name + '-in.gro')
        if frames is not None:
            groname = os.path.join(jobdir, job_name + suffix + '-in.gro')

        if opts[KEYS.mdr_genseed]:
            num = random.randint(0, 65536)
            buildmdp = make_mdp(opts, genseed=num, lmcseed=num)
        else:
            buildmdp = make_mdp(opts)

        if buildmdp is None:
            return
        fir.add_create(CompiledEntry(mdpname, buildmdp))

        eqtime = opts[KEYS.sim_time] / 2
        temp = opts[KEYS.sim_temperature]
        extra = ("python ~/git/alchemical-analysis/alchemical_analysis/" +
            "alchemical_analysis.py -p dhdl -r 8 -s " + str(eqtime) + " -t " +
            str(temp) + " -u kBT -v -w -x &> alchem-" + job_name + ".out\n")

        loc = os.path.dirname(slurmname)
        builder = SlurmGen()
        builder.double_precision = opts[KEYS.sim_precision]
        builder.use_mpi = opts[KEYS.sim_use_mpi]
        builder.calc_wt = opts[KEYS.auto_calc_wt]
        builder.gromacs5 = opts[KEYS.sim_gromacs5]

        slurmworkdir = backup.expandrelpath(jobdir, loc)
        if slurmworkdir == ".":
            slurmworkdir = None
        builder.opt_jobname = job_name + suffix
        builder.opt_ntasks = opts[KEYS.mdr_threads]
        builder.opt_workdir = slurmworkdir
        builder.opt_timens = opts[KEYS.sim_time]
        builder.opt_queuetime = opts[KEYS.mdr_queue_time]
        builder.file_out = job_name + suffix + ".out"

        builder.file_mdp = backup.expandrelpath(mdpname, loc)
        builder.file_gro = backup.expandrelpath(groname, loc)
        builder.file_top = backup.expandrelpath(topname, loc)
        builder.file_ndx = backup.expandrelpath(ndxname, loc)
        builder.file_tpr = backup.expandrelpath(tprname, loc)
        builder.extra = extra

        fir.add_create(CompiledEntry(slurmname, builder))

        tpr_files.append(backup.expandrelpath(tprname, bep_dir))
        xtc_files.append(backup.expandrelpath(xtcname, bep_dir))
        gro_files.append(backup.expandrelpath(groname, bep_dir))
        xvg_files.append(backup.expandrelpath(xvgname, bep_dir))

        frame_ps = opts[KEYS.sim_dt] * opts[KEYS.sim_nstout]  # ps b/w frames
        jobentry = job_utils.SaveEntry()
        jobentry.jobname = job_name + suffix
        jobentry.attr[ATTR.RESNAME] = opts[KEYS.ligand]
        jobentry.attr[ATTR.LOG_NSTEP] = opts[KEYS.sim_nstout]
        jobentry.attr[ATTR.LOG_DELTA] = frame_ps
        jobentry.attr[ATTR.JOB_ID] = ""
        jobentry.files["slurm"] = slurmname
        jobentry.files["mdp"] = mdpname
        jobentry.files["tpr"] = tprname
        jobentry.files["xtc"] = xtcname
        jobentry.files["xvg"] = xvgname
        jobentry.files["log"] = logname
        jobentry.files["out"] = outname
        jobentry.files["top"] = topname
        jobentry.files["ndx"] = ndxname
        jobentry.files["gro"] = groname
        jobentry.files["folder"] = jobdir
        to_submit.append(jobentry)
        if jobsave is not None:
            jobsave.add_job(jobentry)

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
    analysis.fields_output[BEPGen.KEYS.dhdl_files] = xvg_files
    analysis.fields_output[BEPGen.KEYS.check_average_energies] = False
    if opts[KEYS.sim_fixed_lambda]:
        analysis.fields_output[BEPGen.KEYS.single_state] = state_index

    bepname = os.path.join(bep_dir, "analyze-" + opts[KEYS.job_name] + ".bep")
    fir.add_create(WriterEntry(bepname, analysis))

    if opts[KEYS._dryrun]:
        fir.describe()
        for entry in to_submit:
            print("DRYRUN: Would sbatch job " + entry.jobname)
    else:
        fir.execute()
        init_gro = [None] * n_sim
        if frames is not None:
            init_gro = gro_gen.get_gro_frames(n_sim)
            if len(init_gro) != n_sim:
                ratio = len(init_gro) / float(n_sim)
                copy_gro = list(init_gro)
                init_gro = [None] * n_sim
                for i in range(n_sim):
                    init_gro[i] = copy_gro[int(ratio * i)]
        for (entry, gro) in zip(to_submit, init_gro):
            if gro is not None:
                os.system("cp " + gro + " " + os.path.join(wrkdir,
                                                           entry.files["gro"]))

            if opts[KEYS._submit]:
                slurmfile = os.path.join(wrkdir, entry.files["slurm"])
                jid = job_utils.submit_job(slurmfile, entry.jobname)
                entry.attr[ATTR.JOB_ID] = jid
        for gro in init_gro:
            if (gro is not None) and os.path.exists(gro):
                os.remove(gro)

def make_mdp(opts, dir_='.', genseed=10200, lmcseed=10200):
    cur_dir = os.getcwd()
    dir_ = os.path.expandvars(os.path.expanduser(dir_))
    if not os.path.exists(dir_):
        os.mkdir(dir_)
    os.chdir(dir_)

    fep = list(opts[KEYS.sim_fep_values])
    coul = list(opts[KEYS.sim_coul_values])
    vdw = list(opts[KEYS.sim_vdw_values])
    weights = list(opts[KEYS.sim_weight_values])
    if verbose > 2:
        print('Debug, genstates')

    orig_weight = list(weights)
    add_weight = [1.0] * len(weights)

    wgtxcoup = int(opts[KEYS.sim_genxcoupled])
    wgtxuncp = int(opts[KEYS.sim_genxuncupld])
    # coupled
    if verbose > 2:
        print('Value coup {0}'.format(wgtxcoup))
    if type(wgtxcoup) in (int, float):
        if verbose > 2:
            print('weighting coupled state')
        add_weight[0] = add_weight[0] * wgtxcoup
    if type(wgtxcoup) in (list, tuple):
        if verbose > 2:
            print('weighting coupled state')
        for i in range(len(wgtxcoup)):
            if (i < len(add_weight)) and (wgtxcoup[i] > 0):
                add_weight[i] = add_weight[i] * wgtxcoup[i]
    # uncoupled
    if verbose > 2:
        print('Value uncoup {0}'.format(opts[KEYS.sim_wgtxuncupld]))
    if type(wgtxuncp) in (int, float):
        if verbose > 2:
            print('weighting uncoupled state')
        if len(weights) > 0:
            add_weight[-1] = add_weight[-1] * wgtxuncp
    if type(wgtxuncp) in (list, tuple):
        if verbose > 2:
            print('weighting coupled state')
        for i in range(len(wgtxuncp)):
            if (i <= len(add_weight)) and (wgtxuncp[-i] > 0):
                add_weight[-i] = add_weight[-i] * wgtxuncp[-i]
    # apply modifcations
    for i in range(len(weights)):
        wgt = add_weight[i] / add_weight[0]
        weights[i] += math.log(wgt)

    genxcoup = int(opts[KEYS.sim_genxcoupled])
    genxuncp = int(opts[KEYS.sim_genxuncupled])
    if genxcoup > 0:
        fep = [0.0] * genxcoup + fep
        vdw = [0.0] * genxcoup + vdw
        coul = [0.0] * genxcoup + coul
        weights = [weights[0]] * genxcoup + weights
    if genxuncp > 0:
        fep = fep + [0.0] * genxuncp
        vdw = vdw + [1.0] * genxuncp
        coul = coul + [1.0] * genxuncp
        weights = weights + [weights[-1]] * genxuncp
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

    wl_weights = ""
    orig = None
    addw = None
    if opts[KEYS.sim_weights]:
        orig = "original weights: "
        addw = "additional ln(factor): "
        for w in weights:
            wl_weights += '{0:0.5f} '.format(w)
        for w in orig_weight:
            orig += '{0:0.5f} '.format(w)
        for w in add_weight:
            addw += '{0:0.5f} '.format(w)
        if not len(fep) == len(weights):
            print('Number of weights not equal to number of states, please' +
                ' address')
            return None

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
    builder.cmt_weights1 = orig
    builder.cmt_weights2 = addw
    builder.opt_wlweights = wl_weights
    builder.opt_wlincr = opts[KEYS.sim_incrementor]
    builder.opt_wlscale = opts[KEYS.sim_wl_scale]
    builder.opt_wlratio = opts[KEYS.sim_wl_ratio]

    os.chdir(cur_dir)
    return builder

def gen_exit(*args):
    pass

def gen_rand(opts):
    opts[KEYS.mdr_count] = 1
    opts[KEYS.sim_fixed_weights] = False
    opts[KEYS.sim_weights] = False
    opts[KEYS.sim_fixed_lambda] = True
    opts[KEYS.sim_time] = opts[KEYS.auto_time_rand]
    if opts[KEYS._job_name]:
        opts[KEYS.job_name] = opts[KEYS._job_name] + '-rand'

    handle_job(opts)

def gen_equil(opts):
    if opts[KEYS.mdr_randsrc] is None:
        randname = opts[KEYS.job_name] + "-rand"
        opts[KEYS.mdr_randsrc] = os.path.join(randname, randname + ".xtc")
    opts[KEYS.mdr_count] = 1
    opts[KEYS.sim_fixed_weights] = False
    opts[KEYS.sim_weights] = False
    opts[KEYS.sim_fixed_lambda] = False
    opts[KEYS.sim_time] = opts[KEYS.auto_time_equil]
    if opts[KEYS._job_name]:
        opts[KEYS.job_name] = opts[KEYS._job_name] + '-equil'

    handle_job(opts)

def gen_array(opts):
    if opts[KEYS.mdr_randsrc] is None:
        randname = opts[KEYS.job_name] + "-rand"
        opts[KEYS.mdr_randsrc] = os.path.join(randname, randname + ".xtc")
    opts[KEYS.sim_fixed_weights] = True
    opts[KEYS.sim_weights] = True
    opts[KEYS.sim_fixed_lambda] = False
    opts[KEYS.sim_time] = opts[KEYS.auto_time_array]
    if opts[KEYS._job_name]:
        opts[KEYS.job_name] = opts[KEYS._job_name]

    handle_job(opts)

def gen_opt(opts):
    if opts[KEYS.run_mdp]:
        builder = make_mdp(opts)
        if builder is not None:
            job_name = opts[KEYS.job_name]
            if not job_name:
                print('Job name not specified')
                text = builder.compile()
                for line in text.splitlines(keepends=True):
                    print("\t" + line)
            else:
                name = job_name + ".mdp"
                file_name = os.path.realpath(name)
                file_ = open(file_name, 'w')
                file_.write(builder.compile())
                file_.close()

    if opts[KEYS.run_array]:
        handle_job(opts)

def sim_status(jobsave):
    for entry in jobsave.jobs:
        logfile = entry.files["log"]
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
                outfile = entry.files["out"]
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
        print("{0:<16s} (Step {1:>10}) < ".format(status, step) + entry.jobname)
        if warn:
            print ("\t" + warn)
        if extra:
            extras = extra.split("\n")
            for ex in extras:
                print("\t" + ex)

def sim_submit(jobsave):
    for entry in jobsave.jobs:
        slurmname = entry.files["slurm"]
        if os.path.exists(slurmname) and slurmname.endswith('.slurm'):
            job = os.path.basename(slurmname)
            num = job_utils.submit_job(slurmname, job)
            entry.attr[ATTR.JOB_ID] = num

def sim_cancel(jobsave):
    for entry in jobsave.jobs:
        jobname = entry.jobname
        jobid = entry.attr[ATTR.JOB_ID]
        if jobid.isdigit():
            if verbose > 0:
                print('Cancelling ' + jobname + ', job #' + jobid)
            os.system("scancel " + jobid)
        elif verbose > 0:
                print(jobname + " appears to not have been submitted.")

def sim_clean(jobsave):
    for entry in jobsave.jobs:
        filelist = []
        filelist.append(entry.files["slurm"])
        filelist.append(entry.files["mdp"])
        filelist.append(entry.files["tpr"])
        filelist.append(entry.files["xtc"])
        filelist.append(entry.files["xvg"])
        filelist.append(entry.files["log"])
        filelist.append(entry.files["out"])
        for filename in filelist:
            if os.path.exists(filename):
                if verbose > 0:
                    print('Removing ' + filename)
                os.remove(filename)
        folder = entry.files["folder"]
        if os.path.isdir(folder):
            if verbose > 0:
                print('Removing folder ' + folder)
            os.system("rm -r " + folder)

def setup(args, opts, parser, cur_dir, save_name):
    global verbose  # ensure that we are talking about the same verbose here
    if args.verbose:
        verbose = args.verbose
    else:
        verbose = opts[KEYS.verbosity]
    opts[KEYS.verbosity] = verbose

    if args.subcommand:
        subcommand = args.subcommand[0]
        if subcommand not in SUBS:
            subcommand = 'exit'
    else:
        subcommand = opts[KEYS.subcommand]

    if subcommand == 'exit':
        parser.print_help()

    if args.n:
        opts[KEYS.job_name] = args.n
    if args.b:
        opts[KEYS.base_name] = args.b
    if args.mdp:
        opts[KEYS.run_mdp] = args.mdp
    if args.slurm:
        opts[KEYS.run_array] = args.slurm
    opts[KEYS.subcommand] = subcommand

    if opts[KEYS.mdr_seedrand] == None:
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
        else:
            save_path = os.path.realpath(save_name)
            print('Reading previous launch files from ' + save_path)
            jobsave = job_utils.SaveJobs()
            jobsave.load(save_path)

            here = os.getcwd()
            os.chdir(os.path.dirname(jobsave.filename))
            SUBS[subcommand](jobsave)
            os.chdir(here)
            jobsave.save()
        return

    # output run options
    if opts[KEYS.job_name]:
        param_name = opts[KEYS.job_name] + param.file_ext
    else:
        param_name = param.file_ext
    wrkdir = os.path.join(backup.expandrelpath(opts[KEYS.work_dir]), "")
    param_name = os.path.join(wrkdir, param_name)

    opts[KEYS._calling_dir] = os.path.join(cur_dir, '')
    opts[KEYS._params] = args.par
    opts[KEYS._params_out] = param_name

    opts[KEYS._job_name] = opts[KEYS.job_name]  # backup the original name
    opts[KEYS._submit] = args.submit
    opts[KEYS._dryrun] = args.dryrun

    if args.dryrun:
        print("DRYRUN: Would write full option file " + param_name)
    else:
        if not os.path.isdir(wrkdir):
            os.mkdir(wrkdir)
        param_out = backup.backup_file('', param_name, verbose=verbose)
        param.write_options(param_out, opts)

    # perform requested file output and job submissions
    SUBS[subcommand](opts)
    return True  # print save files

def main(argv=None):
    if argv == None:
        argv = sys.argv[1:]

    cur_dir = os.getcwd()
    name_ = os.path.basename(__file__).replace('.pyc', '').replace('.py', '')
    print('Invocation of %s:\n\t%s ' % (__file__, name_) + " ".join(argv) +
          '\n')

    opts = option_defaults()
    subhelp = ("indicate a procedure to follow for this call\n" +
        "The following subcommands may be used:\n" +
        "  {0:<12s} - run the Wang-Landau equilibriation\n".format('equil') +
        "  {0:<12s} - run a position randomizing simulation\n".format('rand') +
        "  {0:<12s} - run a batch of simulations\n".format('array') +
        "  {0:<12s} - output a formatted mdp file\n".format('mdp') +
        "  {0:<12s} - perform tasks based on a param file\n".format('opt') +
        "  {0:<12s} - cancel a series of tasks run by {1}\n".format('cancel',
            name_) +
        "  {0:<12s} - delete a series of tasks run by {1}\n".format('clean',
            name_) +
        "  {0:<12s} - check a series of tasks run by {1}\n".format('status',
            name_) +
        "  {0:<12s} - submit a series of tasks run by {1}\n".format('submit',
            name_) +
        "  {0:<12s} - print usage and exit".format('exit'))

    parser = ArgumentParser(description="Generates Rivanna SLURM scripts for" +
        " running Gromacs simulations", prog=name_,
        formatter_class=RawTextHelpFormatter)
    parser.add_argument("subcommand", help=subhelp, nargs='*')
    parser.add_argument("-v", "--verbose", help="Increase output frequency" +
        " and detail. Stacks three times.", action="count")
    parser.add_argument("--par", help="Optional configuration file to" +
        " specify command line parameters and more.",
        default=None, metavar="file.par")
    parser.add_argument("--save", help="File generated by a previous run of" +
        " the code, needed to check status or other similar tasks",
        default=None, metavar=("file.save"))
    parser.add_argument("--slurm", help="""Prepare scripts to run on Rivanna""",
        default=None, action='store_true')
    parser.add_argument("--mdp", help="Output the Molecular Dynamics" +
        " Parameters (.mdp) configuration file",
        default=None, action='store_true')
    parser.add_argument("-n", help="""Job name for output""", metavar="sim")
    parser.add_argument("-b", help="""Base name for job input""", metavar="sim")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--submit", action='store_true', help="Give permission" +
        " to submit the created jobs to the slurm controller using sbatch")
    group.add_argument("--dryrun", action='store_true', help="Overrides" +
        " --submit, make no files but print what would be done")

    if len(argv) < 1:
        parser.print_help()
        return 0

    args = parser.parse_args(argv)

    # handle options, reading from file if requested
    if args.par:
        print('Reading parameters from ' + args.par)
        opt_list = param.parse_options(args.par, opts)
        if not opt_list:
            print('Error reading from parameters file.')
            return 1
    else:
        opt_list = opts

    save_name = None
    if args.save:
        save_name = args.save

    if args.dryrun:
        args.submit = False

    if isinstance(opt_list, list):
        if args.verbose:
            print("List received, parameters are v2+, {0} entries".format(
                len(opt_list)))
        for opts in opt_list:
            setup(args, opts, parser, cur_dir, save_name)
        opts = opt_list[0]
    else:
        if args.verbose:
            print("Dictionary received, parameters are v1.x")
        setup(args, opt_list, parser, cur_dir, save_name)
        opts = opt_list

POST_COMMANDS.extend(['status', 'submit', 'cancel', 'clean'])
SUBS.update({'exit': gen_exit, 'equil': gen_equil,
    'rand': gen_rand, 'array': gen_array, 'mdp': make_mdp,
    'opt': gen_opt, 'status': sim_status, 'submit': sim_submit,
    'cancel': sim_cancel, 'clean': sim_clean})
KEYS.update()
if __name__ == '__main__':
    # initialize the subcommand list
    main(sys.argv[1:])



