'''
Created on Sep 29, 2015

@author: Tyler P.
'''

class MDP:
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Defines all the variables allowed in the MDP file
        '''
        #################################
        # run control
        #--------------------------------
        self.integrator = None
        self.tinit = None
        self.dt = None
        self.nsteps = None
        self.comm_mode = None
        self.nstcomm = None
        self.comm_grps = None
        self.define = None
        self.include = None
        #################################
        # energy minimization criteria
        #--------------------------------
        self.emtol = None
        self.emstep = None
        self.niter = None
        self.nbfgscorr = None
        #################################
        # output control
        #--------------------------------
        self.nstxout = None
        self.nstvout = None
        self.nstfout = None
        self.nstlog = None
        self.nstcalcenergy = None
        self.nstenergy = None
        self.nstxtcout = None
        self.xtc_precision = None
        self.nstxout_compressed = None
        self.compressed_x_precision = None
        #################################
        # neighbor searching
        #--------------------------------
        self.nstlist = None
        self.ns_type = None
        self.pbc = None
        self.rlist = None
        #################################
        # electrostatics
        #--------------------------------
        self.cutoff_scheme = None
        self.coulombtype = None
        self.coulomb_modifier = None
        self.rcoulomb_switch = None
        self.rcoulomb = None
        #################################
        # van der Waals
        #--------------------------------
        self.vdw_type = None
        self.vdw_modifier = None
        self.rvdw_switch = None
        self.rvdw = None
        #################################
        # corrections
        #--------------------------------
        self.DispCorr = None
        #################################
        # particle mesh ewald
        #--------------------------------
        self.fourierspacing = None
        self.fourier_nx = None
        self.fourier_ny = None
        self.fourier_nz = None
        self.pme_order = None
        self.ewald_rtol = None
        self.ewald_geometry = None
        self.epsilon_surface = None
        self.optimize_fft = None
        #################################
        # coupling
        #--------------------------------
        self.tcoupl = None
        self.nsttcouple = None
        self.tc_grps = None
        self.tau_t = None
        self.ref_t = None
        self.pcoupl = None
        self.pcoupltype = None
        self.nstpcouple = None
        self.tau_p = None
        self.compressibility = None
        self.ref_p = None
        #################################
        # free energy control
        #--------------------------------
        self.free_energy = None
        self.init_lambda = None
        self.delta_lambda = None
        self.foreign_lambda = None
        self.calc_lambda_neighbors = None
        self.sc_alpha = None
        self.sc_power = None
        self.sc_sigma = None
        self.couple_moltype = None
        self.couple_lambda0 = None
        self.couple_lambda1 = None
        self.couple_intramol = None
        self.init_lambda_state = None
        self.fep_lambdas = None
        self.coul_lambdas = None
        self.vdw_lambdas = None
        self.symmetrized_transition_matrix = None
        self.nst_transition_matrix = None
        self.nstdhdl = None
        self.dhdl_print_energy = None
        #################################
        # expanded ensemble
        #--------------------------------
        self.nstexpanded = None
        self.lmc_stats = None
        self.lmc_weights_equil = None
        self.lmc_seed = None
        self.lmc_move = None
        self.lmc_gibbsdelta = None
        self.init_lambda_weights = None
        self.weight_equil_wl_delta = None
        #################################
        # wang landau
        #--------------------------------
        self.wl_scale = None
        self.wl_ratio = None
        self.init_wl_delta = None
        self.wl_oneovert = None
        #################################
        # velocities
        #--------------------------------
        self.gen_vel = None
        self.gen_temp = None
        self.gen_seed = None
        #################################
        # bond constraints
        #--------------------------------
        self.constraints = None
        self.constraint_algorithm = None
        self.shake_tol = None
        self.continuation = None
        self.lincs_order = None
        #################################
        #################################
        # comments
        self.comment_coupling1 = None
        self.comment_coupling2 = None
        self.comment_coupling3 = None
        self.comment_velocities1 = None
        self.comment_constraints1 = None
        self.comment_constraints2 = None

    def reset(self):
        self.__init__()

    def compile(self):
        text = self.run_control()
        text += self.em_criteria()
        if text != "":
            text += "\n"
        text += self.output_control()
        text += self.neighborsearching()
        if text != "":
            text += "\n"
        text += self.electrostatics()
        text += self.vanderWaals()
        text += self.corrections()
        text += self.pme_pppm()
        if text != "":
            text += "\n"
        text += self.coupling()
        if text != "":
            text += "\n"
        text += self.free_energy_control()
        if text != "":
            text += "\n"
        text += self.expanded_ensemble()
        text += self.wang_landau()
        if text != "":
            text += "\n"
        text += self.velocities()
        if text != "":
            text += "\n"
        text += self.bond_constraints()
        return text

    def format_section(self, varlist):
        frmt = "{0:<24s} = {1}\n"
        text = ""
        for pair in varlist:
            key = str(pair[0])
            if pair[1] != None:
                value = str(pair[1])
                if key == "":  # comment
                    text += "; " + value + "\n"
                else:
                    text += frmt.format(key, value)
        return text

    def check_empty(self, varlist):
        '''
        Generic checking, please implement custom checking if needed
        '''
        screen = False
        for pair in varlist:
            if pair[1] != None:
                screen = True
        return not screen

    def run_control(self):
        varlist = [("integrator", self.integrator),
                   ("tinit", self.tinit),
                   ("dt", self.dt),
                   ("nsteps", self.nsteps),
                   ("comm-mode", self.comm_mode),
                   ("nstcomm", self.nstcomm),
                   ("comm-grps", self.comm_grps),
                   ("define", self.define),
                   ("include", self.include)]

        if self.check_empty(varlist):
            return ""

        text = "; Run control\n"
        text += self.format_section(varlist)
        return text

    def em_criteria(self):
        varlist = [("emtol", self.emtol),
                   ("emstep", self.emstep),
                   ("niter", self.niter),
                   ("nbfgscorr", self.nbfgscorr)]

        if self.check_empty(varlist):
            return ""

        text = "; EM criteria and other stuff\n"
        text += self.format_section(varlist)
        return text

    def output_control(self):
        varlist = [("nstxout", self.nstxout),
                   ("nstvout", self.nstvout),
                   ("nstfout", self.nstfout),
                   ("nstlog", self.nstlog),
                   ("nstcalcenergy", self.nstcalcenergy),
                   ("nstenergy", self.nstenergy),
                   ("nstxtcout", self.nstxtcout),
                   ("xtc-precision", self.xtc_precision),
                   ("nstxout-compressed", self.nstxout_compressed),
                   ("compressed-x-precision",
                        self.compressed_x_precision), ]

        if self.check_empty(varlist):
            return ""

        text = "; Output control\n"
        text += self.format_section(varlist)
        return text

    def neighborsearching(self):
        varlist = [("nstlist", self.nstlist),
                   ("ns_type", self.ns_type),
                   ("pbc", self.pbc),
                   ("rlist", self.rlist)]

        if self.check_empty(varlist):
            return ""

        text = ("; Neighborsearching and short-range nonbonded" +
            " interactions\n")
        text += self.format_section(varlist)
        return text

    def electrostatics(self):
        varlist = [("cutoff-scheme", self.cutoff_scheme),
                   ("coulombtype", self.coulombtype),
                   ("coulomb-modifier", self.coulomb_modifier),
                   ("rcoulomb-switch", self.rcoulomb_switch),
                   ("rcoulomb", self.rcoulomb)]

        if self.check_empty(varlist):
            return ""

        text = "; Electrostatics\n"
        text += self.format_section(varlist)
        return text

    def vanderWaals(self):
        varlist = [("vdw-type", self.vdw_type),
                   ("vdw-modifier", self.vdw_modifier),
                   ("rvdw-switch", self.rvdw_switch),
                   ("rvdw", self.rvdw)]

        if self.check_empty(varlist):
            return ""

        text = "; van der Waals\n"
        text += self.format_section(varlist)
        return text

    def corrections(self):
        varlist = [("DispCorr", self.DispCorr)]

        if self.check_empty(varlist):
            return ""

        text = ("; Apply long range dispersion corrections for Energy" +
            " and Pressure \n")
        text += self.format_section(varlist)
        return text

    def pme_pppm(self):
        varlist = [("fourierspacing", self.fourierspacing),
                   ("fourier_nx", self.fourier_nx),
                   ("fourier_ny", self.fourier_ny),
                   ("fourier_nz", self.fourier_nz),
                   ("pme_order", self.pme_order),
                   ("ewald_rtol", self.ewald_rtol),
                   ("ewald_geometry", self.ewald_geometry),
                   ("epsilon_surface", self.epsilon_surface),
                   ("optimize_fft", self.optimize_fft)]

        if self.check_empty(varlist):
            return ""

        text = "; Spacing for the PME/PPPM FFT grid\n"
        varlist.insert(1, ("", "EWALD/PME/PPPM parameters"))
        text += self.format_section(varlist)
        return text

    def coupling(self):
        varlist = [("", self.comment_coupling1),
                   ("", self.comment_coupling2),
                   ("tcoupl", self.tcoupl),
                   ("nsttcouple", self.nsttcouple),
                   ("tc_grps", self.tc_grps),
                   ("tau_t", self.tau_t),
                   ("ref_t", self.ref_t),
                   ("", self.comment_coupling3),
                   ("pcoupl", self.pcoupl),
                   ("pcoupltype", self.pcoupltype),
                   ("nstpcouple", self.nstpcouple),
                   ("tau_p", self.tau_p),
                   ("compressibility", self.compressibility),
                   ("ref_p", self.ref_p)]

        if self.check_empty(varlist):
            return ""

        text = ""
        text += self.format_section(varlist)
        return text

    def free_energy_control(self):
        varlist = [("free_energy", self.free_energy),
                   ("init_lambda", self.init_lambda),
                   ("delta_lambda", self.delta_lambda),
                   ("foreign_lambda", self.foreign_lambda),
                   ("calc-lambda-neighbors", self.calc_lambda_neighbors),
                   ("sc-alpha", self.sc_alpha),
                   ("sc-power", self.sc_power),
                   ("sc-sigma", self.sc_sigma),
                   ("couple-moltype", self.couple_moltype),
                   ("couple-lambda0", self.couple_lambda0),
                   ("couple-lambda1", self.couple_lambda1),
                   ("couple-intramol", self.couple_intramol),
                   ("init-lambda-state", self.init_lambda_state),
                   ("fep-lambdas", self.fep_lambdas),
                   ("coul-lambdas", self.coul_lambdas),
                   ("vdw-lambdas", self.vdw_lambdas),
                   ("symmetrized-transition-matrix",
                        self.symmetrized_transition_matrix),
                   ("nst-transition-matrix", self.nst_transition_matrix),
                   ("nstdhdl", self.nstdhdl),
                   ("dhdl-print-energy", self.dhdl_print_energy)]

        if self.check_empty(varlist):
            return ""

        text = "; Free energy control stuff\n"
        text += self.format_section(varlist)
        return text

    def expanded_ensemble(self):
        varlist = [("nstexpanded", self.nstexpanded),
                   ("lmc-stats", self.lmc_stats),
                   ("lmc-weights-equil", self.lmc_weights_equil),
                   ("lmc-seed", self.lmc_seed),
                   ("lmc-move", self.lmc_move),
                   ("lmc-gibbsdelta", self.lmc_gibbsdelta),
                   ("init-lambda-weights", self.init_lambda_weights),
                   ("weight-equil-wl-delta", self.weight_equil_wl_delta)]

        if self.check_empty(varlist):
            return ""

        text = "; Expanded ensemble\n"
        text += self.format_section(varlist)
        return text

    def wang_landau(self):
        varlist = [("wl-scale", self.wl_scale),
                   ("wl-ratio", self.wl_ratio),
                   ("init-wl-delta", self.init_wl_delta),
                   ("wl-oneovert", self.wl_oneovert)]

        if self.check_empty(varlist):
            return ""

        text = "; Wang-Landau options\n"
        text += self.format_section(varlist)
        return text

    def velocities(self):
        varlist = [("", self.comment_velocities1),
                   ("gen_vel", self.gen_vel),
                   ("gen_temp", self.gen_temp),
                   ("gen_seed", self.gen_seed)]

        if self.check_empty(varlist):
            return ""

        text = ""
        text += self.format_section(varlist)
        return text

    def bond_constraints(self):
        varlist = [("constraints", self.constraints),
                   ("constraint-algorithm", self.constraint_algorithm),
                   ("shake-tol", self.shake_tol),
                   ("", self.comment_constraints1),
                   ("", self.comment_constraints2),
                   ("continuation", self.continuation),
                   ("lincs-order", self.lincs_order)]

        if self.check_empty(varlist):
            return ""

        text = "; options for bonds\n"
        if self.lincs_order != None:
            varlist.insert(6, ("", "Highest order in the expansion" +
            " of the constraint coupling matrix"))
        if self.constraint_algorithm != None:
            varlist.insert(1, ("", "Type of constraint algorithm"))
        text += self.format_section(varlist)
        return text


if __name__ == "__main__":
    # test
    mdp = MDP()
    mdp.ns_type = "grid"
    mdp.pbc = "xyz"
    mdp.rlist = "1.0"
    mdp.coulombtype = "PME"
    mdp.rcoulomb = "1.0"
    mdp.vdw_type = "switch"
    mdp.rvdw_switch = "0.8"
    mdp.rvdw = "0.9"
    mdp.DispCorr = "EnerPres"
    mdp.fourierspacing = "0.12"
    mdp.pme_order = "6"
    mdp.ewald_rtol = "1e-06"
    mdp.epsilon_surface = "0"
    mdp.optimize_fft = "no"
    mdp.free_energy = "yes"
    mdp.init_lambda = "0.0"
    mdp.delta_lambda = "0"
    mdp.foreign_lambda = "0.05"
    mdp.sc_alpha = "0.5"
    mdp.sc_power = "1.0"
    mdp.sc_sigma = "0.3 "
    mdp.couple_moltype = "Methane  ; name of moleculetype to decouple"
    mdp.couple_lambda0 = "vdw      ; only van der Waals interactions"
    mdp.couple_lambda1 = ("none     ; turn off everything, in this" +
        " case only vdW")
    mdp.couple_intramol = "no"
    mdp.nstdhdl = "10"
    mdp.constraint_algorithm = "lincs"
    mdp.lincs_order = "12"
    mdp.integrator = "steep "
    mdp.nsteps = "5000"
    mdp.emtol = "100"
    mdp.emstep = "0.01"
    mdp.niter = "20"
    mdp.nbfgscorr = "10"
    mdp.nstlog = "1"
    mdp.nstenergy = "1"
    mdp.nstlist = "1"
    mdp.comment_coupling1 = ("Temperature and pressure coupling are" +
            " off during EM")
    mdp.tcoupl = "no"
    mdp.pcoupl = "no"
    mdp.comment_velocities1 = "Generate velocities to start"
    mdp.gen_vel = "no "
    mdp.constraints = "h-bonds  ; we only have C-H bonds here"
    mdp.comment_constraints1 = ("Do not constrain the starting" +
            " configuration")
    mdp.continuation = "no"
    print(mdp.compile())


























