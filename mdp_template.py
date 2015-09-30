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
        self.define = None
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
        self.nstenergy = None
        self.nstxtcout = None
        self.xtc_precision = None
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
        self.coulombtype = None
        self.rcoulomb = None
        #################################
        # van der Waals
        #--------------------------------
        self.vdw_type = None
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
        self.pme_order = None
        self.ewald_rtol = None
        self.epsilon_surface = None
        self.optimize_fft = None
        #################################
        # coupling
        #--------------------------------
        self.tcoupl = None
        self.tc_grps = None
        self.tau_t = None
        self.ref_t = None
        self.pcoupl = None
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
        self.sc_alpha = None
        self.sc_power = None
        self.sc_sigma = None
        self.couple_moltype = None
        self.couple_lambda0 = None
        self.couple_lambda1 = None
        self.couple_intramol = None
        self.nstdhdl = None
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
            value = str(pair[1])
            if value != None:
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
                   ("define", self.define)]

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
                   ("nstenergy", self.nstenergy),
                   ("nstxtcout", self.nstxtcout),
                   ("xtc-precision", self.xtc_precision)]

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
        varlist = [("coulombtype", self.coulombtype),
                   ("rcoulomb", self.rcoulomb)]

        if self.check_empty(varlist):
            return ""

        text = "; Electrostatics\n"
        text += self.format_section(varlist)
        return text

    def vanderWaals(self):
        varlist = [("vdw-type", self.vdw_type),
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
                   ("pme_order", self.pme_order),
                   ("ewald_rtol", self.ewald_rtol),
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
                   ("tc_grps", self.tc_grps),
                   ("tau_t", self.tau_t),
                   ("ref_t", self.ref_t),
                   ("", self.comment_coupling3),
                   ("pcoupl", self.pcoupl),
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
                   ("sc-alpha", self.sc_alpha),
                   ("sc-power", self.sc_power),
                   ("sc-sigma", self.sc_sigma),
                   ("couple-moltype", self.couple_moltype),
                   ("couple-lambda0", self.couple_lambda0),
                   ("couple-lambda1", self.couple_lambda1),
                   ("couple-intramol", self.couple_intramol),
                   ("nstdhdl", self.nstdhdl)]

        if self.check_empty(varlist):
            return ""

        text = "; Free energy control stuff\n"
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
                   ("", self.comment_constraints1),
                   ("", self.comment_constraints2),
                   ("continuation", self.continuation),
                   ("lincs-order", self.lincs_order)]

        if self.check_empty(varlist):
            return ""

        text = "; options for bonds\n"
        varlist.insert(6, ("", "Highest order in the expansion of the" +
            " constraint coupling matrix"))
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


























