'''
Created on Oct 29, 2015

@author: Tyler P.
'''

class Slurm:
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Defines all the known variables in the Slurm file
        '''
        #################################
        # sbatch directives
        #--------------------------------
        self.job_name = None
        self.partition = None
        self.nodes = None
        self.ntasks = None
        self.time = None
        self.signal = None
        self.output = None
        self.workdir = None
        self.checkpoint = None
        self.checkpoint_dir = None
        self.mail_types = []
        self.mail_user = None
        self.comments = []
        #################################
        # slurm import variables
        #--------------------------------
        self.modulepath_prepends = []
        self.modulepath_appends = []
        self.modules = []
        #################################
        # export environment variables
        #--------------------------------
        self.exports = []
        #################################
        # work load - commands to execute
        #--------------------------------
        self.commands = []
        #################################
        # broiler-plate options
        #--------------------------------
        self.status_filename = None
        self.include_diagnostics = True
        self.include_timestart = True
        self.include_timeend = True

    def reset(self):
        self.__init__()

    def compile(self):
        text = self.sbatch_directives()
        if text != "":
            text += "\n"
        text += self.module_import()
        if text != "":
            text += "\n"
        text += self.build_exports()
        if text != "":
            text += "\n"
        if self.include_diagnostics:
            text += self.build_diagnostic()
            text += "\n"
        if self.include_timestart:
            text += self.build_timestart()
            text += "\n"
        text += self.build_commands()
        if text != "":
            text += "\n"
        if self.include_timeend:
            text += self.build_timeend()
        return text

    def format_sbatch(self, varlist):
        frmt = "{0}={1}\n"
        text = ""
        for pair in varlist:
            key = str(pair[0])
            if pair[1] != None:
                value = str(pair[1])
                if key == "":  # comment
                    text += "# " + value + "\n"
                else:
                    text += frmt.format(key, value)
        return text

    def format_section(self, varlist):
        frmt = "{0}={1}\n"
        text = ""
        for pair in varlist:
            key = str(pair[0])
            if pair[1] != None:
                value = str(pair[1])
                if key == "":  # comment
                    text += "# " + value + "\n"
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

    def sbatch_directives(self):
        varlist = [("#SBATCH --job-name", self.job_name),
                   ("#SBATCH --partition", self.partition),
                   ("#SBATCH --nodes", self.nodes),
                   ("#SBATCH --ntasks", self.ntasks),
                   ("#SBATCH --time", self.time),
                   ("#SBATCH --signal", self.signal),
                   ("#SBATCH --output", self.output),
                   ("#SBATCH --workdir", self.workdir),
                   ("#SBATCH --checkpoint", self.checkpoint),
                   ("#SBATCH --checkpoint-dir", self.checkpoint_dir),
                   ("#SBATCH --mail-user", self.mail_user)]

        if (self.mail_types != None) and isinstance(self.mail_types, list):
            for mt in self.mail_types:
                varlist.append(("#SBATCH --mail-type", mt))
        if (self.comments != None) and isinstance(self.comments, list):
            for com in self.comments:
                varlist.append(("#SBATCH --comment", com))

        if self.check_empty(varlist):
            return ""

        text = "#!/bin/sh\n"
        text += self.format_sbatch(varlist)
        return text

    def module_import(self):
        statements = []
        if ((self.modulepath_prepends != None) and
                isinstance(self.modulepath_prepends, list)):
            for prepend in self.modulepath_prepends:
                statements.append("MODULEPATH=" + str(prepend) + ":$MODULEPATH")
        if ((self.modulepath_appends != None) and
                isinstance(self.modulepath_appends, list)):
            for append in self.modulepath_appends:
                statements.append("MODULEPATH=$MODULEPATH:" + str(append))
        if (self.modules != None) and isinstance(self.modules, list):
            for module in self.modules:
                statements.append("module load " + str(module))

        if not statements:
            return ""

        text = "# Modify module environment\n"
        text += "\n".join(statements)
        text += "\n"
        return text

    def build_exports(self):
        statements = []
        if (self.exports != None) and isinstance(self.exports, list):
            for exp in self.exports:
                statements.append("export " + str(exp))

        if not statements:
            return ""

        text = "# Define variable exports\n"
        text += "\n".join(statements)
        text += "\n"
        return text

    def build_diagnostic(self):
        text = "# Call diagnostics, a script that checks env info\n"
        text += "sleep 1\n"
        text += "diagnostics\n"
        return text

    def build_timestart(self):
        jobstatus = None
        jobid = ""
        if self.job_name != None:
            jobid = str(self.job_name) + " "
        if self.status_filename != None:
            jobstatus = str(self.status_filename) + "\n"

        text = "# Script is starting\n"
        if jobstatus != None:
            text += 'echo "SUBMITTED ' + jobid + '`date`" >> ' + jobstatus
        text += 'echo "Job Started at `date`"\n'
        return text

    def build_commands(self):
        statements = []
        if (self.commands != None) and isinstance(self.commands, list):
            for cmd in self.commands:
                statements.append(str(cmd))

        if not statements:
            return ""

        text = "# Perform desired tasks\n"
        text += "\n".join(statements)
        text += "\n"
        return text

    def build_timeend(self):
        mark = "#" * 40
        jobstatus = None
        jobid = ""
        if self.job_name != None:
            jobid = str(self.job_name) + " "
        if self.status_filename != None:
            jobstatus = str(self.status_filename) + "\n"

        text = "# Script is finishing\n"
        if jobstatus != None:
            text += 'echo "FINISHED ' + jobid + '`date`" >> ' + jobstatus
        text += 'echo\n'
        text += 'echo "Job Ended at `date`"\n'
        text += 'echo "' + mark + '"\n'
        return text

if __name__ == "__main__":
    # test
    slurm = Slurm()
    slurm.status_filename = "jobstatus.txt"

    slurm.job_name = "1meth"
    slurm.partition = "serial"
    slurm.nodes = 1
    slurm.ntasks = 10
    slurm.time = "2-12:00:00"
    slurm.signal = 15
    slurm.output = "howdy/test.out"
    slurm.workdir = "howdy"
    slurm.mail_types.append("REQUEUE")
    slurm.mail_types.append("END")
    slurm.comments.append("#SBATCH --array-0,1,3")

    slurm.modulepath_prepends.append("/h3/n1/shirtsgroup/modules")
    slurm.modulepath_prepends.append("$HOME/modules")
    slurm.modules.append("jtp4kc")
    slurm.modules.append("anaconda-jtp")
    slurm.modules.append("gromacs-shirtsgroup/4.6.7")

    slurm.exports.append('THREADINFO="-nt {ntasks} "')
    slurm.exports.append('GMX_SUPPRESS_DUMP=1 #prevent step file output')

    slurm.commands.append('echo "hi"')
    slurm.commands.append('echo "Hello"')
    slurm.commands.append('echo "world"')
    print(slurm.compile())










