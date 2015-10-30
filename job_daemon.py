#!/usr/local/bin/python2.7
# encoding: utf-8
'''
job_daemon -- monitors, restarts, and gives updates on running jobs on Rivanna

@author:     jtp4kc, Tyler P

@copyright:  2015 Shirts Group, University of Virginia. All rights reserved.

@contact:    jtp4kc@virginia.edu
@deffield    updated: 2015-10-11
'''

import sys
import os
import math
import job_save
import traceback
import time
import datetime

from argparse import ArgumentParser
# from argparse import RawDescriptionHelpFormatters

__all__ = []
__version__ = 0.1
__date__ = '2015-10-11'
__updated__ = '2015-10-11'

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

class Attr():
    def __init__(self):
        self.NUM_STEPS = "num_of_steps"
        self.OLD_STEPS = "oldnum_steps"
        self.STATUS = "run_status"
        self.TIME = "last_checked"
ATTR = Attr()

def submit_self(jobname, command):
    fields = dict()
    fields['job-name'] = jobname
    fields['partition'] = 'serial'
    fields['nnodes'] = 1
    fields['ntasks'] = 1
    fields['queue-time'] = '7-00:00:00'
    fields['outputfile'] = 'daemon.log'
    fields['workdir'] = "."
    text = """#!/bin/sh
#SBATCH --job-name={job-name}-daemon
#SBATCH --partition={partition}
#SBATCH --nodes={nnodes}
#SBATCH --ntasks={ntasks}
#SBATCH --time={queue-time}
#SBATCH --signal=15 --comment="15 = SIGTERM"
#SBATCH --output={outputfile}
#SBATCH --workdir={workdir}
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jtp4kc@virginia.edu

MODULEPATH=/h3/n1/shirtsgroup/modules:$MODULEPATH
MODULEPATH=$HOME/modules:$MODULEPATH
module load jtp4kc
module load anaconda-jtp
module load gromacs-shirtsgroup/4.6.7

sleep 1
diagnostics

""".format(**fields)
    text += command

def reschedule(savefilename, pathtohere):
    cmd = "python " + pathtohere + " --save " + savefilename
    submit_self(cmd)

def _filename(entry, key, index=0):
    filename = entry.files[key]
    if isinstance(filename, list):
        filename = filename[index]
    else:
        filename = str(filename)
    return filename

def timecheck(entry):
    mindelta = 60 * 10  # ten minutes
    timestamp = job_save.SerialDate.deserialize(entry.attr[ATTR.TIME])
    currenttime = datetime.datetime.now()
    relative = currenttime - timestamp
    while relative.total_seconds() < mindelta:
        entryname = entry.name
        timetosleep = mindelta - relative.total_seconds()
        print("Sleep " + str(timetosleep) + "s to wait for entry " + entryname)
        print("\tasleep at " + str(datetime.datetime.now()))
        time.sleep(timetosleep)
        currenttime = datetime.datetime.now()
        print("\tawake at " + str(currenttime))
        relative = currenttime - timestamp
    entry.attr[ATTR.TIME] = job_save.SerialDate.serialize(currenttime)

def check_log(entry):
    logfile = _filename(entry, "log")
    logscan = job_save.LogScan(logfile)
    logscan.scan()
    if ATTR.NUM_STEPS in entry.attr:
        entry.attr[ATTR.OLD_STEPS] = entry.attr[ATTR.NUM_STEPS]
    else:
        entry.attr[ATTR.OLD_STEPS] = 0
    entry.attr[ATTR.NUM_STEPS] = logscan.get_step_number()
    return logscan

def check_errors(entry, logscan):
    outfile = _filename(entry, "out")
    outscan = job_save.OutputScan(outfile)
    outscan.scan()
    errors = []

    # error = (name, status, message)
    errors.append(("Warning", outscan.warning_detected, outscan.warn_statement))
    errors.append(("Finished", logscan.finish_detected, ""))
    errors.append(("Cancelled", logscan.cancel_detected or
                   outscan.cancel_detected, outscan.cancel_statement))
    errors.append(("Fault", outscan.fault_detected, outscan.fault_statement))
    errors.append(("Failed", logscan.failure_detected, logscan.fail_statement))
    return errors

def check_resubmit(entry, errors):
    # a fault should always restart from old cpt, in the hopes that it won't
    #    fault again
    # a cancel should only restart if it was due to time limit, since the job
    #    should reasonably still be able to run, run from current cpt
    # warnings don't require action
    # finish means the job should be removed from the daemon's list, after eval
    # a failure is non-recoverable, in almost all cases

    rsc_cpt = False
    rsc_old = False
    status = "Running"
    for (ename, state, message) in errors:
        if state:
            status = ename
            if status == "Cancelled":
                if "due to time limit" in message:
                    rsc_cpt = True
            if status == "Failed":
                rsc_old = True

    entry.attr[ATTR.STATUS] = status
    return status, rsc_old, rsc_cpt

def daemon(savefilename):
    savemgr = job_save.SaveJobs()
    savemgr.load(savefilename)

    daemon_cancel = False
    entries = savemgr.get_jobs()
    if not entries:
        daemon_cancel = True

    while not daemon_cancel:
        savemgr.save(savefilename)
        for entry in entries:
            timecheck(entry)
            try:
                # get log progress
                logscan = check_log(entry)
                # evaluate errors
                error_list = check_errors(entry, logscan)
                # change and resubmit if needed
                (status, oldcpt, cpt) = check_resubmit(entry, error_list)


            except Exception as e:
                traceback.print_exc()
                print("Error occurred while processing job " + entry.jobname)
                print(e)

        # what should the daemon do for each job?
        #    check log to get progress
        #    check output for errors
        #    if errors, fix, change slurm if needed, and resubmit
        #
        #    if no errors, proceed
        #    copy log, xtc, tpr, and xvg to alternate location
        #    process xtc, using visualizer to produce vmd markup
        #    run alchemical analysis to look at weights
        #
        # at end of each cycle, print a summary file showing status of each
        #    job, such as progress, errors/resets, etc
        # if anything is strange, such as temp, pressure, etc, produce error
        #    file for user to read

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

  Created by jtp4kc on %s.
  Copyright 2015 Shirts Group. All rights reserved.

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
        parser.add_argument('-V', '--version', action='version',
            version=program_version_message)
        parser.add_argument('-s', '--save', help="(required) File generated" +
            " by run of genscript, which indicates the files that were" +
            " generated by that run. These jobs are the ones that will be" +
            " monitored and analyzed.")
        parser.add_argument(dest="paths", help="paths to folder(s) with" +
            " source file(s) [default: %(default)s]", metavar="path", nargs='*')

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
            raise CLIError("include and exclude pattern are equal! Nothing" +
                " will be processed.")

        for inpath in paths:
            ### do something with inpath ###
            print(inpath)

        save_name = None
        if args.save:
            save_name = args.save

        if save_name:
            daemon(save_name)
        else:
            parser.print_usage()

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
        profile_filename = 'job_daemon_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())
