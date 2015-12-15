#!/usr/local/bin/python2.7
# encoding: utf-8
'''
job_daemon -- monitors, restarts, and gives updates on running jobs on Rivanna

@author:     jtp4kc, Tyler P

@copyright:  2015 Shirts Group, University of Virginia. All rights reserved.

@contact:    jtp4kc@virginia.edu
@deffield    updated: 2015-11-05
'''

import os
import sys
import job_utils
import traceback
import tempfile
import time
import datetime
import signal
import backup
import subprocess
import visualizer
import numpy

from argparse import ArgumentParser
from slurm_template import Slurm
from audioop import avg
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

class SignalError(Exception):
    '''Exception to raise when a signal is encountered.'''
    def __init__(self, signum=-1):
        super(SignalError).__init__(type(self))

        self.signum = signum
        msg = "SIGNAL " + str(signum)
        if signum == 14:
            msg = "14:SIGALRM"
        elif signum == 15:
            msg = "15:SIGTERM"
        elif signum == 16:
            msg = "16:SIGUSR1"

        self.msg = "E: %s" % msg
    def __str__(self):
        return self.msg
    def __unicode__(self):
        return self.msg

class Attr():
    def __init__(self):
        # can't contain '-' or other special characters
        self.NUM_STEPS = "num_of_steps"
        self.OLD_STEPS = "oldnum_steps"
        self.VIS_MARK = "steps_visualized"
        self.STATUS = "run_status"
        self.TIME = "last_checked"
        self.IGNORE = "daemon_ignore"
        self.JOB_ID = "job_id"
        self.LOG_NSTEP = "log_steps_per_frame"
        self.LOG_DELTA = "log_simtime_delta"
        self.RESNAME = "residue_name"
ATTR = Attr()

def get_slurm(jobname, command, time="7-00:00:00", output="daemon.log"):
    slurm = Slurm()
    slurm.job_name = jobname + "-daemon"
    slurm.partition = "serial"
    slurm.nodes = 1
    slurm.ntasks = 1
    slurm.time = time
    slurm.output = output
    slurm.signal = 15  # 14 = SIGALRM, 15 = SIGTERM, 16 = SIGUSR1
    slurm.mail_types.extend(["REQUEUE", "END", "FAIL"])
    slurm.mail_user = "jtp4kc@virginia.edu"

    slurm.modulepath_prepends.extend(["/h3/n1/shirtsgroup/modules",
                                      "$HOME/modules"])
    slurm.modules.extend(["jtp4kc", "anaconda-jtp", "gromacs-shirtsgroup/4.6.7"])

    slurm.commands.append(command)
    return slurm

def reschedule_self(jobname, savefilename, pathtohere=None, time=None,
                    live=False):
    """ Submit a slurm script that will run this daemon with the same setup
    @return: (str) job id number from slurm submission
    """
    if pathtohere == None:
        pathtohere = __file__
    if pathtohere.endswith(".pyc"):
        pathtohere = pathtohere.replace(".pyc", ".py")
    fname = savefilename
    if "." in fname:
        fname = "".join(fname.split(".")[:-1]) + ".slurm"
    else:
        fname += ".slurm"
    opdir = os.path.dirname(fname)
    fname = os.path.join(opdir, "daemon-" + os.path.basename(fname))
    count = 1
    outname = "daemon-" + jobname + "-1.log"
    outpath = os.path.join(opdir, outname)
    while os.path.exists(outpath):
        count += 1
        outname = "daemon-" + jobname + "-" + str(count) + ".log"
        outpath = os.path.join(opdir, outname)

    extr = ""
    if live:
        extr = " --live"

    cmd = ("python -u " + pathtohere + " " + os.path.basename(savefilename) +
           extr)
    if time != None:
        slurm = get_slurm(jobname, cmd, time, outname)
    else:
        slurm = get_slurm(jobname, cmd, output=outname)
    return job_utils.submit_slurm(slurm, fname)

def _filename(entry, key, index=0):
    filename = entry.files[key]
    if isinstance(filename, list):
        filename = filename[index]
    else:
        filename = str(filename)
    return filename

def setup_signalhandler():
    def handler(signum, _):  # signum, frame
        print("Signal encountered: " + str(signum))
        raise(SignalError(signum))

    signal.signal(signal.SIGALRM, handler)  # @UndefinedVariable
    signal.signal(signal.SIGTERM, handler)
    signal.signal(signal.SIGUSR1, handler)  # @UndefinedVariable

def print_entry_messages(savefilename, jobname, savemgr, timestamp):
    print("Daemon start - " + timestamp.isoformat())
    print("    file: " + savefilename)
    print("    jobn: " + jobname)
    print("    numj: " + str(len(savemgr.jobs)))
    for job in savemgr.get_jobs():
        print("       -> " + job.jobname)

def timecheck(entry):
    mindelta = 60 * 10  # ten minutes
    currenttime = datetime.datetime.now()
    if ATTR.TIME in entry.attr:
        timestamp = job_utils.SerialDate.deserialize(entry.attr[ATTR.TIME])
        relative = currenttime - timestamp
        if relative.total_seconds() < 0:
            print("Error: timestamp on file is newer than current time")
            print("stamp: " + str(timestamp))
            print("ctime: " + str(currenttime))
            print(" diff: " + str(relative))
        else:
            while relative.total_seconds() < mindelta:
                entryname = entry.jobname
                timetosleep = mindelta - relative.total_seconds()
                print("Sleep " + str(timetosleep) + "s to wait for entry " +
                      entryname)
                print("\tasleep at " + str(datetime.datetime.now()))
                time.sleep(timetosleep)
                currenttime = datetime.datetime.now()
                print("\tawake at " + str(currenttime))
                relative = currenttime - timestamp
    entry.attr[ATTR.TIME] = job_utils.SerialDate.serialize(currenttime)

def check_log(entry):
    logfile = _filename(entry, "log")
    logscan = job_utils.LogScan(logfile)
    logscan.scan()
    if ATTR.NUM_STEPS in entry.attr:
        entry.attr[ATTR.OLD_STEPS] = entry.attr[ATTR.NUM_STEPS]
    else:
        entry.attr[ATTR.OLD_STEPS] = "0"
    entry.attr[ATTR.NUM_STEPS] = str(logscan.get_step_number())
    return logscan

def check_errors(entry, logscan):
    outfile = _filename(entry, "out")
    outscan = job_utils.OutputScan(outfile)
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
    msg = ""
    for (ename, state, message) in errors:
        if state:
            status = ename
            msg = message
            if status == "Cancelled":
                if "DUE TO TIME LIMIT" in message:
                    rsc_cpt = True
            if status == "Fault":
                rsc_old = True

    entry.attr[ATTR.STATUS] = status
    return status, rsc_old, rsc_cpt, msg

def resubmit_job(entry, live, prev=False):
    slurmfile = _filename(entry, "slurm")
    tmpfile = slurmfile + ".tmp"
    outfile = _filename(entry, "out")
    tprfile = _filename(entry, "tpr")

    old_id = entry.attr[ATTR.JOB_ID]
    found = False
    file_ = tempfile.NamedTemporaryFile(mode="w+t", prefix='jdr', dir='.')
    subprocess.call(["squeue"], stdout=file_)
    file_.seek(0)
    for line in file_:
        if old_id in line:
            if " CG " in line:
                # this job appears to have gotten stuck, and is not valid
                print("Job " + old_id + " appears to be invalid, resubmitting.")
            else:
                found = True
            break
    file_.close()
    if found:
        print("Waiting for job already submitted #" + old_id + ", before" +
              "attempting to resubmit.")
        return

    if prev:
        cptfile = _filename(entry, "prev")
    else:
        cptfile = _filename(entry, "cpt")
    cptfpath = os.path.relpath(cptfile, os.path.dirname(slurmfile))

    if live:
        backup.backup_file(".", outfile, copy=True)
        backup.backup_file(".", slurmfile, copy=True)
    else:
        print("DRYRUN: Would backup file " + outfile)
        print("DRYRUN: Would backup file " + slurmfile)

    if live:
        writer = open(tmpfile, "w")
    else:
        class Facade(list):
            def __init__(self):
                super(list)
                self.write = self.append

            def close(self):
                pass
        writer = Facade()  # create a false, file-like object
        #    inspired by xml ElementTree.toString

    wkdr = None
    for line in open(slurmfile, "r"):
        tname = os.path.basename(tprfile)
        if ("SBATCH" in line) and ("workdir=" in line):
            offset = line.find("workdir=")
            half = line.replace("workdir=", "")[offset:]
            wkdr = half.split()[0]
        elif line.strip().startswith("rm ") and (tname in line):
            pass
        elif ("grompp" in line):
            pass
        elif ("mdrun" in line):
            args = line.split()
            start = line.find(args[0])  # we know args[0] is in line
            wspc = line[:start]  # leading whitespace
            if "-cpi" in args:
                i = args.indexOf("-cpi")
                args = args[:i] + args[(i + 2):]
            if wkdr:
                cptfile = os.path.join(wkdr, cptfile)
            cmd = wspc + " ".join(args) + " -cpi " + cptfpath + "\n"
            writer.write(cmd)
        else:
            writer.write(line)
    writer.close()

    if live:
        print("Slurm code copied to temporary file " + tmpfile)
    else:
        print("DRYRUN: Slurm code would be written as the following (tab" +
              " indented):")
        for line in writer:
            print("\t" + line.rstrip())

    if live:
        os.rename(tmpfile, slurmfile)
        entry.attr[ATTR.JOB_ID] = job_utils.submit_job(slurmfile, entry.jobname)
    else:
        print("DRYRUN: Would move " + tmpfile + " to " + slurmfile)
        print("DRYRUN: Would sbatch " + slurmfile + " and store job id")

def analyze_job(entry, logscan, is_being_removed, live):
    tprfile = _filename(entry, "tpr")
    xtcfile = _filename(entry, "xtc")
    xvgfile = _filename(entry, "xvg")

    resname = entry.attr[ATTR.RESNAME]
    dt = float(entry.attr[ATTR.LOG_DELTA])
    frames_per_ns = int(1000 / dt)

    nstep = int(entry.attr[ATTR.LOG_NSTEP])
    last = int(entry.attr[ATTR.OLD_STEPS]) / nstep
    now = int(entry.attr[ATTR.NUM_STEPS]) / nstep
    delta = now - last
    nanoseconds = int(now / frames_per_ns)  # how many ns have passed?
    end = nanoseconds * frames_per_ns  # how many frames for this many ns?
    # this rounds the current end value to the nearest ns

    if ATTR.VIS_MARK not in entry.attr:
        entry.attr[ATTR.VIS_MARK] = 0

    print(("Analysis, step numbers: dt|{0:0.2f} fpn|{1:0.0f} log#steps|" +
           "{2:0.0f} last|{3:0.0f} now|{4:0.0f} ns|{5:0.0f} vis|{6:0.0f}" +
           " end|{7:0.0f}").format(dt, frames_per_ns, nstep, last, now,
           nanoseconds, entry.attr[ATTR.VIS_MARK], end))

    message1 = None
    message2 = None
    if "Count" in logscan.log_entries:
        data = logscan.log_numbers["Count"]
        count = []
        for number in data:
            count.append(int(number))
        if len(count) > 0:
            count = numpy.array(count, numpy.float32)
            message1 = "[" + numpy.array2string(count, precision=3,
                                         separator=", ") + "]"
            avg = numpy.average(count)
            count /= avg
            count /= numpy.min(count)
            message2 = "[" + numpy.array2string(count, precision=3,
                                         separator=", ") + "]"

    if message1 is not None:
        print("Current number of samples:")
        print(message1)
        print("Current ratio of samples:")
        print(message2)

    if delta <= 0:
        print("Entry " + entry.jobname + " appears to have made no progress")
        print("Old: {0}, New: {1}, Delta: {2}".format(last, now, delta))

    if (entry.attr[ATTR.VIS_MARK] < end) or is_being_removed:
        workdir = os.path.join("daemon-vis", "")
        fldrname = os.path.basename(os.path.dirname(xtcfile))
        newfldr = os.path.join(workdir, fldrname, "")

        if not os.path.exists(workdir):
            if live:
                os.mkdir(workdir)
            else:
                print("DRYRUN: Would create folder " + workdir)

        if not os.path.exists(newfldr):
            if live:
                os.mkdir(newfldr)
            else:
                print("DRYRUN: Would create folder " + newfldr)

        for fname in [tprfile, xtcfile, xvgfile]:
            cmd = "cp " + fname + " " + newfldr
            if live:
                print(cmd)
                os.system(cmd)
            else:
                print("DRYRUN: Would " + cmd)

        tpr2 = os.path.join(newfldr, os.path.basename(tprfile))
        xtc2 = os.path.join(newfldr, os.path.basename(xtcfile))
        xvg2 = os.path.join(newfldr, os.path.basename(xvgfile))

        start = entry.attr[ATTR.VIS_MARK]
        if is_being_removed:
            end = now
        length = end - start

        if live:
            visualizer.visualize(tpr2, xtc2, xvg2, resname,
                nstart=start, nlength=length, doCenter=True, doVMD=True,
                timedelta=dt)
        else:
            print("DRYRUN: Would make a call to visualizer to generate vmd" +
                  " instructions using " + xtc2)
        entry.attr[ATTR.VIS_MARK] = end
    else:
        print("Waiting for a nanosecond of simulation to pass before" +
              " visualizing.")

def check_cancel():
    return os.path.exists("daemon-cancel.txt")

def daemon(savefilename, live=False):
    timestamp = datetime.datetime.now()

    savemgr = job_utils.SaveJobs()
    savemgr.load(savefilename)

    if "name" in savemgr.attr:
        jobname = savemgr.attr["name"]
    else:
        jobname = "(unknown)"

    print_entry_messages(savefilename, jobname, savemgr, timestamp)

    daemon_cancel = False
    entries = savemgr.get_jobs()
    if len(entries) == 0:
        print("No jobs to monitor found in save file.")
        daemon_cancel = True
    else :
        setup_signalhandler()

    try:
        while not daemon_cancel:
            if live:
                print("Saving job status to file...")
                savemgr.save()
            else:
                print("DRYRUN: Would save job status")
            for entry in entries[:]:  # make a slice copy, entries might change
                # a previous daemon has already marked this file as ignore
                if ATTR.IGNORE in entry.attr:
                    if entry.attr[ATTR.IGNORE] == "true":
                        entries.remove(entry)
                        continue

                # don't spin our wheels, see if the daemon should sleep
                timecheck(entry)
                try:
                    # get log progress
                    logscan = check_log(entry)
                    # evaluate errors
                    err_list = check_errors(entry, logscan)
                    # change and resubmit if needed
                    (status, oldcpt, cpt, out) = check_resubmit(entry, err_list)

                    msg = " is believed to have encountered a "
                    if status.endswith("ed"):
                        msg = " is believed to have "
                    if status == "Running":
                        msg = " is believed to be "
                    if status == "Cancelled":
                        msg = " is believed to have been "
                    print("Job " + entry.jobname + msg + status)
                    out = out.splitlines()
                    for line in out:
                        print("\t" + line)

                    remove_job = False
                    if (status == "Finished") or (status == "Failed"):
                        # to be removed after analysis
                        remove_job = True

                    # reschedule
                    if oldcpt:
                        resubmit_job(entry, live, prev=True)
                    elif cpt:
                        resubmit_job(entry, live, prev=False)
                    else:  # analyze
                        if status == "Cancelled":  # if not reschedule
                            remove_job = True
                        analyze_job(entry, logscan, remove_job, live)

                    if remove_job:
                        # mark as ignore and remove from current list
                        entry.attr[ATTR.IGNORE] = "true"
                        entries.remove(entry)
                except Exception as e:
                    if isinstance(e, SignalError):
                        raise e
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
                # at end of each cycle, print a summary file showing status of
                #    each job, such as progress, errors/resets, etc
                # if anything is strange, such as temp, pressure, etc, produce
                #    error file for user to read

            if len(entries) == 0:
                print("All jobs have been removed from the daemon queue.")
                daemon_cancel = True

            if check_cancel():
                print("Cancel requested by user.")
                daemon_cancel = True

            runtime = datetime.datetime.now() - timestamp
            if runtime.total_seconds() > (7 * 24 * 60 * 60):  # 1 week
                # As a backup, use the total runtime as a signal to stop,
                #    since on the cluster max runtime is currently 7 days
                # If this is still running, it probably shouldn't be.
                print("Maximum execution time encountered (1wk), stopping.")
                daemon_cancel = True
    except SignalError as sig:
        if sig.signum == 15:
            if check_cancel():
                # user requested cancel
                print("Cancel requested by user")
            elif live:
                print("Rescheduling due to signal.")
                reschedule_self(jobname, savefilename, live=live)
            else:
                print("DRYRUN: Would reschedule self, but dry-mode is active.")
        else:
            print("Daemon stopping due to signal.")
        daemon_cancel = True
    print("Daemon exit - " + datetime.datetime.now().isoformat())

def generate_runentry(entry):
    log = entry.files["log"]
    out = entry.files["out"]
    slurm = entry.files["slurm"]

    writer = open(log, "w")
    writer.write(textfor_log(running=True))
    writer.close()

    writer = open(out, "w")
    writer.write(textfor_out(running=True))
    writer.close()

    writer = open(slurm, "w")
    writer.write(textfor_slurm(entry.jobname.replace("test-", "")))
    writer.close()

def generate_warnentry(entry):
    log = entry.files["log"]
    out = entry.files["out"]
    slurm = entry.files["slurm"]

    writer = open(log, "w")
    writer.write(textfor_log(running=True))
    writer.close()

    writer = open(out, "w")
    writer.write(textfor_out(running=True, warn=True))
    writer.close()

    writer = open(slurm, "w")
    writer.write(textfor_slurm(entry.jobname.replace("test-", "")))
    writer.close()

def generate_finentry(entry):
    log = entry.files["log"]
    out = entry.files["out"]
    slurm = entry.files["slurm"]

    writer = open(log, "w")
    writer.write(textfor_log())
    writer.close()

    writer = open(out, "w")
    writer.write(textfor_out())
    writer.close()

    writer = open(slurm, "w")
    writer.write(textfor_slurm(entry.jobname.replace("test-", "")))
    writer.close()

def generate_cnclentry(entry):
    log = entry.files["log"]
    out = entry.files["out"]
    slurm = entry.files["slurm"]

    writer = open(log, "w")
    writer.write(textfor_log(cancel=True))
    writer.close()

    writer = open(out, "w")
    writer.write(textfor_out(cancel=True))
    writer.close()

    writer = open(slurm, "w")
    writer.write(textfor_slurm(entry.jobname.replace("test-", "")))
    writer.close()

def generate_faltentry(entry):
    log = entry.files["log"]
    out = entry.files["out"]
    slurm = entry.files["slurm"]

    writer = open(log, "w")
    writer.write(textfor_log())
    writer.close()

    writer = open(out, "w")
    writer.write(textfor_out(fault=True))
    writer.close()

    writer = open(slurm, "w")
    writer.write(textfor_slurm(entry.jobname.replace("test-", "")))
    writer.close()


def generate_failentry(entry):
    log = entry.files["log"]
    out = entry.files["out"]
    slurm = entry.files["slurm"]

    writer = open(log, "w")
    writer.write(textfor_log(fail=True))
    writer.close()

    writer = open(out, "w")
    writer.write(textfor_out())
    writer.close()

    writer = open(slurm, "w")
    writer.write(textfor_slurm(entry.jobname.replace("test-", "")))
    writer.close()

def textfor_log(running=False, cancel=False, fail=False):
    text = """DD  step 14096499  vol min/aver 0.855  load imb.: force  1.2%  pme mesh/force 0.095

           Step           Time         Lambda
       14096500    28193.00000        0.00000

             MC-lambda information
  N   FEPL  CoulL   VdwL    Count   G(in kT)  dG(in kT)
  1  0.000  0.000  0.000    65995    0.00000    0.21165
  2  0.000  0.300  0.000    30003    0.21165    0.57929
  3  0.000  0.500  0.000    31519    0.79094    0.97076
  4  0.000  1.000  0.000    31283    1.76170    1.05917
  5  0.000  1.000  0.100    33674    2.82087    1.02449
  6  0.000  1.000  0.200    36832    3.84536    1.01499
  7  0.000  1.000  0.300    38029    4.86035    0.54104
  8  0.000  1.000  0.350    39562    5.40139    0.47087
  9  0.000  1.000  0.400    38420    5.87226    0.47913
 10  0.000  1.000  0.450    37423    6.35139    0.41223
 11  0.000  1.000  0.500    34736    6.76362    0.28777
 12  0.000  1.000  0.530    35423    7.05139    0.20000
 13  0.000  1.000  0.570    31695    7.25139    0.03752
 14  0.000  1.000  0.600    27613    7.28891   -0.03752
 15  0.000  1.000  0.630    23829    7.25139   -0.15000
 16  0.000  1.000  0.670    21001    7.10139   -0.26714 <<
 17  0.000  1.000  0.700    19105    6.83425   -0.48286
 18  0.000  1.000  0.730    16394    6.35139   -0.75000
 19  0.000  1.000  0.770    15445    5.60139   -0.70304
 20  0.000  1.000  0.800    14321    4.89835   -0.64696
 21  0.000  1.000  0.830    14289    4.25139   -0.75000
 22  0.000  1.000  0.870    14835    3.50139   -0.56082
 23  0.000  1.000  0.900    13445    2.94057    0.43249
 24  0.000  1.000  1.000    39953    3.37306    0.00000

   Energies (kJ/mol)
           Bond          Angle    Proper Dih. Ryckaert-Bell.          LJ-14
    2.04772e+03    5.34532e+03    2.89774e+02    3.97123e+03    2.38392e+03
     Coulomb-14        LJ (SR)  Disper. corr.   Coulomb (SR)   Coul. recip.
    2.24794e+04    3.20601e+04   -2.19726e+03   -2.66631e+05   -7.75192e+04
      Potential    Kinetic En.   Total Energy  Conserved En.    Temperature
   -2.77770e+05    5.38398e+04   -2.23931e+05   -2.25970e+05    2.96962e+02
 Pres. DC (bar) Pressure (bar)    dVremain/dl      dVcoul/dl       dVvdw/dl
   -1.72064e+02   -1.68627e+02    0.00000e+00    5.54800e+00    6.58150e+00

Writing checkpoint, step 14096940 at Thu Oct 29 18:46:40 2015


DD  step 14096999  vol min/aver 0.866  load imb.: force  1.0%  pme mesh/force 0.095

           Step           Time         Lambda
       14097000    28194.00000        0.00000

             MC-lambda information
  N   FEPL  CoulL   VdwL    Count   G(in kT)  dG(in kT)
  1  0.000  0.000  0.000    65995    0.00000    0.21165
  2  0.000  0.300  0.000    30003    0.21165    0.57929
  3  0.000  0.500  0.000    31519    0.79094    0.97076
  4  0.000  1.000  0.000    31283    1.76170    1.05917
  5  0.000  1.000  0.100    33674    2.82087    1.02449
  6  0.000  1.000  0.200    36832    3.84536    1.01499
  7  0.000  1.000  0.300    38031    4.86035    0.54104
  8  0.000  1.000  0.350    39565    5.40139    0.47087
  9  0.000  1.000  0.400    38425    5.87226    0.47913 <<
 10  0.000  1.000  0.450    37425    6.35139    0.41223
 11  0.000  1.000  0.500    34741    6.76362    0.28777
 12  0.000  1.000  0.530    35424    7.05139    0.20000
 13  0.000  1.000  0.570    31696    7.25139    0.03752
 14  0.000  1.000  0.600    27616    7.28891   -0.03752
 15  0.000  1.000  0.630    23830    7.25139   -0.15000
 16  0.000  1.000  0.670    21002    7.10139   -0.26714
 17  0.000  1.000  0.700    19106    6.83425   -0.48286
 18  0.000  1.000  0.730    16394    6.35139   -0.75000
 19  0.000  1.000  0.770    15445    5.60139   -0.70304
 20  0.000  1.000  0.800    14321    4.89835   -0.64696
 21  0.000  1.000  0.830    14289    4.25139   -0.75000
 22  0.000  1.000  0.870    14835    3.50139   -0.56082
 23  0.000  1.000  0.900    13445    2.94057    0.43249
 24  0.000  1.000  1.000    39953    3.37306    0.00000

   Energies (kJ/mol)
           Bond          Angle    Proper Dih. Ryckaert-Bell.          LJ-14
    2.04753e+03    5.12697e+03    3.26790e+02    4.00911e+03    2.32533e+03
     Coulomb-14        LJ (SR)  Disper. corr.   Coulomb (SR)   Coul. recip.
    2.24557e+04    3.17309e+04   -2.19836e+03   -2.65506e+05   -7.75150e+04
      Potential    Kinetic En.   Total Energy  Conserved En.    Temperature
   -2.77197e+05    5.40505e+04   -2.23147e+05   -2.25977e+05    2.98124e+02
 Pres. DC (bar) Pressure (bar)    dVremain/dl      dVcoul/dl       dVvdw/dl
   -1.72150e+02   -3.42048e+02    0.00000e+00   -9.36929e-01    3.97799e+01


"""
    if fail:
        text += """Fatal error: Particles are more than 2/3 the cutoff 
distance outside of the domain decomposition.
"""
    if cancel:
        text += """Received the TERM signal, stopping at the next NS step
"""
    if not running:
        text += """
           Step           Time         Lambda
       14097250    28194.50000        0.00000

Writing checkpoint, step 14097250 at Thu Oct 29 18:46:54 2015


   Energies (kJ/mol)
           Bond          Angle    Proper Dih. Ryckaert-Bell.          LJ-14
    2.02788e+03    5.07956e+03    2.54175e+02    3.99711e+03    2.43670e+03
     Coulomb-14        LJ (SR)  Disper. corr.   Coulomb (SR)   Coul. recip.
    2.24827e+04    3.19662e+04   -2.19857e+03   -2.66080e+05   -7.75138e+04
      Potential    Kinetic En.   Total Energy  Conserved En.    Temperature
   -2.77548e+05    5.39620e+04   -2.23586e+05   -2.25984e+05    2.97636e+02
 Pres. DC (bar) Pressure (bar)    dVremain/dl      dVcoul/dl       dVvdw/dl
   -1.72166e+02   -6.06205e+01    0.00000e+00    1.13615e+00    2.29902e+01

        <======  ###############  ==>
        <====  A V E R A G E S  ====>
        <==  ###############  ======>

        Statistics over 14097251 steps using 14097251 frames

   Energies (kJ/mol)
           Bond          Angle    Proper Dih. Ryckaert-Bell.          LJ-14
    2.05924e+03    5.22387e+03    3.10870e+02    4.00705e+03    2.35929e+03
     Coulomb-14        LJ (SR)  Disper. corr.   Coulomb (SR)   Coul. recip.
    2.23773e+04    3.20246e+04   -2.19828e+03   -2.65990e+05   -7.74680e+04
      Potential    Kinetic En.   Total Energy  Conserved En.    Temperature
   -2.77294e+05    5.43909e+04   -2.22903e+05   -2.24607e+05    3.00001e+02
 Pres. DC (bar) Pressure (bar)    dVremain/dl      dVcoul/dl       dVvdw/dl
   -1.72144e+02   -1.56403e+02    0.00000e+00    4.41300e+00    1.11044e+01

   Total Virial (kJ/mol)
    1.91251e+04   -1.76265e+01   -6.91801e+00
   -1.76265e+01    1.91140e+04   -1.15412e+01
   -6.91801e+00   -1.15412e+01    1.91476e+04

   Pressure (bar)
   -1.55342e+02    4.05102e-01    9.97033e-01
    4.05102e-01   -1.56859e+02    1.25685e+00
    9.97033e-01    1.25685e+00   -1.57008e+02


        M E G A - F L O P S   A C C O U N T I N G

 NB=Group-cutoff nonbonded kernels    NxN=N-by-N cluster Verlet kernels
 RF=Reaction-Field  VdW=Van der Waals  QSTab=quadratic-spline table
 W3=SPC/TIP3p  W4=TIP4p (single or pairs)
 V&F=Potential and force  V=Potential only  F=Force only

 Computing:                               M-Number         M-Flops  % Flops
-----------------------------------------------------------------------------
 NB Elec. [V&F]                   3810965569.444568  3810965569.445    55.7
 NB Generic kernel                1061276257.306932  1061276257.307    15.5
 NB Free energy kernel             47064412.994508    47064412.995     0.7
 1,4 nonbonded interactions           96693.044609     8702374.015     0.1
 Calc Weights                        894682.034715    32208553.250     0.5
 Spread Q Bspline                  38173100.147840    76346200.296     1.1
 Gather F Bspline                  38173100.147840   229038600.887     3.3
 3D-FFT                           172527401.915372  1380219215.323    20.2
 Solve PME                            88417.958272     5658749.329     0.1
 NS-Pairs                           3410304.388011    71616392.148     1.0
 Reset In Box                         10319.194320       30957.583     0.0
 CG-CoM                               29822.774685       89468.324     0.0
 Bonds                                18551.982316     1094566.957     0.0
 Angles                               67032.428505    11261447.989     0.2
 Propers                               7274.181516     1665787.567     0.0
 RB-Dihedrals                         73714.525479    18207487.793     0.3
 Virial                              305839.860445     5505117.488     0.1
 Stop-CM                                  0.021155           0.212     0.0
 Calc-Ekin                           298227.344905     8052138.312     0.1
 Shake                               164362.595293     4930877.859     0.1
 Constraint-V                        298452.900921     2387623.207     0.0
 Shake-Init                           37216.745280      372167.453     0.0
 Constraint-Vir                      298452.900921     7162869.622     0.1
 Settle                              174157.451208    56252856.740     0.8
-----------------------------------------------------------------------------
 Total                                              6840109692.100   100.0
-----------------------------------------------------------------------------


    D O M A I N   D E C O M P O S I T I O N   S T A T I S T I C S

 av. #atoms communicated per step for force:  2 x 45373.2

 Average load imbalance: 1.1 %
 Part of the total run time spent waiting due to load imbalance: 0.8 %
 Steps where the load balancing was limited by -rdd, -rcon and/or -dds: X 0 % Y 0 %
 Average PME mesh/force load: 0.085
 Part of the total run time spent waiting due to PP/PME imbalance: 35.3 %

NOTE: 35.3 % performance was lost because the PME nodes
      had less work to do than the PP nodes.
      You might want to decrease the number of PME nodes
      or decrease the cut-off and the grid spacing.


     R E A L   C Y C L E   A N D   T I M E   A C C O U N T I N G

 Computing:         Nodes   Th.     Count  Wall t (s)     G-Cycles       %
-----------------------------------------------------------------------------
 Domain decomp.        12    1    1409726    1559.595    46674.451     0.2
 DD comm. load         12    1    1409725       9.392      281.081     0.0
 DD comm. bounds       12    1    1409725      26.122      781.771     0.0
 Send X to PME         12    1   14097251     325.825     9751.048     0.0
 Neighbor search       12    1    1409726   13027.650   389882.279     1.3
 Comm. coord.          12    1   14097251    1353.895    40518.419     0.1
 Force                 12    1   14097251  512370.427 15333859.409    50.8
 Wait + Comm. F        12    1   14097251   69387.393  2076576.779     6.9
 PME mesh               8    1   14097251   49020.179   978027.471     3.2
 PME wait for PP        8                  555794.152 11088942.613    36.8
 Wait + Recv. PME F    12    1   14097251     270.195     8086.218     0.0
 Write traj.           12    1      28854      35.819     1071.970     0.0
 Update                12    1   14097251     536.188    16046.663     0.1
 Constraints           12    1   28194502    3057.494    91502.506     0.3
 Comm. energies        12    1   28194502    1206.574    36109.506     0.1
 Rest                  12                     988.678    49314.070     0.2
-----------------------------------------------------------------------------
 Total                 20                  604814.352 30167426.251   100.0
-----------------------------------------------------------------------------
-----------------------------------------------------------------------------
 PME redist. X/F        8    1   42291753    4851.154    96787.934     0.3
 PME spread/gather      8    1   56389004   16383.008   326866.040     1.1
 PME 3D-FFT             8    1   56389004   15843.872   316109.452     1.0
 PME 3D-FFT Comm.       8    1   56389004    3090.970    61669.576     0.2
 PME solve              8    1   28194502    8802.028   175613.909     0.6
-----------------------------------------------------------------------------

               Core t (s)   Wall t (s)        (%)
       Time: 12101659.800   604814.352     2000.9
                         7d00h00:14
                 (ns/day)    (hour/ns)
Performance:        4.028        5.959
Finished mdrun on node 0 Thu Oct 29 18:46:54 2015

"""
    return text

def textfor_out(running=False, cancel=False, fault=False, warn=False):
    text = """starting mdrun 'Protein in water'
25000000 steps,  50000.0 ps.

NOTE: Turning on dynamic load balancing
"""
    if warn:
        text += """
WARNING: Shake did not converge in 1000 steps
"""
    if cancel:
        text += """
slurmstepd: *** JOB 61707 CANCELLED AT 2015-10-29T18:46:54 DUE TO TIME LIMIT ***


Received the TERM signal, stopping at the next NS step
"""
    if fault:
        text += """
Segmentation fault- unable to continue
"""
    if not running:
        text += """


 Average load imbalance: 1.1 %
 Part of the total run time spent waiting due to load imbalance: 0.8 %
 Steps where the load balancing was limited by -rdd, -rcon and/or -dds: X 0 % Y 0 %
 Average PME mesh/force load: 0.085
 Part of the total run time spent waiting due to PP/PME imbalance: 35.3 %

NOTE: 35.3 % performance was lost because the PME nodes
      had less work to do than the PP nodes.
      You might want to decrease the number of PME nodes
      or decrease the cut-off and the grid spacing.


               Core t (s)   Wall t (s)        (%)
       Time: 12101659.800   604814.352     2000.9
                         7d00h00:14
                 (ns/day)    (hour/ns)
Performance:        4.028        5.959

"""
    return text

def textfor_slurm(jobname="1meth-24states"):
    text = """#!/bin/sh
#SBATCH --job-name=""" + jobname + """-00
#SBATCH --partition=serial
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=7-00:00:00
#SBATCH --signal=15 --comment="15 = SIGTERM"
#SBATCH --output=/sfs/lustre/scratch/jtp4kc/simulations/1meth-mstates/""" + jobname + """-00/""" + jobname + """.out
#SBATCH --workdir=/sfs/lustre/scratch/jtp4kc/simulations/1meth-mstates/""" + jobname + """-00/
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
module load anaconda-jtp
module load openmpi/1.8.1
module load gromacs-shirtsgroup/4.6.7

export THREADINFO="-nt 20 "
export GMX_SUPPRESS_DUMP=1 #prevent step file output

sleep 1
diagnostics

#Change the job status to 'SUBMITTED'
echo "SUBMITTED """ + jobname + """-00 `date`" >> ../jobstatus.txt
echo "Job Started at `date`"

#PRODUCTION
rm """ + jobname + """.tpr
grompp_d -c ../t41meth-00-in.gro -p ../t41meth.top -n ../t41meth.ndx -f """ + jobname + """-00.mdp -o """ + jobname + """.tpr -maxwarn 15
if [ -f """ + jobname + """.tpr ];
then
    echo "MD Simulation launching..."
    mdrun_d ${THREADINFO} -deffnm """ + jobname + """
else
    echo "Error generating """ + jobname + """.tpr"
fi

mv /sfs/lustre/scratch/jtp4kc/simulations/1meth-mstates/""" + jobname + """.par """ + jobname + """.par



echo "FINISHED """ + jobname + """-00 `date`" >> ../jobstatus.txt
# print end time
echo
echo "Job Ended at `date`"
echo "###################################################################"
"""
    return text

def test_daemon(savename):
    # ensure unique name
    orig = savename.replace(".save", "")
    workdir = orig
    count = 0
    while os.path.exists(savename) or os.path.exists(workdir):
        count += 1
        workdir = orig + "_td" + str(count)
        savename = workdir + ".save"

    os.mkdir(workdir)
    savemgr = job_utils.SaveJobs()
    savemgr.attr["name"] = "test-daemon"

    run_entry = job_utils.SaveEntry()
    warn_entry = job_utils.SaveEntry()
    fin_entry = job_utils.SaveEntry()
    cncl_entry = job_utils.SaveEntry()
    falt_entry = job_utils.SaveEntry()
    fail_entry = job_utils.SaveEntry()

    entries = [run_entry, warn_entry, fin_entry, cncl_entry, falt_entry,
               fail_entry]
    for entry in entries:
        savemgr.add_job(entry)
        entry.attr[ATTR.JOB_ID] = "900001"
        entry.attr[ATTR.LOG_DELTA] = "1.0"
        entry.attr[ATTR.RESNAME] = "TMP"

    run_entry.jobname = "test-run"
    warn_entry.jobname = "test-warning"
    fin_entry.jobname = "test-finish"
    cncl_entry.jobname = "test-cancel"
    falt_entry.jobname = "test-fault"
    fail_entry.jobname = "test-fail"

    exts = ["log", "out", "tpr", "xvg", "xtc", "cpt", "slurm"]
    for entry in entries:
        prefix = entry.jobname.split("-")[1]
        for ext in exts:
            entry.files[ext] = os.path.join(workdir, prefix + "." + ext)
        entry.files["prev"] = os.path.join(workdir, prefix + "_prev.cpt")

    generate_runentry(run_entry)
    generate_warnentry(warn_entry)
    generate_finentry(fin_entry)
    generate_cnclentry(cncl_entry)
    generate_faltentry(falt_entry)
    generate_failentry(fail_entry)

    savemgr.save(savename)
    daemon(savename, live=False)

def mergefiles(mergename, savefiles):
    newsave = job_utils.SaveJobs()

    for sname in savefiles:
        readsave = job_utils.SaveJobs()
        readsave.load(sname)
        newsave.jobs.extend(readsave.jobs)
        newsave.attr.update(readsave.attr)
    newsave.attr["name"] = mergename

    svdir = os.path.dirname(savefiles[0])
    newsave.save(os.path.join(svdir, mergename + ".save"))

def extractinfo(savename):
    savemgr = job_utils.SaveJobs()
    savemgr.load(savename)


    for entry in savemgr.get_jobs():
        logscan = check_log(entry)
        if "Count" in logscan.log_entries:
            nsamples = logscan.log_numbers["Count"]
            print(entry.jobname)
            print("\t".join(nsamples))

def main(argv=None):  # IGNORE:C0111
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    # program_name = os.path.basename(sys.argv[0])
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
        parser = ArgumentParser(description=program_license)
#        parser.add_argument("-v", "--verbose", dest="verbose", action="count",
#            help="set verbosity level [default: %(default)s]")
        parser.add_argument('-V', '--version', action='version',
            version=program_version_message)
        parser.add_argument('--live', action="store_true", help="Use this to" +
            " enable actual daemon control. [default: %(default)s]" +
            " If not activated, the program will run in dry-mode, that is," +
            " no jobs will be submitted and no files created, merely output"
            " detailing what would happen in live-mode.")
        parser.add_argument('save', help="File generated by run of genscript," +
            " which indicates the files that were generated by that run." +
            " These jobs are the ones that will be monitored and analyzed.",
            nargs="+")
        parser.add_argument('-s', '--submit', help="Instead of starting" +
            " immediately, submit a script using the slurm scheduler. The" +
            " value passed will be used as the job name.", default=None)
        parser.add_argument('-t', '--time', help="Time to request on the" +
            " cluster, this will be the max time allowed for this to run." +
            " Format: d-HH:mm:ss", default=None)
        parser.add_argument('--test', action="store_true", help="Perform a" +
            " test run of the job daemon. [default: %(default)s]" +
            " This option will create a set of test files based on the value" +
            " of 'save' and run the daemon in dry-mode to monitor them." +
            " The output from the daemon can be used" +
            " to debug issues. This overrides all other options. The process" +
            " is run in the current workspace (does not submit an sbatch job).")
        parser.add_argument('-m', '--merge', help="merge a number of save" +
            " files together using the name passed to this param. New save" +
            " is created in the same location as first param", default=None)
        parser.add_argument('-e', '--extract', help="extract information" +
            " from the log files of a set of jobs detailed by the save file." +
            " Does not start the daemon.", action="store_true")

        # Process arguments
        args = parser.parse_args()
#         verbose = args.verbose
#
#         if verbose > 0:
#             print("Verbose mode on")

        save_name = None
        if args.save:
            save_name = args.save[0]
        submit = args.submit
        merge = args.merge
        time = "7-00:00:00"
        live = False
        test = False
        extract = False
        if args.time:
            time = args.time
        if args.live:
            live = True
        if args.test:
            test = True
        if args.extract:
            extract = True

        if save_name:
            if merge is not None:
                mergefiles(merge, args.save)
            elif extract:
                extractinfo(save_name)
            elif test:
                test_daemon(save_name)
            elif submit != None:
                reschedule_self(submit, save_name, time=time, live=live)
            else:
                daemon(save_name, live=live)
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
















