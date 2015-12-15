'''
Created on Oct 11, 2015

@author: jtp4kc
'''

import os, sys
import xml.etree.ElementTree as ETree
import xml.dom.minidom as minidom
import json, tempfile
import dateutil.parser as dateparse
import random
import math
from datetime import datetime

ENC = json.JSONEncoder()
DEC = json.JSONDecoder()

class SaveJobs():

    def __init__(self):
        self.jobs = []
        self.attr = dict()
        self.filename = None

    def save(self, filename=None, pp=True):
        if filename != None:
            self.filename = filename

        elem = ETree.Element("save_jobs")
        elem.attrib = self.attr
        for job in self.jobs:
            elem.append(job._element())

        text = ETree.tostring(elem)
        if pp:  # pretty-print
            # reparse in order to pretty-print the xml
            # note: due to security issues with minidom, it is not used to
            #    load xml, except that generated here
            parsedom = minidom.parseString(text)
            text = parsedom.toprettyxml(indent="    ")

        if self.filename == None:
            output = tempfile.NamedTemporaryFile(mode="w", prefix="jobsave-",
                suffix=".save", dir=".", delete=False)
            self.filename = output.name
        else:
            output = open(self.filename, "w")
        for line in text:
            output.write(line)
        output.close()

    def load(self, filename=None):
        if filename != None:
            self.filename = filename
        tree = ETree.parse(self.filename)
        root = tree.getroot()
        if root.tag == "save_jobs":
            self._process(root)

    def _process(self, root):
        self.attr = root.attrib
        for child in root:
            self.add_job(SaveEntry(child))

    def add_job(self, save_entry):
        if isinstance(save_entry, SaveEntry):
            self.jobs.append(save_entry)
            return True
        else:
            return False

    def get_jobs(self):
        retval = []
        retval.extend(self.jobs)
        return retval

    def remove_job(self, job):
        if isinstance(job, SaveEntry):
            self.jobs.remove(job)
            return True
        else:
            return False

class SaveEntry():
    """ Use only strings as names and as keys when saving data to this structure,
    as it is serialized to xml which will not allow integers, floats, etc., as 
    attribute names or values. 
    Certain values are serialized as xml 'text' using JSON, which allows those 
    values to be lists, integers, and most other python data types.
    """

    def __init__(self, _child=None):
        self.jobname = ""
        self.files = dict()
        self.attr = dict()

        if _child != None:
            self._parse(_child)

    def _element(self, check=True):
        elem = ETree.Element("job")
        elem.attrib = self.attr
        elem.set("name", str(self.jobname))
        if check:
            for key in elem.attrib:
                elem.attrib[key] = str(elem.attrib[key])
        for key in self.files:
            e = ETree.Element("files")
            e.set("type", str(key))
            e.text = ENC.encode(self.files[key])
            elem.append(e)
        return elem

    def _parse(self, xml_entry):
        if xml_entry.tag == "job":
            self.jobname = xml_entry.get("name", "NO_NAME")
            self.attr = xml_entry.attrib
            for e in xml_entry:
                if e.tag == "files":
                    key = e.get("type", "")
                    self.files[key] = DEC.decode(e.text)

class LogScan:

    def __init__(self, path):
        self.filepath = path
        self.oldscan = []
        self.midscan = []
        self.newscan = []
        self.weights = None
        self.num_of_steps = None
        self.sim_time = None
        self.log_delta = None
        self.finish_detected = False
        self.cancel_detected = False
        self.failure_detected = False
        self.fail_statement = ""
        self.log_entries = []
        self.log_numbers = {}
        self.energy_items = []
        self.energy_data = []

    def scan(self):
        if not os.path.exists(self.filepath):
            raise Exception('Simulation file not found ' + self.filepath)
        import subprocess
        file_ = tempfile.NamedTemporaryFile(mode="w+t", prefix='lso', dir='.')
        subprocess.call(["tail", "-500", self.filepath], stdout=file_)

        capture_fail = False
        file_.seek(0)  # reset to be able to read
        for line in file_:
            if capture_fail:
                capture_fail = False
                self.fail_statement = line
            if 'Step' in line and 'Time' in line:
                self.oldscan = self.midscan
                self.midscan = self.newscan
                self.newscan = []
            if 'Received the TERM signal' in line:
                self.cancel_detected = True
                break
            if 'Fatal error:' in line:
                self.failure_detected = True
                capture_fail = True
            if 'Finished mdrun' in line:
                self.finish_detected = True
            self.newscan.append(line)
        file_.close()
        if (len(self.newscan) == 0) and (len(self.midscan) == 0):
            raise Exception('Tail command appears to have failed.')

        # assume a short scan indicates the file was cut off
        #     which could get fooled by dumping the state matrix, but it's only
        #     one frame off
        if len(self.midscan) > len(self.newscan):
            self.newscan = self.midscan
            self.midscan = self.oldscan
        self.oldscan = None  # Free up memory, I suppose

        count = 0
        capture_weights = False
        capture_energies = False
        label = False
        for line in self.newscan:
            count += 1
            if count == 2:
                splt = line.split()
                if splt:
                    try:
                        self.num_of_steps = int(splt[0])
                        self.sim_time = float(splt[1])
                    except ValueError as e:
                        tb = sys.exc_info()[2]
                        print(str(tb) + ":" + str(e))
            if line.isspace():
                capture_weights = False
                capture_energies = False
            if capture_weights:
                splt = line.split()
                for i in range(len(splt)):
                    col = self.log_entries[i]
                    self.log_numbers[col].append(splt[i])
                weight = float(splt[5])
                self.weights.append(weight)
            if capture_energies:
                n = len(line)
                for i in range(n / 15 + 1):
                    part = line[i:(i + 15)].strip()
                    if part == "":
                        continue
                    if label:
                        self.energy_items.append(part)
                    else:
                        self.energy_data.append(part)
                label = not label
            if 'N' in line and 'Count' in line and 'G(in kT)' in line:
                for col in line.split():
                    self.log_entries.append(col)
                    self.log_numbers[col] = []
                capture_weights = True
                self.weights = []
            if 'Energies' in line:
                capture_energies = True
                label = True

        if len(self.midscan) > 1:
            line = self.midscan[1]
            splt = line.split()
            try:
                num_steps_prev = int(splt[0])
                self.log_delta = self.num_of_steps - num_steps_prev
            except ValueError as e:
                tb = sys.exc_info()[2]
                print(str(tb) + ":" + str(e))

    def get_wanglandau_weights(self):
        return self.weights

    def get_step_number(self):
        return self.num_of_steps

    def get_time_value(self):
        return self.sim_time

    def get_log_delta(self):
        return self.log_delta

class OutputScan:

    def __init__(self, path):
        self.filepath = path
        self.print_shake = True
        self.warning_detected = False
        self.cancel_detected = False
        self.fault_detected = False
        self.warn_statement = ""
        self.fault_statement = ""
        self.cancel_statement = ""

    def scan(self):
        if not os.path.exists(self.filepath):
            raise Exception('Simulation file not found ' + self.filepath)

        shakecount = 0
        for line in open(self.filepath, "r"):
            if 'slurmstepd:' in line:
                self.cancel_detected = True
                self.cancel_statement = line.replace("slurmstepd:", "")
            if 'Shake did not converge' in line and self.print_shake:
                self.warning_detected = True
                shakecount += 1
                self.warn_statement = 'SHAKE experienced issues x{0}'.format(
                    shakecount)
            if ' fault' in line:
                self.fault_detected = True
                self.fault_statement = line

class MDPReader:

    def __init__(self, path):
        self.filepath = path

    def scan(self):
        if not os.path.exists(self.filepath):
            raise Exception('Input file not found ' + self.filepath)

        results = dict()
        for line in open(self.filepath, "r"):
            line = line.strip()
            if not line.startswith(";"):
                keyval = line.split("=")
                if len(keyval) > 1:
                    key = keyval[0].strip()
                    val = keyval[1].strip()
                    results[key] = val

        return results

class ExtractFrames:

    def __init__(self, xtcpath, tprpath=None, logpath=None, mdppath=None):
        if tprpath is None:
            tprpath = xtcpath.replace(".xtc", ".tpr")
        if logpath is None:
            logpath = xtcpath.replace(".xtc", ".log")
        if mdppath is None:
            mdppath = xtcpath.replace(".xtc", ".mdp")
        self.xtc = xtcpath
        self.tpr = tprpath
        self.log = logpath
        self.mdp = mdppath
        self.num_of_steps = None
        self.log_delta = None
        self.estimate = 0

    def estimate_available_frames(self):
        logscan = LogScan(self.log)
        logscan.scan()
        log_delta = logscan.get_log_delta()
        num_of_steps = logscan.get_step_number()
        sim_time = logscan.get_time_value()

        if num_of_steps is None:
            raise Exception("Number of steps could not be determined from log" +
                            " file")

        if log_delta is None:
            mdpread = MDPReader(self.mdp)
            options = mdpread.scan()
            if "nstxtcout" in options:
                log_delta = int(options["nstxtcout"])
            elif "nstxout-compressed" in options:
                log_delta = int(options["nstxout-compressed"])
            elif "nstxout_compressed" in options:
                log_delta = int(options["nstxout_compressed"])
            else:
                raise Exception("Number of steps per frame could not be" +
                                " determined from log or mdp file")

        self.num_of_steps = num_of_steps
        self.log_delta = log_delta
        self.time_per_step = sim_time / num_of_steps
        self.estimate = num_of_steps / log_delta + 1
        return self.estimate

    def get_gro_frames(self, n_frames=1, method=0, framenums=False):
        """ Extract frames in .gro representation and return filenames
        @param n_frames: (int) number of frames desired out of extraction, >= 1
        @param method: (int) one of four methods:
                        0: (default) segment and return middle frames
                        1: equally spaced by segment, always includes last frame
                        2: random within a segment
                        3: equally spaced by segment, always include first frame
                        4: on segment boundaries, includes first and last frames
        Example: ("-" represents frames, "'" signifies segment boundary, 
                  "X" marks a selected frame, n_frames = 3)
            simulation line 0---------------------------------------------47
            method 0 (cntr) '      X       '       X       '       X      '
            method 1 (rght) '              X               X              X
            method 2 (rand) '          X   '    X          '         X    '
            method 3 (left) X              X               X              '
            method 4 (even) X                      X                      X
        @param framenums: (bool) If only interested in where the frames would
                be chosen, set to True - returns a list of frame numbers
        """
        if n_frames < 1:
            raise ValueError('Cannot extract less than one frame.')
        if self.num_of_steps is None or self.log_delta is None:
            self.estimate_available_frames()

        seg_locator = 0.5  # multiplier within segment, default middle
        if method == 1:
            seg_locator = 1  # right hand side
        elif method == 2:
            seg_locator = -1  # random
        elif method == 3:
            seg_locator = 0  # left hand side
        elif method == 4:
            seg_locator = -2  #  special, then left hand side

        nseg = n_frames
        if seg_locator == -2:
            nseg -= 1
        n_est = self.estimate - 1

        if nseg > 0:
            len_seg = n_est / float(nseg)  # length of a segment
            offset = n_est - int(len_seg * nseg)  # extra frames
        else:
            len_seg = n_est
            offset = 0
        if (method == 3):
            offset = 0  # ensure first frame is chosen

        frames = []
        num_left = n_frames
        if n_est <= n_frames:
            frames = range(n_est + 1)
        else:
            if seg_locator == -2:
                seg_locator = 1
                frames.append(0)
                num_left -= 1

            rand = random.Random()
            for n in range(num_left):
                if seg_locator < 0:
                    place = rand.randint(0, int(len_seg) - 1)
                else:
                    place = int(len_seg * seg_locator)
                frames.append(offset + int(n * len_seg) + place)
            if ((method == 1) or (method == 4)) and (len(frames) > 0):
                frames[-1] = n_est

        if framenums:
            return frames  # exit ############################################

        if (not os.path.exists(self.xtc)) or (not os.path.exists(self.tpr)):
            raise Exception('Simulation files not found ' + self.xtc +
                ' and/or ' + self.tpr)
        cmnd = ('echo "System" | trjconv -f {xtc} -s {tpr} -o {gro}' +
            ' -b {timeb} -e {timee} &>trjconv{index}.log')

        if ('GMXBIN' in os.environ and os.environ['GMXBIN'] and
            os.environ['GMXBIN'] not in os.environ['PATH']):
            print('Adding GMXBIN to PATH...')
            os.environ['PATH'] = os.environ['PATH'] + ':' + os.environ['GMXBIN']

        gro_files = []
        time_delta = 1
        if self.time_per_step is not None:
            time_delta = self.time_per_step
        ndigits = int(math.log(len(frames), 10)) + 1  # for printing
        for (i, n) in zip(frames, range(len(frames))):
            suffix = ("-{0:0>" + str(ndigits) + "}").format(n)
            modif = 0.5 * time_delta
            # frame * step/frame * time/step = time
            timeval = (i * self.log_delta) * time_delta
            tb = timeval - modif
            te = timeval + modif
            print("log delta")
            print(self.log_delta)
            print("time delta")
            print(time_delta)
            print("i")
            print(i)
            print("n")
            print(n)
            print("frames")
            print(frames)
            gro_name = self.xtc.replace(".xtc", suffix + '-in.gro')
            fmt = cmnd.format(xtc=self.xtc, tpr=self.tpr, gro=gro_name,
                timeb=tb, timee=te, index=suffix)
            os.system(fmt)
            if not os.path.exists(gro_name):
                raise Exception('Extraction of starting frames' +
                    ' unsuccessful - Err: ' + gro_name)
            gro_files.append(gro_name)

        return gro_files

def submit_job(filename_of_slurm, jobname, doprint=True):
    import subprocess
    here = os.getcwd()
    fname = os.path.basename(filename_of_slurm)
    slurmdir = os.path.dirname(filename_of_slurm)
    if slurmdir:
        os.chdir(slurmdir)

    file_ = tempfile.NamedTemporaryFile(mode="w+t", prefix='sbo', dir='.')
    subprocess.call(["sbatch", fname], stdout=file_)
    file_.seek(0)  # reset to be able to read
    line = file_.readline().replace("\n", "")
    file_.close()  # file should be deleted shortly hereafter
    num = line.split(" ")[-1]
    if doprint:
        if num.isdigit():
            print("sbatch'd job " + jobname + " as jobid #" + num)
        else:
            print(line)

    os.chdir(here)
    return num

def submit_slurm(slurm_obj, filename, doprint=True):
    """ Submit a slurm script that exists as a Slurm instance.
    @param slurm_obj: Slurm instance - the data to write to file
    @param filename: string - filename as which to save the slurm
    """
    jobname = slurm_obj.job_name
    filename = str(filename)
    if not filename.endswith(".slurm"):
        filename += ".slurm"

    slurmout = open(filename, "w")
    for line in slurm_obj.compile():
        slurmout.write(line)
    slurmout.close()

    return submit_job(filename, jobname, doprint)

class SerialDate:
    SEP = ":"
    TAG = "date"
    QUO = '"'

    @staticmethod
    def serialize(date_like):
        """ Converts datetime like object to string format
        Note: datetime is not time-zone aware by default
        """
        return (SerialDate.TAG + SerialDate.SEP + SerialDate.QUO +
                date_like.isoformat() + SerialDate.QUO)

    @staticmethod
    def deserialize(string_repr):
        """ Attempt to read a date from the given string.
        Returns None if read is not possible.
        """
        retval = None
        if string_repr.startswith(SerialDate.TAG + SerialDate.SEP):
            new_repr = string_repr.replace(SerialDate.TAG + SerialDate.SEP, "")
            iso_stamp = new_repr.strip(SerialDate.QUO)
            retval = dateparse.parse(iso_stamp)

        return retval

if __name__ == "__main__":
    # test
    print("SAVE")
    saver = SaveJobs()
    job1 = SaveEntry()
    job2 = SaveEntry()

    saver.attr["name"] = "test-save"
    job1.jobname = "1meth"
    job1.files["tpr"] = "1meth.tpr"
    job1.files["m>395@#$%@#%&!@)(<>?/"] = "junk"
    job1.files[47] = "junk"
    job1.attr["time"] = SerialDate.serialize(datetime.now())
    job2.jobname = "benz"
    job2.files["tpr"] = "benz.tpr"
    job2.files["top"] = ["benz.top"]
    job2.files["mdp"] = ["benz.mdp", "benz2.mdp"]

    saver.add_job(job1)
    saver.add_job(job2)

    saver.save()
    fname = saver.filename
    for line in open(fname, "r"):
        print(line.rstrip())

    print("")
    print("LOAD")
    saver = SaveJobs()
    saver.load(fname)
    print("fnam: " + saver.filename)
    print("attr: " + str(saver.attr))
    for job in saver.jobs:
        print("  jobn: " + job.jobname)
        print("    file: " + str(job.files))
        print("    attr: " + str(job.attr))












