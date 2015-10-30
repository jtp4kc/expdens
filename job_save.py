'''
Created on Oct 11, 2015

@author: jtp4kc
'''

import os, sys
import xml.etree.ElementTree as ETree
import json, tempfile
import dateutil.parser as dateparse

ENC = json.JSONEncoder()
DEC = json.JSONDecoder()

class SaveJobs():

    def __init__(self):
        self.jobs = []
        self.attr = dict()

    def save(self, filename=None):
        elem = ETree.Element()
        elem.attrib = self.attr
        elem.tag = "save_jobs"
        for job in self.jobs:
            elem.append(job._element())

    def load(self, filename):
        tree = ETree.parse(filename)
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
        return [].extend(self.jobs)

    def remove_job(self, job):
        if isinstance(job, SaveEntry):
            self.jobs.remove(job)
            return True
        else:
            return False

class SaveEntry():

    def __init__(self, _child=None):
        self.jobname = ""
        self.files = dict()
        self.attr = dict()

        if _child:
            self._parse(_child)

    def _element(self):
        elem = ETree.Element()
        elem.tag = "job"
        elem.attrib = self.attr
        elem.set("name", self.jobname)
        for key in self.files:
            e = ETree.Element()
            e.tag = str(key)
            e.text = ENC.encode(self.files[key])
            elem.append(e)
        return elem

    def _parse(self, xml_entry):
        if xml_entry.tag == "job":
            self.jobname = xml_entry.get("name", "NO_NAME")
            self.attr = xml_entry.attrib
            for e in xml_entry:
                self.files[e.tag] = DEC.decode(e.text)

class LogScan:

    def __init__(self, path):
        self.filepath = path
        self.oldscan = []
        self.newscan = []
        self.weights = None
        self.num_of_steps = None
        self.finish_detected = False
        self.cancel_detected = False
        self.failure_detected = False
        self.fail_statement = ""
        self.log_entries = []
        self.log_numbers = {}

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
            if 'Step' in line and 'Time' in line and 'Lambda' in line:
                self.oldscan = self.newscan
                self.newscan = []
            if 'Received the TERM signal' in line:
                self.cancel_detected = True
            if 'Fatal error:' in line:
                self.failure_detected = True
                capture_fail = True
            if 'Finished mdrun' in line:
                self.finish_detected = True
            self.newscan.append(line)
        file_.close()
        if not (self.newscan or self.oldscan):
            raise Exception('Tail command appears to have failed.')

        # assume a short scan indicates the file was cut off
        #     which could get fooled by dumping the state matrix, but it's only
        #     one frame off
        if len(self.oldscan) > len(self.newscan):
            self.newscan = self.oldscan

        count = 0
        capture_weights = False
        for line in self.newscan:
            count += 1
            if count == 2:
                splt = line.split()
                if splt:
                    try:
                        self.num_of_steps = int(splt[0])
                    except Exception as e:
                        tb = sys.exc_info()[2]
                        print(tb + ":" + str(e))
            if line.isspace():
                capture_weights = False
            if capture_weights:
                splt = line.split()
                for i in range(len(splt)):
                    col = self.log_entries[i]
                    self.log_numbers[col].append(splt[i])
                weight = float(splt[5])
                self.weights.append(weight)
            if 'N' in line and 'Count' in line and 'G(in kT)' in line:
                for col in line.split():
                    self.log_entries.append(col)
                    self.log_numbers[col] = []
                capture_weights = True
                self.weights = []

    def get_wanglandau_weights(self):
        return self.weights

    def get_step_number(self):
        return self.num_of_steps

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

def submit_slurm(slurm_path, jobname, doprint=True):
    import subprocess
    file_ = tempfile.NamedTemporaryFile(mode="w+t", prefix='sbo', dir='.')
    subprocess.call(["sbatch", slurm_path], stdout=file_)
    file_.seek(0)  # reset to be able to read
    line = file_.readline().replace("\n", "")
    file_.close()  # file should be deleted shortly hereafter
    num = line.split(" ")[-1]
    if doprint:
        print(line)
        print("Sbatch'd Job " + jobname + " as jobname #" + num)
    return num

class SerialDate:
    SEP = ":"
    TAG = "date"
    QUO = '"'

    @staticmethod
    def serialize(date_like):
        """ Converts datetime like object to string format
        """
        univ = date_like.utcnow()
        return (SerialDate.TAG + SerialDate.SEP + SerialDate.QUO +
                univ.isoformat() + SerialDate.QUO)

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














