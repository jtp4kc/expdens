'''
Created on Oct 11, 2015

@author: jtp4kc
'''

import os, sys
import xml.etree.ElementTree as ETree
import xml.dom.minidom as minidom
import json, tempfile
import dateutil.parser as dateparse
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

        if self.filename == None:
            output = tempfile.NamedTemporaryFile(mode="w+t", prefix="jobsave-",
                suffix=".save", dir=".", delete=False)
            self.filename = output.name
        else:
            output = open(self.filename, "w+t")
        ETree.ElementTree(elem).write(output, None)
        if pp:  # pretty-print
            output.seek(0)  # SOF
            # reparse in order to pretty-print the xml
            # note: due to security issues with minidom, it is not used to
            #    load xml, except that generated here
            parsedom = minidom.parse(output)
            output.seek(0)  # SOF
            output.truncate(0)  # clear file
            parsedom.writexml(output, addindent="    ", newl="\n")
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
        return [].extend(self.jobs)

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

    def _element(self):
        elem = ETree.Element("job")
        elem.attrib = self.attr
        elem.set("name", str(self.jobname))
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


def submit_job(filename_of_slurm, jobname, doprint=True):
    import subprocess
    file_ = tempfile.NamedTemporaryFile(mode="w+t", prefix='sbo', dir='.')
    subprocess.call(["sbatch", filename_of_slurm], stdout=file_)
    file_.seek(0)  # reset to be able to read
    line = file_.readline().replace("\n", "")
    file_.close()  # file should be deleted shortly hereafter
    num = line.split(" ")[-1]
    if doprint:
        print(line)
        print("Sbatch'd Job " + jobname + " as jobname #" + num)
    return num

def submit_slurm(slurm_obj, filename, doprint=True):
    """ Submit a slurm script that exists as a Slurm instance.
    @param slurm_obj: Slurm instance - the data to write to file
    @param filename: string - filename as which to save the slurm
    """
    import subprocess
    jobname = slurm_obj.job_name
    filename = str(filename)
    if not filename.endswith(".slurm"):
        filename += ".slurm"

    slurmout = open(filename, "w")
    for line in slurm_obj.compile():
        slurmout.write(line)
    slurmout.close()

    file_ = tempfile.NamedTemporaryFile(mode="w+t", prefix='sbo', dir='.')
    subprocess.call(["sbatch", filename], stdout=file_)
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











