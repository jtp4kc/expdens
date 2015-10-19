'''
Created on Oct 11, 2015

@author: jtp4kc
'''

from param_versions.version_2_0 import Parameters
from param_versions.version_2_0 import Keys

class SaveKeys(Keys):

    def __init__(self):
        Keys.__init__(self)  # since this overrides, make sure to execute super
        self.name = "_save-name"
        self.jobs = "sim-job-tuples"
        self.files = "sim-files"
        self.folders = "sim-folder-tuples"

        self.add_keys("", self.jobs, self.files, self.folders)

save_keys = SaveKeys()

def save_defaults():
    options = dict()

    options[save_keys.jobs] = []
    options[save_keys.files] = []
    options[save_keys.folders] = []

    return options

class SaveJobs(Parameters):

    def __init__(self):
        Parameters.__init__(self, save_keys)
        self.option_defaults = save_defaults
