#!/usr/local/bin/python2.7
# encoding: utf-8
'''
fixxvg -- shortdesc

fixxvg is a description

It defines classes_and_methods

@author:     user_name

@copyright:  2016 organization_name. All rights reserved.

@license:    license

@contact:    user_email
@deffield    updated: Updated
'''

import sys
import os

import numpy as np
import visualizer_functions as vf
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

__all__ = []
__version__ = 0.1
__date__ = '2016-03-21'
__updated__ = '2016-03-21'

DEBUG = 0
TESTRUN = 1
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

def fix_xvg(xvgfile):
    print('Correcting state information for ' + xvgfile)
    dhdl = vf.read_dhdl(xvgfile)
    length = len(dhdl.time)
    if dhdl.states is None:
        print("Unable to read file or no states found, skipping")
        return

    states = np.array(dhdl.states)

    thermo = -1  # what column does the dhdl info start upon
    deltacount = 0  # how many columns of dhdl info
    energy_info = np.array(dhdl.data)
    c = -1
    for col in dhdl.columns:
        # since there can be other items in the xvg file, try to locate
        # the dhdl info. assume it is sequential, by column.
        c += 1
        if "\\xD\\f{}H \\xl\\f{} " in col:
            if thermo < 0:
                thermo = c
            deltacount += 1

    # ## New: issue with changing state values in GROMACS
    # Need to infer the state by the column which has the value of zero
    #    when the state index is equal to zero, due to a bug.
    #    Use the first state available with the value of 0.00.
    #    If this value is not found exactly, use a tolerance.
    tol = 5e-5
    state_val = 0  # the state which might be having issues
    for frame in range(length):
        s1 = states[frame]
        ene = energy_info[thermo:(thermo + deltacount), frame]
        if s1 == state_val:
            ind = np.where(ene == 0.0)
            if ind[0].shape[0] == 0:  # not found
                ind = np.where(np.abs(ene) <= tol)
            if ind[0].shape[0] > 0:  # found
                s1 = ind[0][0]  # first identified index
                states[frame] = s1

    xvgfile2 = xvgfile.replace(".xvg", "_fix.xvg")
    vf.write_dhdl(dhdl, xvgfile2)

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

  Created by user_name on %s.
  Copyright 2016 organization_name. All rights reserved.

  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-r", "--recursive", dest="recurse", action="store_true", help="recurse into subfolders [default: %(default)s]")
        parser.add_argument("-v", "--verbose", dest="verbose", action="count", help="set verbosity level [default: %(default)s]")
        parser.add_argument("-i", "--include", dest="include", help="only include paths matching this regex pattern. Note: exclude is given preference over include. [default: %(default)s]", metavar="RE")
        parser.add_argument("-e", "--exclude", dest="exclude", help="exclude paths matching this regex pattern. [default: %(default)s]", metavar="RE")
        parser.add_argument('-V', '--version', action='version', version=program_version_message)
        parser.add_argument(dest="paths", help="paths to folder(s) with source file(s) [default: %(default)s]", metavar="path", nargs='+')

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
            raise CLIError("include and exclude pattern are equal! Nothing will be processed.")

        for inpath in paths:
            if inpath.endswith(".xvg"):
                fix_xvg(inpath)
        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception, e:
        if DEBUG or TESTRUN:
            raise(e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2

if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-h")
        sys.argv.append("-v")
        sys.argv.append("-r")
#     if TESTRUN:
#         import doctest
#         doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'fixxvg_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())
