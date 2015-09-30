#!/usr/local/bin/python2.7
# encoding: utf-8
'''
relaunch -- modify .slurm files to run gromacs from checkpoint

@author:    Tyler P

@copyright: 2015 Shirts Group, University of Virginia. 
            All rights reserved.

@contact:   jtp4kc@virginia.edu
@deffield   updated: 30 Sep 2015
'''

import sys
import os

from argparse import ArgumentParser
# from argparse import RawDescriptionHelpFormatter

__all__ = []
__version__ = 0.1
__date__ = '2015-09-30'
__updated__ = '2015-09-30'

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

def update_slurm(filename, cptname=None, ignore=None, resume=None,
        time=None, copyfile=None):
    path = os.path.realpath(filename)
    new = path + ".bak"
    os.system("mv " + path + " " + new)

    timekey = "--time="
    ignore_active = False
    output = open(path, 'w')
    for line in open(new, 'r'):
        line = line.rstrip()
        if (resume != None) and (resume in line):
            ignore_active = False
            if copyfile != None:
                for copyline in open(copyfile, 'r'):
                    copyline = copyline.rstrip()
                    copyfile.write(copyline + "/n")
        if not ignore_active:
            if (timekey in line) and (time != None):
                split = line.split(timekey)
                line = split[0] + timekey + time
            if "grompp" in line:
                continue
            if ("mdrun" in line) and (cptname != None):
                line = line + " -cpi " + cptname
            output.write(line + "\n")
        if (ignore != None) and (ignore in line):
            ignore_active = True
    output.close()

def conv(string):
    if string == None:
        return None
    return str(string)

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
  Copyright 2015 organization_name. All rights reserved.

  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license)  # , formatter_class=RawDescriptionHelpFormatter)
#        parser.add_argument("-r", "--recursive", dest="recurse", action="store_true", help="recurse into subfolders [default: %(default)s]")
        parser.add_argument("-v", "--verbose", dest="verbose", action="count", help="set verbosity level [default: %(default)s]")
#        parser.add_argument("-i", "--include", dest="include", help="only include paths matching this regex pattern. Note: exclude is given preference over include. [default: %(default)s]", metavar="RE")
#        parser.add_argument("-e", "--exclude", dest="exclude", help="exclude paths matching this regex pattern. [default: %(default)s]", metavar="RE")
        parser.add_argument('-V', '--version', action='version', version=program_version_message)
        parser.add_argument('-n', "--name", help="name of cpt file to " +
            " reference when specifying a restart point")
        parser.add_argument("-i", "--ignore", help="ignore lines after" +
            " this str until the resume str is encountered")
        parser.add_argument('-r', '--resume', help="continue copying" +
            " lines after this string is encountered")
        parser.add_argument('-t', '--time', help="replace the time" +
            " request string in the file with this one")
        parser.add_argument('-c', '--copy', help="if an ignore section" +
            "is encountered, when it finishes, copy this file's" +
            " contents into the ignored space")
        parser.add_argument('--restore', action="store_true", help=
            "look for .bak files and move them to their original name")
        parser.add_argument(dest="paths", help="slurm script files to" +
            " modify [default: %(default)s]", metavar="path", nargs='*')

        # Process arguments
        args = parser.parse_args()

        paths = args.paths
        verbose = args.verbose
#         recurse = args.recurse
#         inpat = args.include
#         expat = args.exclude
        restore = args.restore
        name = conv(args.name)
        ignore = conv(args.ignore)
        resume = conv(args.resume)
        time = conv(args.time)
        copy = conv(args.copy)
        if not name.endswith(".cpt"):
            name += ".cpt"

        if verbose > 0:
            print("Verbose mode on")
#             if recurse:
#                 print("Recursive mode on")
#             else:
#                 print("Recursive mode off")

#         if inpat and expat and inpat == expat:
#             raise CLIError("include and exclude pattern are equal! Nothing will be processed.")

        for inpath in paths:
            if restore and inpath.endswith(".bak"):
                filename = os.path.realpath(inpath)
                new = filename.replace(".bak", "")
                os.system("mv " + filename + " " + new)
            if inpath.endswith(".slurm"):
                update_slurm(inpath, name, ignore, resume, time, copy)
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
        profile_filename = 'relaunch_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())
