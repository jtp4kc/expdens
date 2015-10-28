#!/usr/local/bin/python2.7
# encoding: utf-8
'''
topmod -- manipulate a .top file

@author:     jtp4kc, Tyler P

@copyright:  2015 Shirts Group, University of Virginia. All rights reserved.

@contact:    jtp4kc@virginia.edu
@deffield    updated: Updated
'''

import sys
import os
import math

from argparse import ArgumentParser
# from argparse import RawDescriptionHelpFormatters

__all__ = []
__version__ = 0.1
__date__ = '2015-27-10'
__updated__ = '2015-27-10'

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

class Vector():

    def __init__(self):
        self.x = 0
        self.y = 0
        self.z = 0

    def magnitude(self):
        return math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

    def scale(self, scalar, make_new_return=True):
        result = self
        if make_new_return:
            result = Vector()
        result.x = self.x * scalar
        result.y = self.y * scalar
        result.z = self.z * scalar
        return result

    def plus(self, vec, make_new_return=True):
        add = self
        if make_new_return:
            add = Vector()
        add.x = self.x + vec.x
        add.y = self.y + vec.y
        add.z = self.z + vec.z
        return add

    def minus(self, vec, make_new_return=True):
        sub = self
        if make_new_return:
            sub = Vector()
        sub.x = self.x - vec.x
        sub.y = self.y - vec.y
        sub.z = self.z - vec.z
        return sub

    def dot(self, vec):
        return self.x * vec.x + self.y * vec.y + self.z * vec.z

    def cross(self, vec, make_new_return=True):
        x = self.y * vec.z - self.z * vec.y
        y = -self.x * vec.z + self.z * vec.x
        z = self.x * vec.y - self.y * vec.x
        prod = self
        if make_new_return:
            prod = Vector()
        prod.x = x
        prod.y = y
        prod.z = z
        return prod

class TopItem():

    def __init__(self):
        self.data = {}

    def parse(self, line):
        info = line.split()
        cid = 0
        for datum in info:
            cid += 1
            self.data[cid] = datum

    def output(self):
        out = ""
        n = max(self.data.keys())
        for i in range(1, n + 1):
            out += self.data[i]
            out += " "
        return out

class Atom(TopItem):

    def __init__(self):
        TopItem.__init__(self)
        self.nr = 0
        self.col = 0

    def parse(self, line):
        TopItem.parse(self, line)
        #   nr       type  resnr residue  atom   cgnr    charge       mass
        #    1         ca      1    TMP     C1      1  -0.130100    12.0100
        n = max(self.data.keys())
        for i in range(1, n + 1):
            datum = self.data[i]
            if isinstance(datum, str) and datum.isdigit():
                self.nr = int(datum)
                self.col = i
                break

    def output(self):
        self.data[self.col] = str(self.nr)
        return TopItem.output(self)

class Bond(TopItem):

    def __init__(self):
        TopItem.__init__(self)
        self.ai = 0
        self.aj = 0
        self.atom1 = None
        self.atom2 = None

    def parse(self, line):
        TopItem.parse(self, line)
        #    ai     aj funct         c0         c1         c2         c3
        #     5      6     1   0.13870 400325.120000
        self.ai = int(self.data[1])
        self.aj = int(self.data[2])

    def output(self):
        self.data[1] = str(self.ai)
        self.data[2] = str(self.aj)
        return TopItem.output(self)

class Pair(TopItem):

    def __init__(self):
        TopItem.__init__(self)
        self.ai = 0
        self.aj = 0
        self.atom1 = None
        self.atom2 = None

    def parse(self, line):
        TopItem.parse(self, line)
        #    ai     aj funct         c0         c1         c2         c3
        #     3      4     1
        self.ai = int(self.data[1])
        self.aj = int(self.data[2])

    def output(self):
        self.data[1] = str(self.ai)
        self.data[2] = str(self.aj)
        return TopItem.output(self)

class Angle(TopItem):

    def __init__(self):
        TopItem.__init__(self)
        self.ai = 0
        self.aj = 0
        self.ak = 0
        self.atom1 = None
        self.atom2 = None
        self.atom3 = None

    def parse(self, line):
        TopItem.parse(self, line)
        #    ai     aj     ak funct         c0         c1         c2         c3
        #     4      6      5     1   119.97005 562.329600
        self.ai = int(self.data[1])
        self.aj = int(self.data[2])
        self.ak = int(self.data[3])

    def output(self):
        self.data[1] = str(self.ai)
        self.data[2] = str(self.aj)
        self.data[3] = str(self.ak)
        return TopItem.output(self)

class Dihedral(TopItem):

    def __init__(self):
        TopItem.__init__(self)
        self.ai = 0
        self.aj = 0
        self.ak = 0
        self.al = 0
        self.atom1 = None
        self.atom2 = None
        self.atom3 = None
        self.atom4 = None

    def parse(self, line):
        TopItem.parse(self, line)
        #    ai     aj     ak     al funct         c0         c1         c2         c3         c4         c5
        #     3      1      2      4     1  180.00008  15.16700  2
        self.ai = int(self.data[1])
        self.aj = int(self.data[2])
        self.ak = int(self.data[3])
        self.al = int(self.data[4])

    def output(self):
        self.data[1] = str(self.ai)
        self.data[2] = str(self.aj)
        self.data[3] = str(self.ak)
        self.data[4] = str(self.al)
        return TopItem.output(self)

class Top:

    def __init__(self):
        self.title = ""
        self.atoms = []
        self.bonds = []
        self.pairs = []
        self.angles = []
        self.dihedrals = []
        self.comments = {}

    def add_atom(self, line):
        atom = Atom()
        atom.parse(line)
        self.atoms.append(atom)

    def add_bond(self, line):
        bond = Bond()
        bond.parse(line)
        self.bonds.append(bond)

    def add_pair(self, line):
        pair = Pair()
        pair.parse(line)
        self.pairs.append(pair)

    def add_angle(self, line):
        angle = Angle()
        angle.parse(line)
        self.angles.append(angle)

    def add_dihedral(self, line):
        dihedral = Dihedral()
        dihedral.parse(line)
        self.dihedrals.append(dihedral)

    def renumber(self):
        count = 0
        for atom in self.atoms:
            count += 1
            atom.nr = count
        for item in self.bonds:
            item.ai = item.atom1.nr
            item.aj = item.atom2.nr
        for item in self.pairs:
            item.ai = item.atom1.nr
            item.aj = item.atom2.nr
        for item in self.angles:
            item.ai = item.atom1.nr
            item.aj = item.atom2.nr
            item.ak = item.atom3.nr
        for item in self.dihedrals:
            item.ai = item.atom1.nr
            item.aj = item.atom2.nr
            item.ak = item.atom3.nr
            item.al = item.atom4.nr

def read_top(top_filename):
    categories = ["atoms", "bonds", "pairs", "angles", "dihedrals"]
    top = Top()
    current = []
    cat = "pre"
    for line in open(top_filename, 'r'):
        line = line.strip()
        if line.startswith(";"):
            current.append(line)
        elif line.startswith("["):
            category = line.split()[1]
            if category != cat:
                top.comments[cat] = current
                current = []
                cat = category
        elif not (line.isspace() or line == ""):
            if cat == categories[0]:
                top.add_atom(line)
            elif cat == categories[1]:
                top.add_bond(line)
            elif cat == categories[2]:
                top.add_pair(line)
            elif cat == categories[3]:
                top.add_angle(line)
            elif cat == categories[4]:
                top.add_dihedral(line)
    if current:
        top.comments[cat] = current

    atoms = {}
    for atom in top.atoms:
        atoms[atom.nr] = atom
    for item in top.bonds:
        item.atom1 = atoms[item.ai]
        item.atom2 = atoms[item.aj]
    for item in top.pairs:
        item.atom1 = atoms[item.ai]
        item.atom2 = atoms[item.aj]
    for item in top.angles:
        item.atom1 = atoms[item.ai]
        item.atom2 = atoms[item.aj]
        item.atom3 = atoms[item.ak]
    for item in top.dihedrals:
        item.atom1 = atoms[item.ai]
        item.atom2 = atoms[item.aj]
        item.atom3 = atoms[item.ak]
        item.atom4 = atoms[item.al]

    return top

def write_top(top_filename, top):
    top.renumber()

    categories = ["atoms", "bonds", "pairs", "angles", "dihedrals"]
    output = open(top_filename, 'w')

    def _format(topitems, category):
        comments = []
        if category in top.comments:
            comments = top.comments[category]
        output.write("[ {0} ]\n".format(category))
        for cmnt in comments:
            output.write(cmnt + "\n")
        for item in topitems:
            output.write(item.output() + "\n")

    for entry in top.comments["pre"]:
        output.write(entry + "\n")
    _format(top.atoms, categories[0])
    _format(top.bonds, categories[1])
    _format(top.pairs, categories[2])
    _format(top.angles, categories[3])
    _format(top.dihedrals, categories[4])
    output.close()

def modify(args):
    verbose = args.verbose
    top_filename = args.name
    out_filename = args.output

    top = read_top(top_filename)
    if verbose > 0:
        print("Found " + str(len(top.atoms)) + " atom entries")

    if verbose > 0:
        print("Writing " + str(len(top.atoms)) + " atom entries")
        print("To file " + os.path.realpath(out_filename))
    write_top(out_filename, top)

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
        parser = ArgumentParser(description=program_license)
        parser.add_argument('-v', "--verbose", action="count",
            help="verbosity, stacks [default: %(default)s]")
        parser.add_argument('-n', "--name", default="topol.top",
            help="name of top file [default: %(default)s]")
        parser.add_argument('-o', "--output", default="newtop.top",
            help="name of output top file [default: %(default)s]")

        # Process arguments
        args = parser.parse_args()

        verbose = args.verbose

        if verbose > 0:
            print("Verbose mode on")

        modify(args)

        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0

if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-h")
        sys.argv.append("-v")
    if TESTRUN:
        import doctest
        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'solvate_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())
