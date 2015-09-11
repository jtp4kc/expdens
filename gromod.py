#!/usr/local/bin/python2.7
# encoding: utf-8
'''
gromod -- manipulate a .gro file

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
__date__ = '2015-09-10'
__updated__ = '2015-09-10'

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

    def plus(self, vec):
        add = Vector()
        add.x = self.x + vec.x
        add.y = self.y + vec.y
        add.z = self.z + vec.z
        return add

    def minus(self, vec):
        sub = Vector()
        sub.x = self.x - vec.x
        sub.y = self.y - vec.y
        sub.z = self.z - vec.z
        return sub

    def dot(self, vec):
        return self.x * vec.x + self.y * vec.y + self.z * vec.z

    def cross(self, vec):
        prod = Vector()
        prod.x = self.y * vec.z - self.z * vec.y
        prod.y = -self.x * vec.z + self.z * vec.x
        prod.z = self.x * vec.y - self.y * vec.x
        return prod

class Box():

    def __init__(self):
        self.v1 = Vector()
        self.v2 = Vector()
        self.v3 = Vector()

    def reset(self, value=0.0):
        self.v1.x = value
        self.v1.y = value
        self.v1.z = value
        self.v2.x = value
        self.v2.y = value
        self.v2.z = value
        self.v3.x = value
        self.v3.y = value
        self.v3.z = value

    def parse(self, line):
        self.reset()
        # 6.69337    6.69337    4.73316 <- cubic dimensions
        # v1(x)        v2(y)        v3(z)
        #
        # 0.00000    0.00000    0.00000
        # v1(y) Z    v1(z) Z    v2(x)
        #
        # 0.00000    3.34668    3.34668 <-triclinic slant modifiers
        # v2(z) Z    v3(x)        v3(y)
        info = line.strip().split()
        num = len(info)
        if num < 1:
            return
        self.v1.x = float(info[0])
        if num < 2:
            return
        self.v2.y = float(info[1])
        if num < 3:
            return
        self.v3.z = float(info[2])
        if num < 4:
            return
        self.v1.y = float(info[3])
        if num < 5:
            return
        self.v1.z = float(info[4])
        if num < 6:
            return
        self.v2.x = float(info[5])
        if num < 7:
            return
        self.v2.z = float(info[6])
        if num < 8:
            return
        self.v3.x = float(info[7])
        if num < 9:
            return
        self.v3.y = float(info[8])

    def output(self):
        spacing = "   "
        frmt = "{0:7.5f}"
        outlist = []
        outlist.append(frmt.format(self.v1.x))
        outlist.append(frmt.format(self.v2.y))
        outlist.append(frmt.format(self.v3.z))
        outlist.append(frmt.format(self.v1.y))
        outlist.append(frmt.format(self.v1.z))
        outlist.append(frmt.format(self.v2.x))
        outlist.append(frmt.format(self.v2.z))
        outlist.append(frmt.format(self.v3.x))
        outlist.append(frmt.format(self.v3.y))

        line = ""
        for entry in outlist:
            line += entry
            line += spacing

        return line

    def volume(self):
        # https://www.math.ucdavis.edu/~daddel/linear_algebra_appl/
        #    Applications/Determinant/Determinant/node11.html
        # Suppose three vectors v, u and w in three dimensional space
        # R^3 are given so that they do not lie in the same plane.
        # These three vectors form three edges of a parallelepiped.
        # The volume of this parallelepiped (is the product of area
        # of the base and altitude) is equal to the scalar triple
        # product u dot ( v cross w ).
        return self.v1.dot(self.v2.cross(self.v3))

    def convert(self, a, b, c, alpha, beta, gamma):
        #  box vectors     box vector angles
        #  a    b    c     <bc    <ac    <ab
        pass

class Atom():

    def __init__(self):
        self.loc = Vector()
        self.velocity = None
        self.resname = ""
        self.name = ""
        self.resid = 0
        self.atomn = 0
        self.position_precision = 3
        self.velocity_precision = 4

    def parse(self, line):
        # |aaaaabbbbbcccccddddd----|-------|-------|--- ...
        # |  163TMP     N1 2604   4.333   4.508   1.853
        num = len(line)
        if num >= 5:
            self.resid = int(line[0:5])
        if num >= 5 + 5:
            self.resname = line[5:10].strip()
        if num >= 10 + 5:
            self.name = line[10:15].strip()
        if num >= 15 + 5:
            self.atomn = int(line[15:20])
        if num > 20:  # numbers format is different
            x, y, z = (0.0,) * 3
            self.velocity = None
            numbers = line[20:]
            meas = numbers.split('.')
            if len(meas) > 2:
                pos_spacing = len(meas[1]) + 1
                self.position_precision = pos_spacing - 5
                if len(numbers) >= (1 * pos_spacing):
                    x = float(numbers[
                        (0 * pos_spacing):(1 * pos_spacing)])
                if len(numbers) >= (2 * pos_spacing):
                    y = float(numbers[
                        (1 * pos_spacing):(2 * pos_spacing)])
                if len(numbers) >= (3 * pos_spacing):
                    z = float(numbers[
                        (2 * pos_spacing):(3 * pos_spacing)])
                    numbers = numbers[(3 * pos_spacing):]
                    meas = numbers.split('.')
                    if len(meas) > 2:
                        v_spacing = len(meas[1]) + 1
                        vx, vy, vz = (0.0,) * 3
                        self.velocity = Vector()
                        self.velocity_precision = v_spacing - 4
                        if len(numbers) >= (1 * v_spacing):
                            vx = float(numbers[
                                (0 * v_spacing):(1 * v_spacing)])
                        if len(numbers) >= (2 * v_spacing):
                            vy = float(numbers[
                                (1 * v_spacing):(2 * v_spacing)])
                        if len(numbers) >= (3 * v_spacing):
                            vz = float(numbers[
                                (2 * v_spacing):(3 * v_spacing)])
                        self.velocity.x = vx
                        self.velocity.y = vy
                        self.velocity.z = vz
            self.loc.x = x
            self.loc.y = y
            self.loc.z = z

    def output(self):
        out = ""
        frmt = "{0:>5}{1:<5}{2:>5}{3:>5}"
        if self.loc:
            d = self.position_precision
            n = self.position_precision + 5
            frmt += ("{4:" + str(n) + ".0" + str(d) + "}")
            frmt += ("{5:" + str(n) + ".0" + str(d) + "}")
            frmt += ("{6:" + str(n) + ".0" + str(d) + "}")
            print(frmt)
            if self.velocity:
                d = self.velocity_precision
                n = self.velocity_precision + 4
                frmt += ("{7:" + str(n) + ".0" + str(d) + "}")
                frmt += ("{8:" + str(n) + ".0" + str(d) + "}")
                frmt += ("{9:" + str(n) + ".0" + str(d) + "}")
                out = frmt.format(self.resid, self.resname, self.name,
                    self.atomn, self.loc.x, self.loc.y, self.loc.z,
                    self.velocity.x, self.velocity.y, self.velocity.z)
            else:
                out = frmt.format(self.resid, self.resname, self.name,
                    self.atomn, self.loc.x, self.loc.y, self.loc.z)
            print(out)
        else:
            out = frmt.format(self.resid, self.resname, self.name,
                self.atomn)
        return out

class Gro:

    def __init__(self):
        self.title = ""
        self.atoms = []
        self.box = Box()

    def add_atom(self, line):
        atom = Atom()
        atom.parse(line)
        self.atoms.append(atom)

    def resize(self, image_distance):
        # rhombic dodecahedron
        self.box.reset()
        self.box.v1.x = image_distance
        self.box.v2.y = image_distance
        self.box.v3.x = 0.5 * image_distance
        self.box.v3.y = 0.5 * image_distance
        self.box.v3.z = math.sqrt(2) * 0.5 * image_distance

def read_gro(gro_filename):
    gro = Gro()
    atom_line = None
    count = 0
    for line in open(gro_filename, 'r'):
        count += 1
        if count == 1:
            gro.title = line.strip()
        elif count == 2:
            gro.atom_count = int(line.strip())
        else:
            if atom_line != None:  # there was a previous line
                gro.add_atom(atom_line.rstrip())  # left space important
                atom_line = None
            atom_line = line  # save the line to be parsed as box dims
    if atom_line != None:
        gro.box.parse(atom_line.strip())
    return gro

def write_gro(gro_filename, gro):
    output = open(gro_filename, 'w')
    output.write(gro.title + "\n")
    output.write(str(len(gro.atoms)) + "\n")
    for atom in gro.atoms:
        output.write(atom.output() + "\n")
    output.write(gro.box.output() + "\n")
    output.close()

def modify(args):
    verbose = args.verbose
    gro_filename = args.name
    out_filename = args.output
    isolate = args.isolate
    shrink = args.shrink
    move = args.move
    offset = [args.offsetx, args.offsety, args.offsetz]
    exclude = []
    cutoff = float(args.cutoff)
    if type(args.exclude) == list:
        exclude.extend(args.exclude)
    else:
        exclude.append(args.exclude)

    gro = read_gro(gro_filename)
    if verbose > 0:
        print("Read: " + str(gro.title))
        print("Found " + str(len(gro.atoms)) + " atom entries")
        print("Volume of gro file: " + str(gro.box.volume()))
    write_gro(out_filename, gro)

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
        parser.add_argument('-n', "--name", default="ligand.gro",
            help="name of gro file [default: %(default)s]")
        parser.add_argument('-o', "--output", default="output.gro",
            help="name of output gro file [default: %(default)s]")
        parser.add_argument('-i', "--isolate", default=None,
            help="create a file with only the first instance" +
            " of this residue name [default: %(default)s]")
        parser.add_argument('-m', "--move", default=None,
            help="move the residue found by --isolate to either" +
            " 'origin', 'center', or 'offset' [default: %(default)s]")
        parser.add_argument('-x', "--offsetx", default=0.0,
            help="offset for move [default: %(default)s]", type=float)
        parser.add_argument('-y', "--offsety", default=0.0,
            help="offset for move [default: %(default)s]", type=float)
        parser.add_argument('-z', "--offsetz", default=0.0,
            help="offset for move [default: %(default)s]", type=float)
        parser.add_argument('-s', "--shrink", action="store_true",
            help="shrink the gro box to fit the contents" +
            "using the general rule for periodic interactions" +
            " (2xCutoff + Width) [default: %(default)s]")
        parser.add_argument('-e', "--exclude", default=["SOL", "CL-"],
            help="ignore these residues unless requested elsewhere," +
            " can be number or name [default: %(default)s]",
            nargs="*")
        parser.add_argument('-c', "--cutoff", default=1.0,
            help="cutoff (in nm) to use when computing size changes" +
            " [default: %(default)s]", type=float)

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
