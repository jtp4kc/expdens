'''
Created on Jun 19, 2015

@author: Tyler P
'''
import sys, os

def place_vars(base, suffix, dir):
    text_ = """#!/bin/sh
        #SBATCH --job-name=%s%s
        """ % (base, suffix)
    text_ = text_ + """#SBATCH --partition=serial
        #SBATCH --nodes=1
        #SBATCH --ntasks=10
        #SBATCH --time=2-06:00:00
        #SBATCH --signal=15 --comment="15 = SIGTERM"
        """
    text_ = text_ + """#SBATCH --output=%s.out
        """ % (os.path.join(dir, base + suffix, base))
    text_ = text_ + """#SBATCH --workdir=%s
        """ % (os.path.join(dir, base + suffix))
    text_ = text_ + """#SBATCH --checkpoint=06:00:00
        #SBATCH --checkpoint-dir=/scratch/jtp4kc/checkpoints
        #SBATCH --mail-type=ALL
        #SBATCH --mail-user=jtp4kc@virginia.edu
        """
#    text_ = text_ + """#SBATCH --comment="#SBATCH --array-0,1,3" """
    text_ = text_ + """
        # submit a single SLURM job for a free energy calc
        module load gromacs/5.0.2-sse
        export THREADINFO="-nt 10 "
        export GMX_SUPPRESS_DUMP=1 #prevent step file output when parameters don't converge
        
        #export OMP_NUM_THREADS=10
        
        sleep 1
        
        #Change the job status to 'SUBMITTED'
        echo "SUBMITTED %s" >> ../jobstatus.txt
        """ % suffix
    text_ = text_ + """echo "Job Started at `date`"
        
        #NVT PRODUCTION
        """
    text_ = text_ + """INFOLDER="../%s-in/"
        GNAME="%s"
        """ % (base, base)
    text_ = text_ + """GNAMEIN="${INFOLDER}${GNAME}"
        
        cp ${GNAMEIN}%s-in.gro ${GNAME}-in.gro
        """ % (suffix)
    text_ = text_ + """cp ${GNAMEIN}.top ${GNAME}.top
        cp ${GNAMEIN}.ndx ${GNAME}.ndx
        #avoid overwriting file
        if [ ! -f ${GNAME}.mdp ];
        then
            cp ${GNAMEIN}.mdp ${GNAME}.mdp
        fi
        
        grompp_d -c ${GNAME}-in.gro -p ${GNAME}.top -n ${GNAME}.ndx -f ${GNAME}.mdp -o ${GNAME}.tpr -maxwarn 15
        mdrun_d ${THREADINFO} -deffnm ${GNAME}
        
        echo "FINISHED %s" >> ../jobstatus.txt
        """ % suffix
    text_ = text_ + """# print end time
        echo
        echo "Job Ended at `date`"
        echo "###################################################################"
        """
    text = ''
    for line in text_.splitlines():
        text = text + line.lstrip() + '\n'
    
    return text

def generate(base_name, submit):
    dir = '/scratch/jtp4kc/simulations/t41meth'
#    dir = os.getcwd()
    for i in range(20):
        suffix = '_{0:0>2}'.format(i)
        file_name = base_name + suffix + '.slurm'
        file_ = open(file_name, 'w')
        file_.write(place_vars(base_name, suffix, dir))
        file_.close()
        if submit:
            dir_name = os.path.join(dir, base_name + suffix)
            if not os.path.exists(dir_name):
                os.mkdir(dir_name)
            os.system("sbatch " + file_name)
            print("sbatch Job" + suffix)

def main(argv=None):
    if not argv:
        argv = sys.argv[1:]

#    parser = argparse.ArgumentParser(description=("""Generates Rivanna SLURM
#        scripts"""))
#    parser.add_argument("-v", "--verbose", help="""Increase output frequency
#        and detail. Stacks three times.""", action="count")
#     parser.add_argument("-bep", help="""Optional configuration file to specify
#         command line parameters and more.""", default=None,
#         metavar="params.bep")
#    parser.add_argument("-n", help="""Base name for output""", default="sim",
#        metavar="sim")
#    parser.add_argument("--submit", action='store_true')
#
#    args = parser.parse_args(argv)
    
    generate("t41meth", True)

if __name__ == '__main__':
    main(sys.argv[1:])
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
