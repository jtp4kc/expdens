'''
Created on Jun 19, 2015

@author: Tyler P
'''
import sys, os

def place_vars(fields):
    fields['outputfile'] = os.path.join(fields['path'], fields['base'] +
        fields['suffix'], fields['base'])
    fields['workdir'] = os.path.join(fields['path'], fields['base'] +
        fields['suffix'])
    fields['infolder'] = os.path.join('..', fields['base'] + '-in', '')
    fields['jobstatus'] = os.path.join('..', 'jobstatus.txt')
    
    text_ = """#!/bin/sh
#SBATCH --job-name={base}{suffix}
#SBATCH --partition=serial
#SBATCH --nodes=1
#SBATCH --ntasks={ntasks}
#SBATCH --time=2-06:00:00
#SBATCH --signal=15 --comment="15 = SIGTERM"
#SBATCH --output={outputfile}.out
#SBATCH --workdir={workdir}
#SBATCH --checkpoint=06:00:00
#SBATCH --checkpoint-dir=/scratch/jtp4kc/checkpoints
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtp4kc@virginia.edu
#SBATCH --comment="#SBATCH --array-0,1,3"

# submit a single SLURM job for a free energy calc
module load gromacs/5.0.2-sse
export THREADINFO="-nt {ntasks} "
export GMX_SUPPRESS_DUMP=1 #prevent step file output when parameters don't converge

#export OMP_NUM_THREADS={ntasks}

sleep 1

#Change the job status to 'SUBMITTED'
echo "SUBMITTED {suffix}" >> {jobstatus}
echo "Job Started at `date`"

#NVT PRODUCTION
cp {infolder}{gro-in}-in.gro {base}-in.gro
cp {infolder}{base}.top {base}.top
cp {infolder}{base}.ndx {base}.ndx
#avoid overwriting file
if [ ! -f {base}.mdp ];
then
    cp {infolder}{mdp-in}.mdp {base}.mdp
fi

grompp_d -c {base}-in.gro -p {base}.top -n {base}.ndx -f {base}.mdp -o {base}.tpr -maxwarn 15
mdrun_d ${{THREADINFO}} -deffnm {base}

echo "FINISHED {suffix}" >> {jobstatus}
# print end time
echo
echo "Job Ended at `date`"
echo "###################################################################"
""".format(**fields)
    
    return text_

def generate(base_name, submit, randseed):
#    path = '/scratch/jtp4kc/simulations/t41meth'
    path = os.getcwd()
    for i in range(20):
        suffix = '_{0:0>2}'.format(i)
        fields = dict()
        fields['base'] = base_name
        fields['suffix'] = suffix
        fields['path'] = path
        fields['ntasks'] = 10
        fields['gro-in'] = base_name + suffix
        fields['mdp-in'] = base_name
        if randseed:
            fields['gro-in'] = base_name
            fields['mdp-in'] = base_name + suffix
            try: 
                mdp_out = open(os.path.join(path, base_name + "-in", 
                    base_name + suffix + ".mdp"), 'w')
                for line in open(os.path.join(path, base_name + "-in", 
                    base_name + ".mdp")):
                    
                    if 'seed' in line:
                        halves = line.split("=")
                        num = int(halves[1].strip())
                        num += i
                        line = halves[0] + "= " + str(num) + "\n"
                    mdp_out.write(line)
                mdp_out.close()
            except IOError as ioex:
                print(ioex)
        
        file_name = base_name + suffix + '.slurm'
        file_ = open(file_name, 'w')
        file_.write(place_vars(fields))
        file_.close()
        if submit:
            dir_name = os.path.join(path, base_name + suffix)
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
    
    generate("t41meth", False, False)

if __name__ == '__main__':
    main(sys.argv[1:])
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
