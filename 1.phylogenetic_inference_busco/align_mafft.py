import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--path")
parser.add_argument("--pathout")
parser.add_argument("--scripts")
parser.add_argument("--batchsize")
args = parser.parse_args()

commands={}
count=0
for myfile in os.listdir(args.path):
    count+=1
    batch=count//int(args.batchsize)
    if batch not in commands:
        commands[batch]=[]
    commands[batch].append('mafft '+args.path+myfile+' > '+args.pathout+myfile+'\n')

for batch in commands:
    with open(args.scripts+'mafft_batch'+str(batch)+'.sh','w') as SH:
        SH.write('#!/bin/bash\n')
        SH.write('#SBATCH -J YoshimiBattlesThePinkRobots\n')
        SH.write('#SBATCH -o '+args.scripts+'mafft_batch'+str(batch)+'.out\n')
        SH.write('#SBATCH -e '+args.scripts+'mafft_batch'+str(batch)+'.err\n')
        SH.write('module load mafft/7.525\n')
        SH.write('\n'.join(commands[batch])+'\n')
    os.system('sbatch '+args.scripts+'mafft_batch'+str(batch)+'.sh')