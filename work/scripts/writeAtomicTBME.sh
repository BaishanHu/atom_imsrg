#!/bin/bash

# decide if we should run in the current shell or submit to the cluster
# interactive=true means run in the current shell.
interactive=false
NTHREADS=12
WORKDIR=$PWD

exe=../compiled/writeAtomicTBME

#vnn=/itch/exch/me2j/chi2b_srg0625_eMax14_lMax10_hwHO024.me2j.gz
#vnn=none
v3n=none

systemtype=atomic
# 1.0 Hartree ~= 27.21138505(60) eV (according to Wikipedia)
hw=27.21138505 
valence_space=He2
reference=He2
smax=200
emax=8
systemBasis=hydrogen
Lmax=4
e3max=0
lmax3=2
method=magnus
basis=HF
omega_norm_max=0.25 
#file3='file3e1max=12 file3e2max=28 file3e3max=12'

jobname="reference_${reference}_basis_${systemBasis}_emax_${emax}_lmax_${Lmax}_ME2J"

#Operators=KineticEnergy,InverseR,CorrE2b
#Operators=KineticEnergy,InverseR,ElectronTwoBody
Operators=
#scratch=

all_the_flags="emax=${emax} e3max=${e3max} method=${method} valence_space=${valence_space} hw=${hw} smax=${smax} omega_norm_max=${omega_norm_max} reference=${reference} Operators=${Operators} Lmax=${Lmax} systemtype=${systemtype} use_brueckner_bch=false systemBasis=${systemBasis}"
#command="valgrind ${exe} ${all_the_flags}"
#command="valgrind --trace-children=yes --leak-check=full --track-origins=yes -v ${exe} ${all_the_flags}"
# --read-inline-info=yes # This looks like it would be awesome.
command="${exe} ${all_the_flags}"
echo command is $command

if [ $interactive = 'true' ]; then
cd $WORKDIR 
export OMP_NUM_THREADS=${NTHREADS}
$command
else
echo jobname is ${jobname}
qsub -N ${jobname} -q batchmpi -d $PWD -l walltime=192:00:00 -l nodes=1:ppn=${NTHREADS} -l vmem=60gb -m ae -M davidedwardlivermore@gmail.com -j oe -o pbslog/${jobname}.o.`date +"%g%m%d%H%M"` << END
cd $WORKDIR
export OMP_NUM_THREADS=${NTHREADS}
$command
END
fi


# set the command with all the flags, commented out as I already have
#command=${exe} ${all_the_flags}

