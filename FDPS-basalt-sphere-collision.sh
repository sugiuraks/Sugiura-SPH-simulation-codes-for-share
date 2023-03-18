#PBS -N ANEOS-q=1_50
#PBS -l nodes=1
#PBS -q bulk-a
#PBS -j oe
echo Working Directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=20
aprun -n 2 -N 2 -d ${OMP_NUM_THREADS} -cc depth ./sph.out FDPS-bsc-Rt=50km-q=1_50-theta=0-v=5000ms-Ntot=55000-relaxed-basaltANEOS.bin FDPS-bsc-Rt=50km-q=1_50-theta=0-v=5000ms-Ntot=55000-relaxed-basaltANEOS > tmp${PBS_JOBID%.*}.txt
