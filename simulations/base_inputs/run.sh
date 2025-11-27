#----- This is the name of the job.We can use this name to qdel.This also appears in qstat
#$ -N springs

#----- This is the name of the output file(the stdout redirection)
#$ -o scr.$JOB_ID

#----- This is the name of the error file.Possibles errors are reported here
#$ -e err.$JOB_ID

#----- cwd means that the qsub will be run  in current working directory
#$ -cwd

#----- Shell to use to run qsub default is /bin/csh that's why use the -S flag
#$ -S /bin/bash

#----- Tell the SGE that this is an array job, with "tasks" to be numbered 1 to 10000
##$ -t 1-10000

#----- When a single command in the array job is sent to a compute node,
#----- its task number is stored in the variable SGE_TASK_ID,
#----- so we can use the value of that variable to get the results we want:
#----- ~/programs/program -i ~/data/input.$SGE_TASK_ID -o ~/results/output.$SGE_TASK_ID
#----- NOTICE Exactly the same script will be run for 1-100000.so directories must use the SGE_TASK_ID
#----- ALSO IMPORTANT if we want to execute the jobs sequentially(one after the other) we use the 
#----- -hold_jid flag ($-hold_jid) using the -N name of the job 
#----- example qsub -hold_jid test_papis new_job will not start new_job until test_papis has finishe until test_papis has finisheid

#----- IN CASE OF OPEN MP threaded run:The parallel enviroment smp will run the np process on the same node

#----- use np procs of the parallel enviroment orte.This parallel enviroment is defined in qmon.The rule of this pe is
#----- fill up. Which means that a node is taken up fully before moving to another one.pe orte = openmpi

#$ -pe smp 1
export OMP_NUM_THREADS=np
ulimit -s unlimited

/home/luiz/hGAST/hGAST
