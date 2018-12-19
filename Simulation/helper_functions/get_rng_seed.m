function myseed = get_rng_seed()

%courtesy francesco

mypid=feature('GetPid');

myclock=sum(clock);

myjobid=getenv('SLURM_JOBID');
myjobid=0.005*str2num(myjobid);

mytaskid=getenv('SLURM_ARRAY_TASK_ID');
mytaskid=10*str2num(mytaskid);

myseed= mypid + myclock + myjobid + mytaskid;

