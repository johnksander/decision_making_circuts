function myseed = get_rng_seed()

%courtesy francesco

mypid=feature('GetPid');

myclock=sum(clock);

myjobid=getenv('JOB_ID');
myjobid=0.005*str2num(myjobid);

mytaskid=getenv('SGE_TASK_ID');
mytaskid=10*str2num(mytaskid);

myseed= mypid + myclock + myjobid + mytaskid;

