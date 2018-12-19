function progress = worker_progress_tracker(tracker_fn)
%txtappend functionality, plus a little extra so the workers don't step on
%each other's toes. 

f = fopen(tracker_fn, 'a+');
while f == -1 %if at first you don't succeed.. 
    f = fopen(tracker_fn, 'a+');
end
fprintf(f,'1\n'); %add a tally
frewind(f);
[~,progress] = fscanf(f,'%i'); %progress is count of ints read 
fclose(f);
