function Cpoints = find_crossover_points(data)
%takes input data matrix traces x timepoints, finds where traces cross over
%useful for plotting & network analysis 



%for each trace, find where it intersects other traces. For each
%intersection point, record whether the trace is increasing/decreasing
%relative to the intersecting trace. 
%for example:
%I-stay, X=24, 'down' 


Ntraces = numel(data(:,1));
T = 1:numel(data(1,:)); %use for intersections() X coordinates

Cpoints = cell(Ntraces,1); %use legend labels as ref if you get confused
for idx = 1:Ntraces
    curr_trace = idx == 1:Ntraces;
    trace_data = data(curr_trace,:);
    other_traces = data(~curr_trace,:);
    
    for refidx = 1:Ntraces-1
        Jtrace = other_traces(refidx,:);
        Cx = intersections(T,trace_data,T,Jtrace);
        Cx = round(Cx); %crossover time
        %find directionality by just comapring Cx+1
        heading = sign(trace_data(Cx+1) - Jtrace(Cx+1));
        heading = num2cell(heading)';
        heading = cellfun(@num2str,heading,'UniformOutput',false);
        heading = cellfun(@(x) strrep(x,'-1','down'),heading,'UniformOutput',false);
        heading = cellfun(@(x) strrep(x,'1','up'),heading,'UniformOutput',false);
        Cx = num2cell(Cx);
        Cy = cellfun(@(x) trace_data(x),Cx,'UniformOutput',false); %this also forces col vec
        Cxyh = [Cx,Cy,heading];
        Cxyh = num2cell(Cxyh,2);
        Cpoints{idx} = cat(1,Cpoints{idx},Cxyh);
    end
end

