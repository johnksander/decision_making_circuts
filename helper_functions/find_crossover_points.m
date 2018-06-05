function Cpoints = find_crossover_points(data)
%takes input data matrix traces x timepoints, finds where traces cross over
%useful for plotting & network analysis


%find each unique intersection point between two traces. For each intersection
%0) record the x,y coordinate
%1) record what traces intersect
%2) record which trace is increasing, which is decreasing
%for example:
%X=24, Y=10, I-stay, 'down', E-switch, 'up'


Ntraces = numel(data(:,1));
T = 1:numel(data(1,:)); %use for intersections() X coordinates
Cpoints = {}; %use legend labels as ref if you get confused
for i = 1:Ntraces
    trace_data = data(i,:);
    for j = 1:Ntraces
        if j > i %don't repeat comparisons
            Jtrace = data(j,:);
            Cx = intersections(T,trace_data,T,Jtrace);
            Cx = round(Cx); %crossover time
            Cx(Cx == max(T)) = max(T)-1; %don't let it be last timepoint for heading calc
            Cy = trace_data(Cx)'; %Y coordinate
            num_ints = numel(Cx);
            %find directionality by just comapring Cx+1
            heading = sign(trace_data(Cx+1) - Jtrace(Cx+1))';
            %trace index (duplicate to match n intersections) & trace direction
            i_info = [num2cell(repmat(i,num_ints,1)),label_direction(heading)];
            % j index & opposite direction 
            j_info = [num2cell(repmat(j,num_ints,1)),label_direction(-heading)];
            coords = num2cell([Cx,Cy]);
            Cpoints = cat(1,Cpoints,[coords,i_info,j_info]);
        end
    end
end

function h = label_direction(h)
h = num2cell(h);
h = cellfun(@num2str,h,'UniformOutput',false);
h = cellfun(@(x) strrep(x,'-1','down'),h,'UniformOutput',false);
h = cellfun(@(x) strrep(x,'1','up'),h,'UniformOutput',false);


%----old way I did it-----
%for each trace, find where it intersects other traces. For each
%intersection point, record whether the trace is increasing/decreasing
%relative to the intersecting trace.
%for example:
%I-stay, X=24, 'down'

%how I did it before

% Ntraces = numel(data(:,1));
% T = 1:numel(data(1,:)); %use for intersections() X coordinates
%
% Cpoints = cell(Ntraces,1); %use legend labels as ref if you get confused
% for idx = 1:Ntraces
%     curr_trace = idx == 1:Ntraces;
%     trace_data = data(curr_trace,:);
%     other_traces = data(~curr_trace,:);
%
%     for refidx = 1:Ntraces-1
%         Jtrace = other_traces(refidx,:);
%         Cx = intersections(T,trace_data,T,Jtrace);
%         Cx = round(Cx); %crossover time
%         %find directionality by just comapring Cx+1
%         heading = sign(trace_data(Cx+1) - Jtrace(Cx+1));
%         heading = num2cell(heading)';
%         heading = cellfun(@num2str,heading,'UniformOutput',false);
%         heading = cellfun(@(x) strrep(x,'-1','down'),heading,'UniformOutput',false);
%         heading = cellfun(@(x) strrep(x,'1','up'),heading,'UniformOutput',false);
%         Cx = num2cell(Cx);
%         Cy = cellfun(@(x) trace_data(x),Cx,'UniformOutput',false); %this also forces col vec
%         Cxyh = [Cx,Cy,heading];
%         Cxyh = num2cell(Cxyh,2);
%         Cpoints{idx} = cat(1,Cpoints{idx},Cxyh);
%     end
% end
%
