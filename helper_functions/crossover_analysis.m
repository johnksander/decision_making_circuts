function coords = crossover_analysis(data,celltype,Xgroups,varargin)
%takes input data matrix traces x timepoints, celltype logical struct,
%cell array referring to struct fields

if isempty(varargin)
    N = 1; %regular analysis, just find intersection
else
    N = varargin{1}; %num straps
    if numel(varargin) > 1
        comp_mode = varargin{2}; %parallel or not
    else
        comp_mode = 'serial'; %default if not specified  
    end
end

Xgroups = strrep(Xgroups,'-','');
T = 1:numel(data(1,:)); %use for intersections() X coordinates
coords = NaN(N,2);
pool_inds = cellfun(@(x) celltype.(x),Xgroups,'UniformOutput',false);

if N == 1 %regular analysis
    traces = cellfun(@(x) mean(data(x,:)),pool_inds,'UniformOutput',false);
    [Cx,Cy] = intersections(T,traces{1},T,traces{2});
    coords(1,:) = [Cx(1),Cy(1)]; %only take first if there's multiple for whatever reason
else %bootstrap distribution
    traces = cellfun(@(x) data(x,:),pool_inds,'UniformOutput',false);
    traces = cellfun(@(x) bootstrp(N,@mean,x),traces,'UniformOutput',false);
    switch comp_mode
        case 'serial'
            for idx = 1:N
                [Cx,Cy] = intersections(T,traces{1}(idx,:),T,traces{2}(idx,:));
                if ~isempty(Cx) && ~isempty(Cy)
                    coords(idx,:) = [Cx(1),Cy(1)];
                end
            end
        case 'parallel'
            parfor idx = 1:N
                [Cx,Cy] = intersections(T,traces{1}(idx,:),T,traces{2}(idx,:));
                if ~isempty(Cx) && ~isempty(Cy)
                    coords(idx,:) = [Cx(1),Cy(1)];
                end
            end
    end
    
end
