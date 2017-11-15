function dataIDs = make_dataIDs(event_data,sim_names,cell_labels)

num_events = numel(event_data);
%num_models = numel(model_names);
num_celltypes = numel(cell_labels);

full_modelIDs = cellfun(@(x) strsplit(x,'_'), sim_names,'UniformOutput',false);
stimIDs =  cellfun(@(x) x{2}, full_modelIDs,'UniformOutput',false);
modelIDs =  cellfun(@(x) x{1}, full_modelIDs,'UniformOutput',false);
cellIDs = repmat({cell_labels'},numel(full_modelIDs),1);

dataIDs = cell(size(event_data)); %make labeling easier and more explicit
for idx = 1:num_events
   currIDset = cell(size(event_data{idx}));
   
   currIDset = cellfun(@(x,y)... %cellfun setup 
       [cellstr(repmat(x,num_celltypes,1)),cellstr(repmat(y,num_celltypes,1))],... %function
        modelIDs,stimIDs,'UniformOutput',false)'; %input
    
   currIDset = cellfun(@(x,y) [x,y],currIDset,cellIDs,'UniformOutput',false); %add celltypes
   dataIDs{idx} = currIDset;
end
