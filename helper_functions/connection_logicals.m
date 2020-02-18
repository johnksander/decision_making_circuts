function [EtoE,EtoI,ItoE] = connection_logicals(celltype,num_cells)
%make logical matricies for different connections based on 
%celltype logical vectors
%inputs: number of total cells, logical vectors for pool type and
%excitatory/inhibitory cell type. 
%outputs: logical matricies for excitatory-excitatory,
%excitatory-inhibitory, and inhibitory-excitatory based on scheme specified
%in this function. 


EtoE =  repmat(celltype.excit,1,num_cells) & repmat(celltype.excit',num_cells,1);
EtoI =  repmat(celltype.excit,1,num_cells) & repmat(celltype.inhib',num_cells,1);
ItoE =  repmat(celltype.inhib,1,num_cells) & repmat(celltype.excit',num_cells,1);
self_connection = repmat(celltype.pool_stay,1,num_cells) == repmat(celltype.pool_stay',num_cells,1);
EtoE = EtoE & self_connection; %only self-celltype.excitation within group
EtoI = EtoI & ~self_connection; %E to I connection is across group
ItoE = ItoE & self_connection; %I to E connection is within group

%maybe figure out if these need to be used later
%celltype.pool_stay_mat = repmat(celltype.pool_stay,1,num_cells) | repmat(celltype.pool_stay',num_cells,1); 
%celltype.pool_switch_mat = repmat(celltype.pool_switch,1,num_cells) | repmat(celltype.pool_switch',num_cells,1);
