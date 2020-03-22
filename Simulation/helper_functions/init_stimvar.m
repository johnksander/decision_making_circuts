function stim_info = init_stimvar(celltype,pool_options,options)
%---stimuli info--------------
%the structure of this thing should be... 
%targ_cells is a Ncells x Ntargs matrix
%stim_info.(stimulus label) is a 1 x Ntargs vector 

Ncells = pool_options.num_cells; %only reason this gets passed here..
timestep = options.timestep;
stim_labels = {'stim_A','stim_B'}; %matches init_statevar() and timepoint_stimulus()
stim_info = struct();
stim_info.num_targs = numel(options.stim_targs);
stim_info.targ_cells = false(Ncells,stim_info.num_targs);
for idx = 1:numel(stim_labels)
    stim_info.(stim_labels{idx}) = NaN(1,stim_info.num_targs);
end

for idx = 1:stim_info.num_targs
    switch options.stim_targs{idx}
        case 'Eswitch'
            stim_info.targ_cells(:,idx) = celltype.excit & celltype.pool_switch; %Eswitch cells
        case 'Estay'
            stim_info.targ_cells(:,idx) = celltype.pool_stay & celltype.excit; %Estay
        case 'baseline'
            stim_info.targ_cells(:,idx) = false(Ncells,1); %no targets
    end
    Rstim = options.trial_stimuli{idx};
    %careful, each cell in trial_stimuli contains 1 x Nstims vector (not 1 X Ntargets !!)
    for j = 1:numel(stim_labels)
        lambda = Rstim(j) * timestep; %poisson lambda for stimulus conductance;
        %index as stim_info.stim_A(target) or stim_info.(label{i})(target) etc
        stim_info.(stim_labels{j})(idx) = lambda;
    end
end

stim_info.num_cells = Ncells; %just so I don't have to pass pool_options as well later somewhere... 
if all(~isnan(options.stim_pulse))
    stim_info.delivery = 'pulse';
    stim_info.pulse = options.stim_pulse ./ timestep;
    stim_info.sample_schedule = options.stim_schedule;
else
    stim_info.delivery = 'constant';
end

