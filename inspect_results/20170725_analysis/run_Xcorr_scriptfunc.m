clear
clc
close all
format compact

cell_labels = {'E-stay','I-stay','E-switch','I-switch'};
num_cells = numel(cell_labels);
targIDs{1} = {{'stable','stimulus'};{'pre-switch','stimulus'}};
targIDs{2} = {{'stable','stimulus'};{'post-switch','stimulus'}};
targIDs{3} = {{'pre-switch','stimulus'};{'post-switch','stimulus'}};
num_comparisons = numel(targIDs);
Xcor_method{1} = 'coeff';
Xcor_method{2} = 'unbiased';
num_methods = numel(Xcor_method);


for methidx = 1:num_methods
    for cellidx = 1:num_cells
        for compidx = 1:num_comparisons
            Xcorr_analysis_scriptfunc(cell_labels{cellidx},targIDs{compidx},Xcor_method{methidx});
        end
    end
end