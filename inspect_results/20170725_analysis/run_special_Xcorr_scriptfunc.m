clear
clc
close all
format compact

cell_labels = {'E-stay','I-stay','E-switch','I-switch'};
num_cells = numel(cell_labels);
%targIDs{1} = {{'pre-switch','stimulus'};{'pre-switch','stimulus'}};
% targIDs{2} = {{'stable','stimulus'};{'post-switch','stimulus'}};
% targIDs{3} = {{'pre-switch','stimulus'};{'post-switch','stimulus'}};
targIDs = {{'pre-switch','stimulus'}};
num_comparisons = numel(targIDs);
Xcor_method{1} = 'coeff';
Xcor_method{2} = 'unbiased';
num_methods = numel(Xcor_method);

cellpairs = nchoosek(cell_labels,2);
num_cellpairs = numel(cellpairs(:,1));


for methidx = 1:num_methods
    for pairidx = 1:num_cellpairs
        for compidx = 1:num_comparisons
            special_Xcor_scriptfunc(cellpairs(pairidx,:),targIDs,Xcor_method{methidx});
        end
    end
end