function Wnew = reorder_weightmat(W,celltype)
%reorder weight matrix so you can calculate cell inputs like this:
%I = I + (unique(Erev(celltype.excit)) - V(:,idx-1)).*(W(:,celltype.excit)*Sg(celltype.excit,idx-1)).*unique(Gg(celltype.excit));
%That means columns of W should give you all cell inputs from columntype.
%for example, W(:,celltype.excit) needs to give an Ncells by excitatory cells
%weight matrix with excitatory cells' output weights to every cell 

Wnew = NaN(size(W));
Wnew(celltype.excit,celltype.excit) = W(celltype.excit,celltype.excit);
Wnew(celltype.inhib,celltype.excit) = W(celltype.excit,celltype.inhib)';
Wnew(celltype.inhib,celltype.inhib) = W(celltype.inhib,celltype.inhib);
Wnew(celltype.excit,celltype.inhib) = W(celltype.inhib,celltype.excit)';
