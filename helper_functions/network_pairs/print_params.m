function print_params(FN)
%This takes the network pair file & prints the table to a text file.
%Intended for setting up function to make get_network_params formatted
%switch statement


params = load(FN);
params = params.network_pairs;
nettypes = cellfun(@(x) x.Row,params,'UniformOutput',false);
%set the rows as a new var
params = cellfun(@(x,y) [x,cell2table(y,'VariableNames',{'type'})],...
    params,nettypes,'UniformOutput',false);

for idx = 1:numel(params)
    T = params{idx};
    T.Properties.RowNames(:) = ''; %delete row names
    params{idx} = T;
end

params = cat(1,params{:});
num_nets = size(params,1);

outFN = sprintf('code4%s.txt',FN);
outFN = strrep(outFN,'.mat','');

f = fopen(outFN,'w');
for idx = 1:num_nets
    
    P = params(idx,:);
    
    switch P.type{:}
        case 'fast'
            stimtargs = 'Estay';
        case 'slow'
            stimtargs = 'Eswitch';
    end
    
    pstr = ...
        sprintf('options.ItoE = %.4f; options.EtoI = %.4f; Rstim = []; options.stim_targs = ''%s'';',...
        P.ItoE,P.EtoI,stimtargs);
    
    fprintf(f,'case %i\n',idx);
    fprintf(f,'%s\n',pstr);
    
end


fclose(f);


% FN = [FN '.csv']; %append file ext.
% Nrow = size(data,1);
% Ncol = size(data,2);
%
% f = fopen(FN,'w');
% for i = 1:Nrow
%     for j = 1:Ncol
%         entry = data{i,j};
%         if ischar(entry)
%             fmt = '%s';
%         else
%             fmt = '%f';
%         end
%         if j == Ncol
%             dlm = '\n';
%         else
%             dlm = ',';
%         end
%         fprintf(f,[fmt,dlm],entry);
%     end
% end
% fclose(f);
%



%now write in this format...

% switch do_config
%     case 1
%         options.ItoE = []; options.EtoI = []; Rstim = []; options.stim_targs = [];
%     case 2
%         options.ItoE = []; options.EtoI = []; Rstim = []; options.stim_targs = [];
%     case 3
%         options.ItoE = []; options.EtoI = []; Rstim = []; options.stim_targs = [];
%     case 4
%         options.ItoE = []; options.EtoI = []; Rstim = []; options.stim_targs = [];
%     case 5
%         options.ItoE = []; options.EtoI = []; Rstim = []; options.stim_targs = [];
%     case 6
%         options.ItoE = []; options.EtoI = []; Rstim = []; options.stim_targs = [];
%     case 7
%         options.ItoE = []; options.EtoI = []; Rstim = []; options.stim_targs = [];
%     case 8
%         options.ItoE = []; options.EtoI = []; Rstim = []; options.stim_targs = [];
%     case 9
%         options.ItoE = []; options.EtoI = []; Rstim = []; options.stim_targs = [];
%     case 0
%         options.ItoE = []; options.EtoI = []; Rstim = []; options.stim_targs = [];
%     otherwise
%         error('config disaser')
% end



%outFN = sprintf('%s_table.txt',FN);
%writetable(params,outFN)