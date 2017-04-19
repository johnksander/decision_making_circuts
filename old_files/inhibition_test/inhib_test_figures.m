clear
clc
format compact 


basedir = 'C:\Users\jksander.000\Desktop\rotation\project';
addpath(fullfile(basedir,'helper_functions'))
testdir = fullfile(basedir,'inhibition_test');
vardir = fullfile(testdir,'output_vars');

fontsz = 16;
modelnames = {'PM','JK'};
num_models = numel(modelnames);
weight_params = [.5:.02:1.5];
num_params = numel(weight_params);

mdl_Irates = cell(num_models,1);
mdl_Erates = cell(num_models,1);

for mdlidx = 1:num_models
    
    Irates = NaN(num_params,1);
    Erates = NaN(num_params,1);
    for idx = 1:num_params
        
        resultfile = [modelnames{mdlidx} '_' strrep(num2str(weight_params(idx)),'.','p')];
        resultfile = fullfile(vardir,resultfile);
        results = load(resultfile);
        [Eon,Ioff] = inhib_test_vectors(results);
        Erates(idx) = nanmean(Eon);
        Irates(idx) = nanmean(Ioff);
        disp(sprintf('model-%s: weight %.2f rates calculated',modelnames{mdlidx},weight_params(idx)))
    end
    mdl_Irates{mdlidx} = Irates;
    mdl_Erates{mdlidx} = Erates;
end



for mdlidx = 1:num_models
    plot(weight_params,mdl_Irates{mdlidx},'linewidth',2)
    hold on
    plot(weight_params,mdl_Erates{mdlidx},'linewidth',2)
    legend(char({'inactive pool','inhibitory cells'}),...
        char({'active pool','excitatory cells'}),'location','northwest')
    ylabel('mean spike rate')
    xlabel('Current (proportion)')
    title([modelnames{mdlidx} ' model'])
    set(gca,'Fontsize',fontsz)
    hold off 
    orient landscape
    fig_fn = [modelnames{mdlidx} '_spikerates'];
    print(fullfile(testdir,fig_fn),'-djpeg')
end














