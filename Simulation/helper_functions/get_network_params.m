function options = get_network_params(do_config,options)
%takes network ID number, gives back parameter set

switch options.netpair_file
    
    case 'slowD'
        
        switch do_config
            case 1
                options.ItoE = 7.7057; options.EtoI = 0.2640; Rstim = []; options.stim_targs = 'Eswitch';
            case 2
                options.ItoE = 7.7374; options.EtoI = 0.1112; Rstim = []; options.stim_targs = 'Estay';
            case 3
                options.ItoE = 6.8143; options.EtoI = 0.3092; Rstim = []; options.stim_targs = 'Eswitch';
            case 4
                options.ItoE = 6.8697; options.EtoI = 0.1207; Rstim = []; options.stim_targs = 'Estay';
            case 5
                options.ItoE = 6.3216; options.EtoI = 0.7430; Rstim = []; options.stim_targs = 'Eswitch';
            case 6
                options.ItoE = 3.4540; options.EtoI = 0.7420; Rstim = []; options.stim_targs = 'Estay';
            case 7
                options.ItoE = 5.9138; options.EtoI = 0.5039; Rstim = []; options.stim_targs = 'Eswitch';
            case 8
                options.ItoE = 3.2364; options.EtoI = 0.5128; Rstim = []; options.stim_targs = 'Estay';
            case 9
                options.ItoE = 5.9226; options.EtoI = 0.3920; Rstim = []; options.stim_targs = 'Eswitch';
            case 10
                options.ItoE = 3.4651; options.EtoI = 0.2852; Rstim = []; options.stim_targs = 'Estay';
            otherwise
                error('config disaser')
        end
        
    otherwise
        error('config disaser')
end

options.trial_stimuli = [Rstim,Rstim];

%historical log
%original parameter set:
% case 1
%     options.ItoE = 1.2904; options.EtoI = 0.1948; Rstim = 75 *.25; options.stim_targs = 'Eswitch';
% case 2
%     options.ItoE = 1.2877; options.EtoI =  0.1734; Rstim = 100; options.stim_targs = 'Estay';
% case 3
%     options.ItoE = 0.8069; options.EtoI = 0.2067; Rstim = 49.7632 *.25; options.stim_targs = 'Eswitch';
% case 4
%     options.ItoE = 0.7936; options.EtoI = 0.1878; Rstim = 59.2903; options.stim_targs = 'Estay';
% case 5
%     options.ItoE = 0.3679; options.EtoI = 0.2737; Rstim = 30.0436/2; options.stim_targs = 'Eswitch';
% case 6
%     options.ItoE = 0.3161; options.EtoI = 0.2482; Rstim = 175; options.stim_targs = 'Estay';
% case 7
%     options.ItoE = 0.2800; options.EtoI = 0.4228; Rstim = 47.0955 *.25; options.stim_targs = 'Eswitch';
% case 8
%     options.ItoE = 0.2355; options.EtoI = 0.4250; Rstim = 175; options.stim_targs = 'Estay';
% case 9
%     options.ItoE = 0.2921; options.EtoI = 0.6927; Rstim = 45/2; options.stim_targs = 'Eswitch';
% case 0
%     options.ItoE = 0.2119; options.EtoI = 0.6799; Rstim = 175; options.stim_targs = 'Estay';

