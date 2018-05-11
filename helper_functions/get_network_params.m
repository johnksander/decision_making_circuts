function options = get_network_params(do_config,options)
%takes network ID number, gives back parameter set

switch do_config
    case 1
        options.ItoE = 1.2904; options.EtoI = 0.1948; Rstim = 75 *.25; options.stim_targs = 'Eswitch';
    case 2
        options.ItoE = 1.2877; options.EtoI =  0.1734; Rstim = 100; options.stim_targs = 'Estay';
    case 3
        options.ItoE = 0.8069; options.EtoI = 0.2067; Rstim = 49.7632 *.25; options.stim_targs = 'Eswitch';
    case 4
        options.ItoE = 0.7936; options.EtoI = 0.1878; Rstim = 59.2903; options.stim_targs = 'Estay';
    case 5
        options.ItoE = 0.3679; options.EtoI = 0.2737; Rstim = 30.0436/2; options.stim_targs = 'Eswitch';
    case 6
        options.ItoE = 0.3161; options.EtoI = 0.2482; Rstim = 175; options.stim_targs = 'Estay';
    case 7
        options.ItoE = 0.2800; options.EtoI = 0.4228; Rstim = 47.0955 *.25; options.stim_targs = 'Eswitch';
    case 8
        options.ItoE = 0.2355; options.EtoI = 0.4250; Rstim = 175; options.stim_targs = 'Estay';
    case 9
        options.ItoE = 0.2921; options.EtoI = 0.6927; Rstim = 45/2; options.stim_targs = 'Eswitch';
    case 0
        options.ItoE = 0.2119; options.EtoI = 0.6799; Rstim = 175; options.stim_targs = 'Estay';
    otherwise
        error('config disaser')
end
options.trial_stimuli = [Rstim,Rstim];



%historical log
%original parameter set:
%     case 1
%         options.ItoE = 1.2904; options.EtoI = 0.1948; Rstim = 75; options.stim_targs = 'Eswitch';
%     case 2
%         options.ItoE = 1.2877; options.EtoI =  0.1734; Rstim = 100; options.stim_targs = 'Estay';
%     case 3
%         options.ItoE = 0.8069; options.EtoI = 0.2067; Rstim = 49.7632; options.stim_targs = 'Eswitch';
%     case 4
%         options.ItoE = 0.7936; options.EtoI = 0.1878; Rstim = 59.2903; options.stim_targs = 'Estay';
%     case 5
%         options.ItoE = 0.3679; options.EtoI = 0.2737; Rstim = 30.0436; options.stim_targs = 'Eswitch';
%     case 6
%         options.ItoE = 0.3161; options.EtoI = 0.2482; Rstim = 175; options.stim_targs = 'Estay';
%     case 7
%         options.ItoE = 0.2800; options.EtoI = 0.4228; Rstim = 47.0955; options.stim_targs = 'Eswitch';
%     case 8
%         options.ItoE = 0.2355; options.EtoI = 0.4250; Rstim = 175; options.stim_targs = 'Estay';
%     case 9
%         options.ItoE = 0.2921; options.EtoI = 0.6927; Rstim = 45; options.stim_targs = 'Eswitch';
%     case 0
%         options.ItoE = 0.2119; options.EtoI = 0.6799; Rstim = 175; options.stim_targs = 'Estay';

