function options = get_network_params(do_config,options)
%takes network ID number, gives back parameter set

switch options.netpair_file
    
    case 'D2t'
        
        switch do_config
            case 1
                options.ItoE = 12.4541; options.EtoI = 0.2478; Rstim = 70.2; options.stim_targs = 'Eswitch';
            case 2
                options.ItoE = 12.4900; options.EtoI = 0.0892; Rstim = 366.7; options.stim_targs = 'Estay';
            case 3
                options.ItoE = 8.7866; options.EtoI = 0.3584; Rstim = 47.1; options.stim_targs = 'Eswitch';
            case 4
                options.ItoE = 8.7794; options.EtoI = 0.0996; Rstim = 359.6; options.stim_targs = 'Estay';
            case 5
                options.ItoE = 7.7681; options.EtoI = 0.7487; Rstim = 53.4; options.stim_targs = 'Eswitch';
            case 6
                options.ItoE = 3.9050; options.EtoI = 0.7500; Rstim = 650.4; options.stim_targs = 'Estay';
            case 7
                options.ItoE = 8.0928; options.EtoI = 0.4535; Rstim = 48.9; options.stim_targs = 'Eswitch';
            case 8
                options.ItoE = 3.7532; options.EtoI = 0.4543; Rstim = 475.8; options.stim_targs = 'Estay';
            case 9
                options.ItoE = 8.3703; options.EtoI = 0.3849; Rstim = 41.1; options.stim_targs = 'Eswitch';
            case 10
                options.ItoE = 4.6428; options.EtoI = 0.1632; Rstim = 315.5; options.stim_targs = 'Estay';
            otherwise
                error('config disaser')
        end
        
        
    case 'D2t_20s'
        
        switch do_config
            case 1
                options.ItoE = 12.4541; options.EtoI = 0.2478; Rstim = 19.3; options.stim_targs = 'Eswitch';
            case 2
                options.ItoE = 12.4900; options.EtoI = 0.0892; Rstim = 492.7; options.stim_targs = 'Estay';
            case 3
                options.ItoE = 8.7866; options.EtoI = 0.3584; Rstim = 5.4; options.stim_targs = 'Eswitch';
            case 4
                options.ItoE = 8.7794; options.EtoI = 0.0996; Rstim = 465.6; options.stim_targs = 'Estay';
            case 5
                options.ItoE = 7.7681; options.EtoI = 0.7487; Rstim = 8.9; options.stim_targs = 'Eswitch';
            case 6
                options.ItoE = 3.9050; options.EtoI = 0.7500; Rstim = 753.8; options.stim_targs = 'Estay';
            case 7
                options.ItoE = 8.0928; options.EtoI = 0.4535; Rstim = 7.6; options.stim_targs = 'Eswitch';
            case 8
                options.ItoE = 3.7532; options.EtoI = 0.4543; Rstim = 550.4; options.stim_targs = 'Estay';
            case 9
                options.ItoE = 8.3703; options.EtoI = 0.3849; Rstim = 3.7; options.stim_targs = 'Eswitch';
            case 10
                options.ItoE = 4.6428; options.EtoI = 0.1632; Rstim = 377.6; options.stim_targs = 'Estay';
            otherwise
                error('config disaser')
        end
        
    case 'slowD'
        error('old simulation, check what you''re doing')
        switch do_config
            case 1
                options.ItoE = 7.9608; options.EtoI = 0.2401; Rstim = 97.9; options.stim_targs = 'Eswitch';
            case 2
                options.ItoE = 7.9605; options.EtoI = 0.1066; Rstim = 169.3; options.stim_targs = 'Estay';
            case 3
                options.ItoE = 5.9226; options.EtoI = 0.3920; Rstim = 90.3; options.stim_targs = 'Eswitch';
            case 4
                options.ItoE = 5.9009; options.EtoI = 0.1240; Rstim = 182.5; options.stim_targs = 'Estay';
            case 5
                options.ItoE = 6.3216; options.EtoI = 0.7430; Rstim = 145.8; options.stim_targs = 'Eswitch';
            case 6
                options.ItoE = 3.4540; options.EtoI = 0.7420; Rstim = 378.2; options.stim_targs = 'Estay';
            case 7
                options.ItoE = 6.1513; options.EtoI = 0.4612; Rstim = 121.7; options.stim_targs = 'Eswitch';
            case 8
                options.ItoE = 3.3772; options.EtoI = 0.4647; Rstim = 308.7; options.stim_targs = 'Estay';
            case 9
                options.ItoE = 5.8175; options.EtoI = 0.3958; Rstim = 84.7; options.stim_targs = 'Eswitch';
            case 10
                options.ItoE = 3.7917; options.EtoI = 0.1792; Rstim = 286.5; options.stim_targs = 'Estay';
            otherwise
                error('config disaser')
        end
        
        
    case 'fastD'
        error('old simulation, check what you''re doing')
        switch do_config
            case 1
                options.ItoE = 3.6690; options.EtoI = 0.0740; Rstim = 47.6; options.stim_targs = 'Eswitch';
            case 2
                options.ItoE = 3.7135; options.EtoI = 0.0503; Rstim = 130.1; options.stim_targs = 'Estay';
            case 3
                options.ItoE = 2.7421; options.EtoI = 0.0919; Rstim = 35.6; options.stim_targs = 'Eswitch';
            case 4
                options.ItoE = 2.7116; options.EtoI = 0.0619; Rstim = 109.4; options.stim_targs = 'Estay';
            case 5
                options.ItoE = 1.9414; options.EtoI = 0.3474; Rstim = 30.7; options.stim_targs = 'Eswitch';
            case 6
                options.ItoE = 1.5844; options.EtoI = 0.3497; Rstim = 156.2; options.stim_targs = 'Estay';
            case 7
                options.ItoE = 1.9242; options.EtoI = 0.2273; Rstim = 30.7; options.stim_targs = 'Eswitch';
            case 8
                options.ItoE = 1.5063; options.EtoI = 0.2268; Rstim = 197.7; options.stim_targs = 'Estay';
            case 9
                options.ItoE = 2.0132; options.EtoI = 0.1429; Rstim = 19.8; options.stim_targs = 'Eswitch';
            case 10
                options.ItoE = 1.7233; options.EtoI = 0.1043; Rstim = 214.5; options.stim_targs = 'Estay';
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

