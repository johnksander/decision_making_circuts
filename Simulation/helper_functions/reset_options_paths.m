function options = reset_options_paths(options)
%convert bigmem filepaths in options file to linus filpathes

bigmem_basedir = '/data/netapp/jksander/rotation/Simulation/';
bender_basedir = '/Users/ksander/Desktop/work/ACClab/rotation/project/';
%bender_basedir = 'C:/Users/jksander.000/Desktop/rotation/project/';

optfields = fieldnames(options);
for idx = 1:numel(optfields)
    
    data = getfield(options,optfields{idx});
    if isstr(data)
        path2fix = strfind(data,bigmem_basedir);
        if ~isempty(path2fix)
            data = strrep(data,bigmem_basedir,bender_basedir);
        end
        options = setfield(options,optfields{idx},data);
    end
end

