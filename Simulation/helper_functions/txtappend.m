function txtappend(fname,message)
% Append to a text file.
% Particularly useful for keeping track of progress in parfor loops.
% 
% CRM 20150117

if nargin == 1
    message = [datestr(now,31) '\n'];
end

f = fopen(fname, 'a');
fprintf(f, message);
fclose(f);

%disp(sprintf('Wrote to %s: %s',fname,message))