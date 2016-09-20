%% matlab script for determining gistic dependencies

% use depfun to get dependencies for GISTIC2 gene pattern modules
[list,builtins,classes] = depfun('gp_gistic2_from_seg');

% separate our gistic m-files from matlab's 
mathworks_fun_idx = strmatch('/broad/software/nonfree',list);
gistic_fun_idx = setdiff(1:length(list),mathworks_fun_idx);
gistic_list = regexprep(list(gistic_fun_idx),'/home/radon0./[A-Za-z0-9_]+/','~/');

% write out results
script_name = 'gen_gistic2_depend_links';
fid = fopen(script_name,'w');
if fid ~= -1
    cellfun(@(s) fprintf(fid,'ln -s %s\n',s),gistic_list);
    fclose(fid);
end
unix(sprintf(['chmod 775 ',script_name]));

%{
manual post-processing
create @SegArray directory, cd to it, and put @SegArray links there
remove:
    mike/*.m
    xunit/Contents.m
%}
