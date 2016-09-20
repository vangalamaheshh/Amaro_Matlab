function output_list_of_arglists = sg_feval(func,input_list_of_arglists,workingdir)
% SG_FEVAL  Native MATLAB wrapper around Scatter-Gather module
%
% Implements a simple but common use case.
%
% Inputs:
% func - string containing function name to run.
% input_list_of_arglists - length n cell array of cell arrays, where n jobs
%    will be run on LSF.  Each of the n cell arrays contains the arguments to
%    be passed into func on one node.  Huge variables are possible but will be
%    slow.
% workingdir - directory which is assumed to exist and be empty at the
%    start, and will be left full of disposable files at the end.
%
% Outputs:
% output_list_of_arglists - length n cell array of cell arrays, each one
%    containing the return values from one node.
%
% Example:
% !mkdir /tmp/foo
% out = sg_feval('sum',{ {[1 2 3]} , {[4 5 6]} },'/tmp/foo');
% out
%   [6]
%   [15]



% todo: ought to add flags for whether graphics used (ie jvm in compiler, graphics flag for
% xml file)
% todo ought to add items out of an existing config file, eg queue etc etc

use_lsf = 1;

if ~use_lsf
    % Loop through the items sequentially and locally.
    scatter_count = length(input_list_of_arglists);
    output_list_of_arglists = cell(1,scatter_count);
    fh = str2func(func);
    nout = nargout(fh);
    for i=1:scatter_count
        output_list_of_arglists{i} = cell(1,nout);
        [output_list_of_arglists{i}{:}] = feval(fh,input_list_of_arglists{i}{:});
    end
    
else
    % Write input args to working directory.  This will be the parent directory
    % for the scatter jobs.
    [scatter_count,scatter_indicies] = write_input_args(input_list_of_arglists,workingdir);
    
    % Write function wrapper to working directory.
    func_wrapped = generate_matfileargs_func (func,workingdir);
    % Compile function, store in working directory.
    % Would be nice to cache the compilation...
    %   (hash the files output by depfun and store it with executable, or
    %   compare dates of files output by depfun with the compiled file date.
    %   Also... store compiled version + hash or file timestamps in central location.)
    compile_mfile(func_wrapped,workingdir);
    
    % Write file to be used for Prepare
    % Add an item for the Gather stage.
    prepare_path = fullfile(workingdir,'prepare_output.txt');
    prepare_file_contents = [scatter_indicies,{'gather'}];
    write_cellstrs (prepare_path, prepare_file_contents);
    
    
    % Write XML config file
    xml_config_fn = generate_config_file (workingdir);
    
    % Write MATLAB helper file
    copy_matlab_helper(workingdir);
    
    % Launch Scatter-Gather module
    run_scatter_gather(xml_config_fn,workingdir, prepare_path);
    
    % Load results into cell array
    output_list_of_arglists = get_results(workingdir,scatter_count);
    
end

end

function [scatter_count,scatter_indicies] = write_input_args(input_list_of_arglists,workingdir)
    % The filenames used agree with what is expected by the wrapper.
    scatter_count = length(input_list_of_arglists);
    scatter_indicies = cell(1,scatter_count);
    for i = 1:scatter_count
        scatter_indicies{i} = num2str(i);
        fn = fullfile(workingdir,['inargs_' scatter_indicies{i} '.txt']);
        save(fn,input_list_of_arglists{i});
    end

end

function out_exe = compile_mfile(fn,outdir)

% This name aligns with the xml file.
out_exe = fullfile(outdir,'scatter_func');

% Find non-MATLAB directories to add to the compilers path
current_path = path;
pathlist = textscan(current_path,'%s','delimiter',':');
pathlist = pathlist{1};
patharg = '';
for i=1:length(pathlist)
    if isempty(strfind(pathlist{i},'toolbox'))
        patharg = sprintf('%s -I %s',patharg,pathlist{i});
    end
end

% Note: patharg must come before func_fn.
comp_exe = ['mcc -mv -d ' out_dir ' ' patharg ' -o ' out_exe ' ' fn  ];
% Prepend use's to make sure we are on the latest versions
cmd = ['reuse -q matlab;reuse -q GCC-4.4;' comp_exe];

% Call using command line, so that this module itself can be compiled.
system_check(cmd, 'compile failed');

end

function fn = generate_matfileargs_func (funcname_str,outdir)
% generate_matfileargs_func  Create file that wraps a matlab func to use matfiles for args.

wrapper_funcname_str = [funcname_str '_matfileargs'];
fn = fullfile(outdir, [wrapper_funcname_str '.m']);

% Generate contents of output file.
% Note '' indicates a single '
file_contents={
    'function ' wrapper_funcname_str '(scatter_index)'
    '% ' upper(wrapper_funcname_str) '(scatter_index)  Calls ' funcname_str ' with input and output args from matfiles'
    '% Loads inargs_<scatter_index>.mat from the parent directory.  This is expected to be a 1D'
    '%   cell array.  It is passed in as input args to the function.'
    '% The child function name is specified as a function handle in the '
    '%   first line of code, so the compiler can resolve it as a dependent.  '
    '% The output arguments are stored as a cell array in outargs.mat in the'
    '%   current directory.'
    ''
    '% Specify child function'
    'fh = @' funcname_str ';'
    ''
    '% Load input args'
    'inargs = load(''../inargs_' scatter_index '.mat'');'
    'if nargin(fh) ~= length(inargs)'
    '    error(''scatter_wrapper:InputArgNumMismatch'',...'
    '        ''Wrong number of input args.'');'
    'end'
    '% Reserve space for the max number of output args'
    'outargs = cell(1,nargout(fh));'
    ''
    '% Call user function'
    '[outargs{:}] = feval(fh,inargs{:});'
    ''
    '% Store output args.'
    'save(''outargs.mat'',outargs);'
    ''
    };

write_cellstrs(fn,file_contents);

end

function fn = generate_config_file (outdir)
% generate_matfileargs_func  Create file that wraps a matlab func to use matfiles for args.

fn = fullfile(outdir, 'configfile.xml');

% Generate contents of output file.
% Note '' indicates a single '
file_contents={
    '<scatter-gather>'
    '<prepare exe="sh cat"/>'
    '<scatter exe="scatter_func" memory="1" matlab="true" matlabDisplay="true"/>'
    '<gather exe="sh echo" memory="1" matlab="false"/>'
    '</scatter-gather>'
    };

write_cellstrs(fn,file_contents);

end

function copy_matlab_helper(workingdir)
    cmd = ['cp ~gsaksena/CancerGenomeAnalysis/trunk/analysis_pipeline/scripts/run_matlab.py ' workingdir ];
    system_check(cmd,'Copy Failed');

end

function run_scatter_gather(xml_config_fn,workingdir, prepare_output_fn)
    sg_exe = ['python scatter-gather.py broad ' xml_config_fn ' ' workingdir ' ' prepare_output_fn ];
    cmd = ['reuse -q python;' sg_exe];
    system_check(cmd,'Scatter-Gather failed: ');
    
end

function output_list_of_arglists = get_results(workingdir,scatter_count)
    d = dir ([workingdir '/scatter*']);
    scatterdirs = sort(cellstr(char(d.name)));
    dir_count = length(scatterdirs);
    if dir_count ~= scatter_count
        throw(MException('get_results:resultCountMismatch',['Expected ' num2str(scatter_count) ' results but found only ' num2str(dir_count)]));
    end
    output_list_of_arglists = cell(1,dir_count);
    for i=1:dir_count
        % filename aligns with wrapper script
        outarg_fn = fullfile(workingdir,scatterdirs{i},'outargs.mat');
        output_list_of_arglists{i} = load(outarg_fn);
    end 
end

function system_check(cmd,err_str)
    [status,result] = system(cmd);
    if status ~= 0
        throw(MException('system_check:CmdFailed',[err_str ' ' cmd ' ' result ]));
    end
    disp(result)

end

function write_cellstrs(file_path, file_contents)
% dump cellstr to output file
[fid,msg] = fopen (file_path,'w');
if fid == -1
    throw(MException('write_cellstrs:CantOpenWriteFile', ...
        ['Cannot open file ' file_path ' for writing' msg ]));
end
for i = 1:length(file_contents)
    fprintf(fid,[file_contents{i} '\n']);
end
fclose(fid);

end

