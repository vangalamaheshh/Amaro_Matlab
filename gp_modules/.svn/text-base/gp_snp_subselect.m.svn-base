function gp_snp_subselect(varargin)
% gp_snp_subselection -b base_dir -i struct_filename -o output_dir_ext
% -of output_file_name -p parameter_file [-x should_remove_x -cnv cnv_file] 
% ---
% $Id$
% $Date$
% $LastChangedBy$
% $Rev$ 
  
addpath ~/CancerGenomeAnalysis/trunk/matlab
addpath ~/CancerGenomeAnalysis/trunk/matlab/snp
addpath ~/CancerGenomeAnalysis/trunk/matlab/gp_modules

a=handle_args({'b','i','o','of','x','cnv','p'},varargin);
method_st = 'subselect';

if ~isempty(a.b)
    base_dir = a.b;
else
    if ~isempty(a.p)
        params=read_params_file(a.p);
        pidx=grep('^output_dir$',{params.param},1);
        if isempty(pidx)
          error('No base_dir in params file');
        else
          base_dir=params(pidx(1)).value;
        end
    else
        error('Need to provide either a base directory or parameter file!');
    end
end

if base_dir(end)~='/'
  base_dir=[base_dir '/'];
end

if ~isempty(a.i)
    input_struct = a.i;
else
    error('Need to provide input structure name!')
end

if ~isempty(a.o)
    output_dir_ext = a.o;
else 
    output_dir_ext = 'output';
end

output_dir = [method_st '_' output_dir_ext '/'];
output_path = [base_dir output_dir];
if ~exist(output_path)
    mkdir(output_path);
else
    error('Output directory already exists.  Cannot overwrite.');
end

if ~isempty(a.of)
    outfile = a.of;
else
    outfile = [method_st '.struct.mat'];
end

out = [output_path outfile];

if ~isempty(a.x)
    should_use_X = str2num(a.x);
else 
    should_use_X = 0;
end

if isempty(a.cnv)
    cnv_file = [];
else
    cnv_file = a.cnv;
end
    
input_file = [ base_dir input_struct];
tmp=load(input_file);
nms=fieldnames(tmp);
if length(nms)>1
  error('The input file has more than one variable');
else
  CL=getfield(tmp,nms{1});
end

if ~should_use_X
  CL=reorder_D_rows(CL,find(CL.chrn~=23));
else
  disp('Are you sure you want to use X (think of normalization)?');
end

if ~isempty(cnv_file)
    CL = remove_cnv(CL,cnv_file);
end

disp(['Writing output to file:' out]);
eval([ nms{1} '=CL;']); % save a variable in the same name as input
save(out,nms{1},'-v7.3');

param_file = [output_path 'parameter.txt']

disp(['Writing parameter file to:' param_file]);

param_struct = struct('output_dir',output_path,'output_file_name',outfile,'base_dir',base_dir,'input_file_name',input_struct,'should_use_x',num2str(should_use_X),'cnv_file',cnv_file,'parameter_input_file',a.p);
gp_write_params(method_st,param_struct,params,param_file);
