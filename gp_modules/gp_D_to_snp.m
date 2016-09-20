function gp_D_to_snp(varargin)
% gp_D_to_snp -i input.mat -ui /xchip/.... -o output.snp [ -name my_struct ]
% ---
% $Id$
% $Date: 2007-09-18 13:20:16 -0400 (Tue, 18 Sep 2007) $
% $LastChangedBy: rameen $
% $Rev$
%Barbara Weir revision: 10-25-07
%allows user to input flag for output type using write_as_dchip function:  
%   outflag of 2 corresponds to gp_format(in write_as_dchip fxn)  of -2 , 1 corresponds to gp_format of -1.


addpath ~/CancerGenomeAnalysis/trunk/matlab/
addpath ~/CancerGenomeAnalysis/trunk/matlab/gp_modules
addpath ~/CancerGenomeAnalysis/trunk/matlab/snp

a=handle_args({'i','ui','o','name','outflag'},varargin);

if isempty(a.i) && isempty(a.ui)
  error('Missing input file name');
  end

% at least one of them has a filename
if ~isempty(a.i) % a.i overrides a.ui
  infile=a.i;
else
  infile=a.ui;
end
  outflag=str2num(a.outflag);
if outflag==2
    flag=-2;
else if outflag==1
        flag=-1;
    else
        flag=outflag;
    end
end
disp(['Reading input file: ' infile]);
x=load(infile);

if isempty(x)
   error('No variables in input file');
end
 
if isempty(a.name)
  fn=fieldnames(x);
  D=getfield(x,fn{1});
elseif isfield(x,a.name)
  D=getfield(x,a.name);
else
   error('Could not find varaible in the input file');  
end

if isempty(a.o)
   error('Missing output file name');
end

disp(['Writing output file: ' a.o]);
write_as_dchip(a.o,D,1,[],flag);

% for compilation:
% cd ~/matlab/gp_modules
% mcc -m gp_D_to_snp
