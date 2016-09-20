function GISTIC_armfigures(a)
% a= struct of : 
% -ui path to input.mat -GISTIC
% GISTICout.mat -p param file -[-basedir base_directory] 
% gp_GISTIC_outfiles  -ui path to input.mat -GISTIC
% GISTICout.mat -p param file -[-basedir base_directory] -n number of samples
% ---
% $Id$
% $Date: 2007-09-17 15:18:13 -0400 (Mon, 17 Sep 2007) $
% $LastChangedBy: barbara$
% $Rev$

addpath ~/CancerGenomeAnalysis/trunk/matlab/
addpath ~/CancerGenomeAnalysis/trunk/matlab/gp_modules
addpath ~/CancerGenomeAnalysis/trunk/matlab/snp

% at least one of them has a base_dir (output_dir from previous module)
if ~isempty(a.basedir)
    base_dir = a.basedir;
    params=[];
else
    if ~isempty(a.p)
        params=read_params_file(a.p);
        pidx=grep('^output_dir$',{params.param},1);
        if isempty(pidx)
          error('No previous output_dir in params file');
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
  
if ~isempty(a.ui)
	infile=a.ui;
    params=[];
else
    if ~isempty(a.p)
        params=read_params_file(a.p);
        pidx=grep('^input_file_name$',{params.param},1);
		pbd=grep('^base_dir$',{params.param},1);
        if isempty(pidx)
          error('No previous struct filename in params file');
        else
          struct_dir=params(pbd(1)).value;
		  filename=params(pidx(1)).value;
		  infile=[struct_dir filename];
        end
    else
        error('Need to provide either a full unix path to struct or parameter file!');
    end
end
disp(['Reading input file: ' infile]);
x=load(infile);

if isempty(x)
   error('No variables in input file');
end
 
  fn=fieldnames(x);
  D=getfield(x,fn{1});

% --> 371 samples
refgene_file  ='~gadgetz/projects/snp/data/Refgene/hg17/ucsc_20070131/hg17_20070131.mat';
load(refgene_file);% read rg and cyto
rg=add_chrn(rg);
rg=mark_futreal_genes(rg);

% plot histograms
% FIGURE 1B
[qvs1,Y_1,pvs1,ampdel1]=plot_arm_histograms([base_dir 'chromosome_hists'],D,cyto,'median',[-0.5 0.5],0.025,[0.1 0.1]);
sig_arms=find(qvs1<=1e-2);


%%% calculate standard deviation USING ALL ARMS
Y2_1=Y_1;
Y2_1.dat=2.^(Y_1.dat+1); % move to raw space
%Y2_1s=reorder_D_rows(Y2_1,sig_arms);
st=std(Y2_1.dat,0,1);
[sst,ssti]=sort(st,'descend');

D.scores=st;
CL21x=reorder_D_cols(D,ssti);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define thirds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%n=str2num(a.n);
n=size(CL21x.dat,2)
top=floor(n/3)
top_idx=1:top
mid_idx=top+1:(2*top)
bot_idx=(2*top)+1:n


save thirds.mat top_idx mid_idx bot_idx ssti

%CL21x.supacc=strvcat({'Tertile: 1-Top/2-Middle/3-Bottom','UsedN: 1-1/2-5'});
CL21x.supdesc='Tertile: 1-Top/2-Middle/3-Bottom';
CL21x.supdat=[];
CL21x.supdat(1,top_idx)=1;
CL21x.supdat(1,mid_idx)=2;
CL21x.supdat(1,bot_idx)=3;
%CL21x.supdat(2,:)=1+double(cellfun('length',CL21x.used_normals)>1);
%CL21x=add_supmark(CL21x);

%%%%% FIGURE 1A
display_D(CL21x,[],[],'snpsup');
print_D([base_dir 'Figure1A'],{{'fig'},{'pdf'},{'png','-r180'}});

