function GISTIC_outfiles(a)
% a= struct of : 
% -ui path to input.mat -GISTIC
% GISTICout.mat -p param file -[-basedir base_directory] 
% -lout all_lesions.txt  -sout scores.txt -add_vals 1expandsall_lesion
%  -[ -name my_struct ] -[-refgene_file  hg17_20070131.mat] -[-annot_file
%  geneannotations] -[-mark annotation_marker]
% ---
% $Id$
% $Date: 2007-12-07 19:38:07 -0500 (Fri, 07 Dec 2007) $
% $LastChangedBy: jdobson $
% $Rev$
  
% addpath ~/CancerGenomeAnalysis/trunk/matlab/
% addpath ~/CancerGenomeAnalysis/trunk/matlab/gp_modules
% addpath ~/CancerGenomeAnalysis/trunk/matlab/snp

method_st='outfiles';

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

  
if ischar(infile)
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
        error('Could not find variable in the input file');
    end
else
   D=infile;
end


if isempty(a.GISTIC) 
    error('No GISTIC output file');
end
disp('Reading output file from GISTIC');

Gcore=[base_dir a.GISTIC];
  G=load(Gcore);
 
 if isfield(G,'cyto')
    cyto=G.cyto;
 end
 if isfield(G,'regs')
    regs=G.regs;
 end
 if isfield(G,'ts')
    ts=G.ts;
 end
 if isfield(G,'pvs')
    pvs=G.pvs;
 end
 if isfield(G,'p')
    p=G.p;
 end
 if isfield(G,'q')
    q=G.q;
 end
 if isfield(G,'ext')
    ext=G.ext;
 end
 if isfield(G,'ads')
    ads=G.ads;
 end

 if length(D.pos)~=length(q{1})
  error('Structure and GISTIC output file lengths do not match');
end

if exist('a.refgene_file','var')
  load(a.refgene_file);
else
    disp('Loading hg17');
    load('~gadgetz/projects/snp/data/Refgene/hg17/ucsc_20070131/hg17_20070131.mat');
end

if isempty(a.lout)
   error('Missing output file name for all_lesions file');
end

lout=[base_dir a.lout '_' ext '.txt'];
disp(['Writing output all_lesions file: ' lout]);

 start_at=[];
 
 if ~exist('no_call_thresh','var')
     no_call_thresh=ts;
 end
 
 add_vals=str2num(a.add_vals);
 
 qv_thresh=0.25;
 for k=1:2
  score_thresh(k)=min(ads{k}(find(q{k}<=qv_thresh)));
end
for k=1:2
  min_cutoff(k)=min(ads{k}(find(q{k}<1))); % qv < 1
end
if ~exist('p_arm','var')
     p_arm=struct('method','scorecutoff','ads',{ads},'p_arm',0.5,'score_thresh',score_thresh,...
        'score_thresh_focal',score_thresh,'min_cutoff',min_cutoff);
end		
	
regs1=find_broad_regs(D,cyto,regs,p_arm);
f=fopen([base_dir 'focal_broad_' ext '.txt'],'w');
ampdel={'Amp','Del'};
for k=1:2
  for i=1:length(regs1{k})
  [st,chr,bp]=genomic_location(D,{[regs1{k}(i).st regs1{k}(i).en]},cyto,1,0);
  fprintf(f,'%s\n',[ ampdel{k} num2str(i) ')' st ' focal:' num2str(regs1{k}(i).focal) ' broad:' num2str(regs1{k}(i).broad)]);
  end
end
fclose(f);
	
 make_all_lesions_file(lout,D,regs,pvs,cyto,start_at,ts,no_call_thresh,add_vals,p_arm);

disp('Adding gene annotations'); 
 rg=add_chrn(rg);
 annot_file=a.annot_file;
 mark=a.mark;
 rg=add_annotation(rg,annot_file,mark);
 calls=call_regs(D,regs,{ts});
    
disp('Writing output gene table files'); 
 A=[];
 u95ll=[];
 genetables(D,A,u95ll,rg,cyto,regs,calls,ts,ext,lout,base_dir);
 
disp('Writing q value, p value and amp/del score file');
 if length(D.pos)==length(q{1})
     
 sout=[base_dir a.sout '_' ext '.txt'];
 write_score_file(sout,D,p,q,ads,ts);
 else
   disp(['Length of struct does not match length of scores!' size(D.pos,2) length(q{1})]);
 %    Score=[q,p,ads];
 end
copy_number_statistics(D, ext, base_dir)
  
param_file = [base_dir 'parameter_' method_st, '.txt'];
disp(['Writing parameter file to:' param_file]);
local_copy=1;
param_struct = struct('output_dir', base_dir,'output_lesions_filename',a.lout,'output_scores_filename', a.sout, 'input_struct_filename',infile ,'base_dir', base_dir,'input_GISTICcore_filename', a.GISTIC ,'parameter_input_file',a.p, 'add_vals', add_vals, 'my_struct_name', a.name, 'refgene_file', a.refgene_file, 'annotation_file', a.annot_file, 'annotation_marker_symbol', a.mark);
gp_write_params(method_st,param_struct,params,param_file,local_copy);
