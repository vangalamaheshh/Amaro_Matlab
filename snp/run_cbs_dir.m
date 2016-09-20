function run_cbs_dir(dir_name,prefix_name,lsf_queue,redo,undoSD)

% run cbs through LSF

warning(['This function has been deprecated. Please use run_cbs ' ...
         'instead. NS 2009.02.09'] )

if ~exist('redo','var') || isempty(redo)
  redo=-1;
end

cd(dir_name);

if redo==-1
  dr=dir([ prefix_name '*Sample*']);
  tmp=grep('.*Sample[0-9]+$',{dr.name},1);
  dr=dr(tmp);
  tok=regexp({dr.name},'.*Sample([0-9]+)$','tokens');
  
  ne=~cellfun('isempty',tok);
  idx=find(ne);
  
  num=[];
  for i=1:length(idx)
    num(i)=str2num(tok{idx(i)}{1}{1});
  end

  [tmp,si]=sort(num);
  if length(si)~=max(num)
    error('not all samples exist');
  end
  
  sidx=idx(si);
  dr=dr(sidx);

  for i=1:length(dr)
    R_fname=[dr(i).name '.R'];
    disp(['Creating ' R_fname]);
    if exist('undoSD','var') && ~isempty(undoSD)
      run_R_file(R_fname,'/xchip/tcga/Tools/CBS/gp_CBSsample.new_output.R','gpCBS',[ '"' dr(i).name '","' dr(i).name ...
                          '.cbs",nPerm=10000,alpha=0.01,undo.splits="sdundo",undo.SD=' num2str(undoSD)  ',nacheck=T']);
    else
      run_R_file(R_fname,'/xchip/tcga/Tools/CBS/gp_CBSsample.new_output.R','gpCBS',[ '"' dr(i).name '","' dr(i).name ...
                          '.cbs",nPerm=10000,alpha=0.01,nacheck=T']);
    end
  end
else
  dr=dir([ prefix_name '*Sample*']);
  tmp=grep('.*Sample[0-9]+$',{dr.name},1);
  dr=dr(tmp);
  dr_names={dr.name};
  dr2=dir([ prefix_name '*Sample*.cbs']);
  dr2_names=regexprep({dr2.name},'\.cbs','');
  [Mt,m1,m2]=match_string_sets_hash(dr2_names,dr_names);
  delta=setdiff(1:length(dr),m2);
  if length(delta)>0
    dr=dr(delta);
  else
    return
  end
  catn({dr.name})
end

job=[];
if exist('undoSD','var') && ~isempty(undoSD)
  name=['_' prefix_name '_' num2str(undoSD)];
else
  name=['_' prefix_name ];
end
if ~exist('lsf_queue','var') || isempty(lsf_queue)
  lsf_queue='long';
end
for i=1:length(dr)
  R_fname=[dr(i).name '.R'];
  unixstr=['bsub ',...
           '-E ~/CancerGenomeAnalysis/trunk/shell/chk_lsf_cga ',...
           '-R "rusage[mem=4]" ',... 
           '-mig 5 ',...
           '-R "select[cpuf>100]" ',... 
           '-Q "EXCLUDE(127)" ',...
           '-P cbs_cancerfolk ',... 
           '-J CBS_' num2str(i) ' -q ' lsf_queue ' -o ' [dr(i).name '.log.txt'] ' -e ' [dr(i).name '.err.txt '] ,...
           '-r R CMD BATCH --no-save --no-restore ' R_fname ];
  disp(unixstr);
  [r1,r2]=unix(unixstr); 
  tok=regexp(r2,'\<([0-9]+)\>','tokens');
  if ~isempty(tok)
    job(i)=str2num(tok{1}{1});
  else
    job(i)=-1;
  end  
end
% wait_for_lsf_jobs(job);

