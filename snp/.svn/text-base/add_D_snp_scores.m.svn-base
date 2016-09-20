function [D,supids]=add_D_snp_scores(D,cyto,gsupid,prefix)
%[D,supids]=add_D_snp_scores(D,cyto,gsupid,prefix)

if ~exist('gsupid','var') 
  gsupid=[];
end

if ~exist('prefix','var') 
  prefix='';
end

if ~exist('cyto','var') 
  cyto=[];
end

supids=[];

if ~isempty(gsupid)
  [D,sid]=add_D_snp_scores(D,cyto,[],prefix);
  supids=[supids sid];
  if ~iscell(D.gsupacc)
    D.gsupacc=as_column(cellstr(D.gsupacc));
  end
  if ~iscell(D.gsupdesc)
    D.gsupdesc=as_column(cellstr(D.gsupdesc));
  end
  
  [ttl,types]=break_sup_names(D.gsupacc{gsupid});
  for j=1:length(types)
      D1 = copyD(D);
    D1=reorder_D_rows(D1,find(D.gsupdat(gsupid,:)==j));
    [D1,sid]=add_D_snp_scores(D1,cyto,[],[ prefix types{j} '|']);
    [D,sid2]=add_D_sup(D,D1.supacc(sid),D1.supdesc(sid),D1.supdat(sid,:),'cols',1);
    supids=[supids sid2];
  end
else
  %% median
  verbose('Calculating median of autosomes',30);
  autosomes=find(D.chrn<=22);
    D1 = copyD(D);
 %   disp('Copying D')
    if strcmp(class(D),'datastruct')
        D1 = setwriteprotect(D1,'off');
    end
  D1=reorder_D_rows(D1,autosomes);   % get rid of chromosome 23
  
  %Potential memory hog
  %s=nanmedian(D1.dat,1);             % find median value for each sample
  s = itrfcn1(D1,'dat',1,@nanmedian,1);
  
  
  [D,sid]=add_D_sup(D,[prefix 'MED'], [prefix 'Median on autosomes'],s,'cols',1);
  supids=[supids sid];
  
  %% copy quality score
  if exist('cyto','var') && ~isempty(cyto)
    verbose('Calculating copy quality signal and noise',30);
    [q,st,s,Y2,cn_vs_snp_noise]=copy_signal_to_noise(D1,[],cyto);
    if (0)
      armst={'p','q'};
      Y2.gacc=strcat(Y2.chr,armst(Y2.armn)',cellstr(repmat('-',size(Y2.dat,1),1)),cellstr(num2str(Y2.n_collapse)));
      Y2.gdesc=Y2.gacc;
%      write_eisen_dat('/xchip/tcga/gbm/analysis/pipeline/plates/TRIBE/QSdetail.mat.txt',Y2.gacc,Y2.gdesc,Y2.sdesc,'ARMS',Y2.dat,[],[],0);
%      write_eisen_dat('stdout',Y2.gacc,Y2.gdesc,Y2.sdesc,'ARMS',Y2.dat,[],[],0);
    end
    
    %%%  TEMPORARY FIX is cn_vs_snp_noise is returned empty
    if isempty(cn_vs_snp_noise) 
        cn_vs_snp_noise = nan(1,length(q));
    end
    %%%
    
    [D,sid]=add_D_sup(D,{[prefix 'SIG'],[prefix 'NOISE'],[prefix 'S2N'],[prefix 'CNVSSNPN']},...
                      {[prefix 'Signal (based on autosome arms)'],[prefix 'Noise'],[prefix 'Signal to noise'],[prefix 'CN vs. SNP Noise']},...
                      [st;s;q;cn_vs_snp_noise] ,'cols',1);
    supids=[supids sid];
  else
    verbose('Calculating copy quality noise',30);
    s=calc_copy_quality_score(D1,'medianabs',0);
    [D,sid]=add_D_sup(D,[prefix 'NOISE'],...
                      [prefix 'Noise (on autosomes)'],...
                      s ,'cols',1);
    supids=[supids sid];
  end
end
