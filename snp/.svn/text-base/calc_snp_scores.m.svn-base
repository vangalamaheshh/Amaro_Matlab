function scores=calc_snp_scores(D,cyto,gsupid);


%% add scores per probe type

% median
%%%% FIXME
if 

disp('median');
autosomes=find(D.chrn<=22);
scores.median=nanmedian(D.dat(autosomes,:),1);

% copy quality
if exist('cyto','var') && ~isempty(cyto)
  disp('copy quality and signal to noise');
  [q,st,s,Y2]=copy_signal_to_noise(D,[],cyto);
  scores.copy_quality=s;
  scores.signal_to_noise=q;
  scores.signal=st;
else
  disp('copy quality');
  scores.copy_quality=calc_copy_quality_score(D,'medianabs',0);
end

