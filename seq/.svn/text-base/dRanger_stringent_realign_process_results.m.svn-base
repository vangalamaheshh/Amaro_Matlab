function dRanger_stringent_realign_process_results(sample,P)
% dRanger_stringent_realign_process_results(sample,P)
%
% Mike Lawrence 2009-07-23

if ~exist('sample','var'), error('<sample> is required'); end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P = impose_default_value(P,'cancer_sample',true);
P = impose_default_value(P,'tumor_only',false);
P = impose_default_value(P,'normal_only',false);
P = impose_default_value(P,'vicinity_forward',10000);
P = impose_default_value(P,'vicinity_backward',100);
P = impose_default_value(P,'dRangerPreprocess_output_dir_suffix','dR2');
P = impose_default_value(P,'dRangerPreprocess_nums_file','all.weird.pairs.nums');
P = impose_default_value(P,'working_subdir','sa_4');
P = impose_default_value(P,'dRangerRun_input_matfile','dRangerRun_input.mat');
P = impose_default_value(P,'draw_heatmaps',false);
P = impose_default_value(P,'error_if_not_all_results_available',true);
%P = impose_default_value(P,'identity_cutoff',0.8);
%P = impose_default_value(P,'cutoff_mode','or');    % 'and'/'or'


%%
%%
fprintf('2009-10-15\n');
fprintf('dRanger_stringent_realign needs to be modified before use:\n');
fprintf('Input and output files need to be changed\n');
keyboard;
%%
%%


if P.cancer_sample, tn = {'normal','tumor'}; else tn = {'sample'}; end

for i=1:length(tn)
  if strcmp(tn{i},'normal') & P.tumor_only, continue; end
  if strcmp(tn{i},'tumor') & P.normal_only, continue; end

  fprintf('PROCESSING %s\n',upper(tn{i}));
  direc = ['/xchip/tcga_scratch/lawrence/' sample '/' tn{i} '_' P.dRangerPreprocess_output_dir_suffix];

  if exist([direc '/' P.dRanger_stringent_realign_output_file])
    fprintf('stringent_realign output file already exists\n');
    continue;
  end

  sadirec = [direc '/' P.working_subdir];

  % load input file
  infile = [direc '/' P.dRangerPreprocess_nums_file];
  if ~exist(infile,'file'), error('%s does not exist',infile); end
  din = dir(infile);
  fprintf('Loading %s\n',infile);
  tic, X = load_struct(infile,repmat('%f',1,13),0); toc
  X = rename_fields(X,colx(1:13),{'rgrp','namenumber','chr1','start1','end1','strand1','qual1',...
                    'chr2','start2','end2','strand2','qual2','flip'});
  nx = slength(X);

  % load batchlist file
  bfile = [sadirec '/batchlist.txt'];
  B = load_struct(bfile,'%f%f%f%f%f%f%f');
  nb = slength(B);

  z = nan(nx,1);
  X.farscore1 = z; X.nearscore1 = z; X.nearpos1 = z;
  X.farscore2 = z; X.nearscore2 = z; X.nearpos2 = z;
  
  % process results files
  all_found = true;
  for b=1:nb
    rfile = [sadirec '/results_batch' num2str(b) '.txt'];
    if exist(rfile,'file')
      fprintf('  Processing file %d/%d.\n',b,nb);
      % format:      (1)       (2)    (3)  (4) (5)     (6)     (7)    (8)   (9)  (10)
      %         whichrecord whichend self? len id prerawscore gaps rawscore pos score
      R = load_matrix(rfile);
      R(:,10) = round(100 * R(:,8) ./ R(:,4));
      % farscore1 (i.e. self realignment)
      idx = find(R(:,2)==1 & R(:,3)==1);
      X.farscore1(R(idx,1)) = R(idx,10);
      % farscore2 (i.e. self realignment)
      idx = find(R(:,2)==2 & R(:,3)==1);
      X.farscore2(R(idx,1)) = R(idx,10);
      % nearscore1 & nearpos1
      idx = find(R(:,2)==1 & R(:,3)==0);
      X.nearscore1(R(idx,1)) = R(idx,10);
      X.nearpos1(R(idx,1)) = R(idx,9) - P.vicinity_backward;
      % nearscore2 & nearpos2
      idx = find(R(:,2)==2 & R(:,3)==0);
      X.nearscore2(R(idx,1)) = R(idx,10);
      X.nearpos2(R(idx,1)) = R(idx,9) - P.vicinity_backward;
    else
      fprintf('  File %s not found.\n',rfile);
      all_found = false;
    end
  end
  if ~all_found
    if P.error_if_not_all_results_available
      error('Not all results files were found!');
    else
      fprintf('Not all results files were found:\n');
      fprintf('Missing results will not count as "discard" votes.\n');
    end
  end

  % draw heatmaps if requested
  if P.draw_heatmaps, subfunction_draw_heatmaps; end

%  % make cutoffs
%  if isempty(P.identity_cutoff) | ~isnumeric(P.identity_cutoff) | ...
%             P.identity_cutoff < 0.5 | P.identity_cutoff > 2.0 %<--- temporary kludge for control experiment
%    error('P.identity_cutoff must be between 0.5 and 1.0');
%  end  
%  if strcmpi(P.cutoff_mode,'or')
%    discard = find(X.score1>=P.identity_cutoff | X.score2>=P.identity_cutoff);
%  elseif strcmpi(P.cutoff_mode,'and')
%    discard = find(X.score1>=P.identity_cutoff & X.score2>=P.identity_cutoff);
%  else
%    error('Unknown setting for P.cutoff_mode: should be "and" or "or".');
%  end
%
%  Xkeep = rmfield(X,{'seq1','seq2','score1','score2','pairmatepos1','pairmatepos2'});
%  Xkeep = reorder_struct(Xkeep, setdiff(1:nx,discard));

  % save file
  outfile = [direc '/' P.dRanger_stringent_realign_output_file];
  fprintf('Saving file %s\n',outfile);
  tic,save(outfile,'X','-v7.3');toc

end % next i  (T/N)



  function subfunction_draw_heatmaps

    fprintf('Making score heatmap\n');
    H = zeros(100,100);
    x = round(100*X.nearscore1); y = round(100*X.nearscore2);
    x(isnan(x)|x<1) = 1; y(isnan(y)|y<1) = 1;
    for idx=1:length(x)
      if ~mod(idx,100000), fprintf('%d/%d\n',idx,length(x)); end
      H(x(idx),y(idx)) = H(x(idx),y(idx)) + 1;
    end
    figure(1);imagesc(H);

    fprintf('Making nearpos heatmap\n');
    Q = zeros(1000,1000);
    x = round(X.nearpos1/10); y = round(X.nearpos2/10);
    x(isnan(x)|x<1) = 1; y(isnan(y)|y<1) = 1;
    for idx=1:length(x)
      if ~mod(idx,100000), fprintf('%d/%d\n',idx,length(x)); end
      Q(x(idx),y(idx)) = Q(x(idx),y(idx)) + 1;
    end
    figure(2);imagesc(Q);

    fprintf('Making nearscore/nearpos heatmap\n');
    nx = slength(X);
    Z = zeros(10000,100);
    x = [X.nearpos1;X.nearpos2];
    y = [X.nearscore1;X.nearscore2];
    x(isnan(x)|x<1) = 1; y(isnan(y)|y<1) = 1; y(y>100)=100;
    for idx=1:length(x)
      if ~mod(idx,100000), fprintf('%d/%d\n',idx,length(x)); end
      Z(x(idx),y(idx)) = Z(x(idx),y(idx)) + 1;
    end
    cap = 2000;figure(3);imagesc(min(cap,Z));ylim([1 500]);
    set(gcf,'position',[ 172    80   836   720]);


    fprintf('Type "return" to continue.\n');
    keyboard
  end  % subfunction


end % main funcion
