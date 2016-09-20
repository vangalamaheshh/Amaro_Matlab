function collapse_RegionCovPerLane_to_genes(samples,cov_filestem,genecov_filestem)
% Mike Lawrence 2009-07-01

if ~iscell(samples), samples = { samples }; end

basedir = '/xchip/tcga_scratch/lawrence';
for s=1:length(samples)
  fprintf('Processing file %d/%d',s,length(samples));
  sample = samples{s};
  direc = [basedir '/' sample];
  for i=1:2
    if i==1, tn = 'tumor'; else tn = 'normal'; end
    fprintf(' %s',tn);
    infile = [direc '/' cov_filestem '_' tn '.txt'];
    outfile = [direc '/' genecov_filestem '_' tn '.txt'];

    if exist(outfile), continue; end

    % load coverage-by-targets file
    if ~exist(infile), error('%s not found!',infile); end
    X = load_struct_specify_string_cols(infile,[1],0);
    nlanes = length(fields(X))-4;
    X.len = X.col4-X.col3+1;
    nx = slength(X);
    Ct = nan(nx,nlanes);
    for lane=1:nlanes
      Ct(:,lane) = getfield(X,['col' num2str(4+lane)]);
    end

    % collapse to genes
    G = [];
    [G.name ui uj] = unique(X.col1);
    ng = length(G.name);
    Cg = nan(ng,nlanes);
    for g=1:ng, idx = find(uj==g);
      chr = unique(X.col2(idx));
      if length(chr)~=1, error('Gene %s is on %d different chromosomes',G.name{g},length(chr)); end
      G.chr(g,1) = chr;
      G.min(g,1) = min(X.col3(idx));
      G.max(g,1) = max(X.col4(idx));
      G.len(g,1) = sum(X.len(idx));
      Cg(g,:) = sum(Ct(idx,:),1); 
    end
    for lane=1:nlanes
      G = setfield(G,['lane' num2str(lane-1)],Cg(:,lane));
    end

    % save collapsed-to-genes file
    save_struct(G,outfile);

  end  % next T/N
  fprintf('\n');

end  % next sample

fprintf(' Done.\n');

