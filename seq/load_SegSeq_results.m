function X = load_SegSeq_results(infile,P)
if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'log_base',2);

X = load_struct(infile);
f = fieldnames(X);
if length(f)==6
  newflds  = {'sample';'chromosome';'start';'end';'nprobes';'copyratio'};
  X = rename_fields(X,f,newflds);
elseif length(f)==5
  newflds  = {'chromosome';'start';'end';'nprobes';'copyratio'};
  X = rename_fields(X,f,newflds);
else
  fprintf('unknown segfile format: attempting to interpret\n');
  if isfield(X,'chrom'), X = rename_field(X,'chrom','chromosome'); end
  if isfield(X,'locstart'), X = rename_field(X,'locstart','start'); end
  if isfield(X,'locend'), X = rename_field(X,'locend','end'); end
  if isfield(X,'segmean'), X = rename_field(X,'segmean','copyratio'); end
  demand_fields(X,{'start','end','chromosome','copyratio'});
end

X.chr = convert_chr(X.chromosome);
X = make_numeric(X,{'start','end'});
X.len = X.end - X.start + 1;
X.chromosome = regexprep(X.chromosome,'^23$','X');
X.chromosome = regexprep(X.chromosome,'^24$','Y');

X = make_numeric(X,'copyratio');
if any(X.copyratio<=0) || ~any(X.copyratio>=1)
  X = rename_field(X,'copyratio','log_copyratio');
  X.copyratio = P.log_base.^X.log_copyratio;
end
