function C = load_capture_stats(C,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'skip_bam_stats',true);

try
C.sample.tlanes = nan(C.ns,1);
C.sample.nlanes = nan(C.ns,1);
C.sample.treads = nan(C.ns,1);
C.sample.nreads = nan(C.ns,1);
C.sample.treads_ontarget = nan(C.ns,1);
C.sample.nreads_ontarget = nan(C.ns,1);
for i=1:C.ns
  [C.sample.tlanes(i) C.sample.nlanes(i)] = get_lanecounts(C.sample.dir{i});
  if ~P.skip_bam_stats
    [t n] = get_bam_stats(C.sample.dir{i});
    if ~isempty(t), C.sample.treads(i) = t.total - t.non_pf - t.duplicates; end
    if ~isempty(n), C.sample.nreads(i) = n.total - n.non_pf - n.duplicates; end
  else
    %%%%%
  end
  if isfield(C.file,'globalstats')
    fname = C.file.globalstats;
  else
    fname = ['coverage_by_' lower(C.captype) '.stats'];
  end
  fullname = ['/xchip/tcga_scratch/lawrence/' C.sample.dir{i} '/' fname];
  if exist(fullname,'file')
    X = load_struct(fullname,'%f%s%f%f%f%f%f%f%f%f');
    C.sample.treads_ontarget(i) = X.tseqreads(X.categ_num==1);    
    C.sample.nreads_ontarget(i) = X.nseqreads(X.categ_num==1);
  end
end

catch me, excuse(me); end
