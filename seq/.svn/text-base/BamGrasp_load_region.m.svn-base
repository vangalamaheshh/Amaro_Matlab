function [R B S ri bi si rname] = BamGrasp_load_region(x,chr,st,en,params)
% "x" is handle to BamGrasp object which has already had openFile called

if nargin==4 && isstruct(en)
  % ALTERNATE CALLING SYNTAX: BamGrasp_load_region(x,chr,st,params)
  params = en;
  clear en;
end

if ~exist('params','var'), params=[]; end
params=impose_default_value(params,'quiet',1);
params=impose_default_value(params,'maxreads',[]);
params=impose_default_value(params,'include_unmapped_reads',true);
params=impose_default_value(params,'include_duplicate_reads',false);

if ~exist('en','var'), en=st; end
if any(size(st)>1) || any(size(en)>1)   % vector mode
  error('pull_from_bam: vectorized mode not yet implemented');
end

if ~isnumeric(chr), chr = convert_chr(chr); end
if ~isnumeric(st), st = str2double(st); end
if ~isnumeric(en), en = str2double(en); end

if ~isempty(params.maxreads), x.set_maxReads(params.maxreads); end

if params.include_unmapped_reads
  subfprintf('Including unmapped reads\n');
  x.set_unmapped_reads_included();
else
  subfprintf('Excluding unmapped reads\n');
  x.set_unmapped_reads_excluded();
end

if params.include_duplicate_reads
  subfprintf('Including duplicate reads\n');
  x.set_duplicate_reads_included();
else
  subfprintf('Excluding duplicate reads\n');
  x.set_duplicate_reads_excluded();
end

subfprintf('Bam file read:\n');
tt=tic;
x.loadRegion(chr,st,en);
subfprintf('  Total time for memory allocation / disk read = %.2f sec\n',toc(tt));

subfprintf('Java->MATLAB accessor functions...');

if ~params.include_aux_cols
  R = double([x.getReadgroup x.getNamenumber x.getWhichpairmate ...
              x.getReadstart x.getReadend x.getReadstrand ...
              x.getNmismatches x.getReadmapqual ...
              x.getBaseindex+1 x.getPairmatechr x.getPairmatestart x.getPairmatestrand]);
else
  R = double([x.getReadgroup x.getNamenumber x.getWhichpairmate ...
              x.getReadstart x.getReadend x.getReadstrand ...
              x.getNmismatches x.getReadmapqual ...
              x.getBaseindex+1 x.getPairmatechr x.getPairmatestart x.getPairmatestrand ...
              x.getInsertsize x.getIsdup]);
end

B = double([x.getBase x.getBasequal x.getBaseReadIndex x.getBasePosition]);
if isempty(B), B = zeros(0,4); end
S = x.getReference;
rname = x.getReadname;

% (NOW DONE IN JAVA)
%subfprintf('\nComputing additional B columns...');
%tmp=[zeros(size(B,1),1) ones(size(B,1),1)];
%tmp(R(:,9),1)=1;
%tmp(R(:,9),2)=R(:,4)-[0; R(1:end-1,5)];
%B=[B cumsum(tmp,1)];

% partitions (trivial in non-vectorized mode)
ri = size(R,1);
bi = size(B,1);
si = length(S);

subfprintf('\nDone: total time for pull_from_bam = %.2f sec\n',toc(tt));

  function subfprintf(str,varargin)
    if ~(params.quiet), fprintf(str,varargin{:}); end
  end

end

