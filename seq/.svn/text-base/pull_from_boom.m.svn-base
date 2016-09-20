function [R B S ri bi si] = pull_from_boom(boomdir,chr,st,en,params)
% [R B S ri bi si] = pull_from_boom(boomdir,chr,st,en,params)
%
% "st" and "en" may be vectors of the same length, 
%   in which case results will be concatenated in the order requested
%
% "chr" may not be a vector
%
% R = reads
% B = bases
% S = sequence of reference genome
%
% ri,bi,si = partitions by requested-region (useful in vectorized format)
%            each tells the *last* item belonging to that requested-region
%
% params.quiet -- quiet flag [default = 0]
% params.refdir -- reference directory [default = '/xchip/tcga/gbm/analysis/lawrence/genome/hg18']
%                  (non-default setting currently ignored)
%
% the columns of R are:
%  (1) readgroup (lane)
%  (2) readnumber (name hash)
%  (3) whichpairmate (1 or 2--or 0 for unpaired)
%  (4) readstart (1-based in chr)
%  (5) readend (1-based in chr)
%  (6) strand (0=plus, 1=minus)
%  (7) number of mismatches
%  (8) mapping quality
%  (9) index into B (1-based)
%      Note: index of first base = R(i,9)
%            index of last base  = R(i+1,9)-1 = R(i,9)+R(i,5)-R(i,4)
%  (10) pairmate chr (0=unmapped, -1=unpaired)
%  (11) pairmate start (1-based, -1=unmapped/unpaired)
%  (12) pairmate strand (0=plus, 1=minus, -1=unmapped/unpaired)
%  (13) sum of qualities of non-reference bases
%
% the columns of B are:
%  (1) base    -100=deletion
%              -1=N   1=A  2=C  3=G  4=T   base=reference
%              63=N  65=A 66=C 67=G 68=T   base=non-reference
%  (2) base quality (-100=deletion)
%  (3) read index back to R
%  (4) base position (1-based in chr)
%
% Mike Lawrence and Gaddy Getz 2009-08

if ~exist('params','var'), params=[]; end
params=impose_default_value(params,'quiet',0);
params=impose_default_value(params,'refdir','/xchip/cga1/lawrence/xchip_tcga_gbm_analysis_lawrence/genome/hg18');
params=impose_default_value(params,'maxreadlen',130);

if ~exist(boomdir,'dir'), error('directory %s not found',boomdir); end

if any(size(chr)>1), error('multiple chromosomes not supported'); end
if ~isnumeric(chr) || chr<1 || chr>24, error('chr should be 1-24'); end
chrlen = load_chrlen();
chrlen = chrlen(chr);

if isempty(st), error('pull_from_boom: empty request');end

if ~exist('en','var'), en=st; end
if any(size(st)>1) | any(size(en)>1)   % vector mode
  if length(st)~=length(st(:)) | length(en)~=length(en(:)), error('"st" and "en" can be vectors but not matrices'); end
  if length(st)~=length(en), error('"st" and "en" must be same length'); end
end  

if any(st<1), error('"st" cannot be less than 1'); end
if any(en>chrlen), error('"en" cannot be greater than chrlen'); end

% pull data from files

tt=tic;

% load reference sequence
subfprintf('\nLoading reference sequence...');
st=as_column(st);
en=as_column(en);
S = upper(genome_region(chr,st,en,P.refdir))';
if length(st)>1, S = cat(2,S{:})'; end   % concatenate reference sequences
si = cumsum(en-st+1);

% load up region(s) in index
subfprintf('\nBoom file read:');
subfprintf('\n  Index lookup...');
stem = [boomdir '/chr' num2str(chr) '.'];
fname = [stem 'boomindex'];
d = dir(fname);
if isempty(d), error('boomindex not found!'); end
idxlen = d.bytes/8;
sti = get_block([stem 'boomindex'],'long',min(idxlen-1,max(1,(st-1)-params.maxreadlen)));
eni = get_block([stem 'boomindex'],'long',min(idxlen-1,en-1));

subfprintf('\n  Loading start/end...');
readstart = get_block([stem 'readstart'],'int',sti,eni);
readend = get_block([stem 'readend'],'int',sti,eni);

% trim
subfprintf('\n  Trimming to exact overlap range(s)...');

keepread = true(length(readstart),1);
keepregion = true(length(st),1);
pos1 = 1;
for i=1:length(st)
  pos2 = pos1+eni(i)-sti(i);

  firstok = find(readend(pos1:pos2)>=st(i),1,'first');
  lastok = find(readstart(pos1:pos2)<=en(i),1,'last');

  if isempty(firstok) || isempty(lastok) || firstok>lastok
    % no reads overlap this region
    keepread(pos1:pos2) = false;
    keepregion(i) = false;
  else
    keepread(pos1:pos1+firstok-2) = false;
    keepread(pos1+lastok:pos2) = false;
    sti(i)=sti(i)+(firstok-1);
    eni(i)=eni(i)-(pos2-pos1+1-lastok);
  end

  pos1 = pos2+1;
end

readstart = readstart(keepread);
readend = readend(keepread);
sti = sti(keepregion);
eni = eni(keepregion);
st = st(keepregion);
en = en(keepregion);

if isempty(st)    % no data
  R=[];
  B=[];
  ri=[];
  bi=[];
  return
end

subfprintf('\n  Allocating read storage...');
sz = eni-sti+1;
cumsz = cumsum(sz);
R = zeros(cumsz(end),13);

subfprintf('\n  Loading other read data...');
R(:,1) = get_block([stem 'readgroup'],'short',sti,eni);
R(:,2) = get_block([stem 'namenumber'],'int',sti,eni);
R(:,3) = get_block([stem 'whichpairmate'],'byte',sti,eni);
R(:,4) = readstart;
R(:,5) = readend;
R(:,6) = get_block([stem 'readstrand'],'byte',sti,eni);
R(:,7) = get_block([stem 'nmismatches'],'byte',sti,eni);
R(:,8) = get_block([stem 'readmapqual'],'byte',sti,eni);
R(:,9) = get_block([stem 'baseindex'],'long',sti,eni);
R(:,10) = get_block([stem 'pairmatechr'],'byte',sti,eni);
R(:,11) = get_block([stem 'pairmatestart'],'int',sti,eni);
R(:,12) = get_block([stem 'pairmatestrand'],'byte',sti,eni);
R(:,13) = get_block([stem 'nonrefsumq'],'short',sti,eni);

subfprintf('\n  Allocating base storage...');
if length(st)>1
  bsti = R([1;cumsz(1:end-1)+1],9);
  beni = R(cumsz,9)+R(cumsz,5)-R(cumsz,4);
else
  bsti = R(1,9);
  beni = R(end,9)+R(end,5)-R(end,4);
end

bsz = beni-bsti+1;
bcumsz = cumsum(bsz);
B = zeros(bcumsz(end),4);

subfprintf('\n  Loading base data...');
B(:,1) = get_block([stem 'base'],'byte',bsti,beni);
B(:,2) = get_block([stem 'basequal'],'byte',bsti,beni);

% adjust baseindex to the subselection of "base" and "basequal"

pos1=1;
offset=0;
for i=1:length(st)
  pos2 = cumsz(i);
  R(pos1:pos2,9) = R(pos1:pos2,9) - (bsti(i)-1) + offset;
  pos1 = pos2+1;
  offset = bcumsz(i);
end

subfprintf('\n  Total time for memory allocation / disk read = %.2f sec\n',toc(tt));

% convert to pull_from_bam conventions
subfprintf('Converting to pull_from_bam conventions...');
ump = (R(:,10)==-2);   % unmapped pairmate
R(ump,10) = 0;
R(ump,11:12) = -1;

% add additional B columns
subfprintf('\nComputing additional B columns...');
B(:,4) = 1;
B(R(:,9),3)=1;
B(R(:,9),4)=R(:,4)-[0; R(1:end-1,5)];
B(:,[3 4])=cumsum(B(:,[3 4]),1);

% partitions between requested-regions   (including empty regions)
tmp = [0;cumsz];
ri = tmp(cumsum(keepregion)+1);
tmp = [0;bcumsz];
bi = tmp(cumsum(keepregion)+1);

subfprintf('\nDone: total time for pull_from_boom = %.2f sec\n',toc(tt));



  function subfprintf(str,varargin)
    if ~(params.quiet), fprintf(str,varargin{:}); end
  end



end % main function







