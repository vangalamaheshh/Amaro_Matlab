function pos2 = exonic_liftover(varargin)
% pos2 = exonic_liftover(chr1,pos1,build1,build2)

pos2 = liftover(varargin{:});
return





%%%%%% ORIGINAL METHOD BELOW

if nargin==2
  build1=18; build2=19;
  fprintf('Assuming you want hg18->hg19\n');
elseif nargin~=4
  fprintf('Need chr1,pos1,build1,build2\n');
end

if ~isnumeric(chr1), chr1 = convert_chr(chr1); end
if ~isnumeric(pos1), pos1 = str2double(pos1); end

if length(chr1)~=length(pos1), error('length(chr)~=length(pos)'); end

b1 = interpret_build(build1);
b2 = interpret_build(build2);
if ~((b1==18 && b2==19) || (b1==19 && b2==18))
  error('Unsupported combination of builds\n');
end

e18 = load_target_file('/cga/tcga-gsc/home/lawrence/capture/Refseq_exons_good_20101221_hg18.txt');
e19 = load_target_file('/cga/tcga-gsc/home/lawrence/capture/Refseq_exons_good_20101221_hg19.txt');

if b1==18 && b2==19
  fr = e18; to = e19;
elseif b1==19 && b2==18
  fr = e19; to = e18;
else
  error('wha?');
end

flank_trust_size = 200;   % how far outward from the exon do we assume the mapping holds true?
fr.start = fr.start - flank_trust_size;
fr.end = fr.end + flank_trust_size;
to.start = to.start - flank_trust_size;
to.end = to.end + flank_trust_size;

pos2 = nan(length(chr1),1);

for chr=1:24, fprintf('chr%d ',chr);
  idx1 = find(chr1==chr);
  idx2 = find(fr.chr==chr);
  for j=1:length(idx1), i=idx1(j);
    zz = idx2(fr.start(idx2)<=pos1(i) & fr.end(idx2)>=pos1(i));
    if length(zz)>1, zz=zz(1); end
    if ~isempty(zz)
      pos2(i) = to.start(zz) - fr.start(zz) + pos1(i);
    end
  end
end, fprintf('\n');


