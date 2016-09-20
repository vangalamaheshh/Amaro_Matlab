function pon = get_pon(chr,pos,ponfile)

%       #bams
% pon3  2966   old FH normals
% pon4  4214   new FH normals
% pon5   445   external normals (mostly TCGA from Baylor/WashU)
% pon6   376   FH WGS and 15 ExomePlus normals
% pon7  4513   (Petar) = 395 ExomePlus normals + pon4
% pon8   342   PCAWG342 normals
% pon9  8334   TCGA FH normals

if ~exist('ponfile','var')
  ponfile = '/cga/fh/pancan_data/pon/pon9/final_summed_tokens.hist.bin';
end

if ~isnumeric(chr), chr = convert_chr(chr); end
if ~isnumeric(pos), pos = str2double(pos); end
if length(chr) ~= length(pos), error('chr and pos should be same length'); end
n = length(chr);

chrlen = get_chrlen('hg19');
sz = get_filesize(ponfile);
if sz ~= 16*sum(chrlen)
  error('ponfile is not correct length for hg19');
end

coord = nan(n,1);   % 0-based
offset = 0;
for c=1:24
  idx = find(chr==c & pos>=1 & pos<=chrlen(c));
  coord(idx) = pos(idx)-1+offset;
  offset = offset + chrlen(c);
end

f = fopen(ponfile,'rb');
bof = -1;   % beginning of file

pon = nan(n,8);

[tmp ord] = sort(coord);
tt = tic;
for j=1:n, i = ord(j);
  if ~mod(j,1e5), fprintf('%d/%d ',j,n); toc(tt); end
  if isnan(coord(i)), continue; end
  fseek(f,16*coord(i),bof);
  pon(i,:) = fread(f,8,'uint16');
end
 
fclose(f);
