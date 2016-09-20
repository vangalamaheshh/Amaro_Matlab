function [terr cov] = tabulate_coverage(covfwb,categfwb,regionlist,outbinfile_terr,outbinfile_cov)
% [terr cov] = tabulate_coverage(covfwb,categfwb,regionlist[,outbinfile_terr,outbinfile_cov])
%
% covfwb = fwb with coverage to tabulate
% categfwb = fwb with categories to tabulate coverage by
% regionlist = list of regions to get coverage on (no header, chr start end)
% outbinfile_terr = binary file that will be written with the territory (per region per category)
% outbinfile_cov = binary file that will be written with the coverage (per region per category)

if nargin>=4
  if exist(outbinfile_terr,'file')
    if nargin==4 || exist(outbinfile_cov,'file')
      fprintf('output files already exist.\n');
      return
end,end,end

% load region list
T = load_struct_noheader(regionlist);
nf = length(fieldnames(T));
if nf==3
  T = rename_fields(T,{'col1','col2','col3'},{'chr','start','end'});
elseif nf==4
  T = rename_fields(T,{'col1','col2','col3','col4'},{'gene','chr','start','end'});
else
  error('not sure how to parse target list');
end
T = make_numeric(T,{'start','end'});
T.chr = convert_chr(T.chr);
nt = slength(T);
if any(isnan(T.chr)), error('invalid chr(s) in region_list'); end
if any(isnan(T.start)), error('invalid start(s) in region_list'); end
if any(isnan(T.end)), error('invalid end(s) in region_list'); end
if any(T.end<T.start), error('end<start in region_list'); end

% get number of categories
fn = regexprep(categfwb,'all.fwb','categs.txt');
demand_file(fn);
cats = load_struct(fn);
nk = slength(cats);

% open FWBs
C = org.broadinstitute.cga.tools.seq.FixedWidthBinary(covfwb);
K = org.broadinstitute.cga.tools.seq.FixedWidthBinary(categfwb);
cov = zeros(nt,nk);
terr = zeros(nt,nk);

% load it
step = ceil(nt/100);
for i=1:nt, if ~mod(i,step), fprintf('%d/%d ',i,nt); end
  k = double(K.get(T.chr(i),T.start(i),T.end(i)));
  c = double(C.get(T.chr(i),T.start(i),T.end(i)));
  terr(i,:) = histc(k,1:nk);
  mc = max(c); if mc < 1, continue; end
  h = hist2d_fast(c,k,1,mc,1,nk);
  h = bsxfun(@times,h,(1:mc)');    
  cov(i,:) = sum(h,1);
end, fprintf('\n');


keyboard
% save
if nargout>=4
  out = fopen(outbinfile_terr,'w');
  fwrite(out,terr,'uint16');
  fclose(out);
end
if nargin>=5
  out = fopen(outbinfile_cov,'w');
  fwrite(out,cov,'uint16');
  fclose(out);
end


