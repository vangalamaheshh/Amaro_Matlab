function create_summed_FWB_coverage_file_2(fwbfiles,region_list,fwbOutName,width,shared_fwifile)
% create_summed_FWB_coverage_file_2(fwbfiles,target_list,fwbOutName)
%
% fwbfiles = list of fwb's to sum: can be any length or width
% region_list = filename of region list to sum across; will be copied to fwiOutName
% fwbOutName = filename of fwb to write
% width = width of fwb to write: can be 8, 16, or 32.  default=16
%
% Note: Works with any arbitrary collection of FWB files,
%       i.e. they don't have to share the same target_list or have width=1

if exist(fwbOutName,'file')
  fprintf('output file already exists.\n');
  return
end

if ~exist('width','var'), width=16; end
if ~ismember(width,[8 16 32]), error('unsupported width'); end
fmt = ['uint' num2str(width)];
nf = length(fwbfiles);
if nf>(2.^width), error('width is too small to sum that many files'); end
demand_files(fwbfiles);
fwiOutName = regexprep(fwbOutName,'.fwb$','.fwi');
if strcmpi(fwiOutName,fwbOutName), error('fwbOutName should end in ".fwb"'); end

% load region list
T = load_struct_noheader(region_list);
T = rename_fields(T,{'col1','col2','col3'},{'chr','start','end'});
T = make_numeric(T,{'start','end'});
T.chr = convert_chr(T.chr);
nt = slength(T);
if any(isnan(T.chr)), error('invalid chr(s) in region_list'); end
if any(isnan(T.start)), error('invalid start(s) in region_list'); end
if any(isnan(T.end)), error('invalid end(s) in region_list'); end
if any(T.end<T.start), error('end<start in region_list'); end

% open input files
F = cell(nf,1);
if exist('shared_fwifile','var')
  for f=1:nf, F{f} = org.broadinstitute.cga.tools.seq.FixedWidthBinary(fwbfiles{f},shared_fwifile); end
else
  for f=1:nf, F{f} = org.broadinstitute.cga.tools.seq.FixedWidthBinary(fwbfiles{f}); end
end

% open output file
out = fopen(fwbOutName,'w');
endianness = 'b';

% sum it
step = ceil(nt/100);
fprintf('Summing over %d targets: ',nt);
for i=1:nt, if ~mod(i,step), fprintf('%d/%d ',i,nt); end
  cov = zeros(T.end(i)-T.start(i)+1,1);
  for f=1:nf
    cov = cov + max(0,double(F{f}.get(T.chr(i),T.start(i),T.end(i))));
  end
  fwrite(out,cov,fmt,endianness);
% (noticed the missing "endianness" on 12/1/12: earlier files written with this function may be wrong)
end, fprintf('\n');

% finish
fclose(out);
save_struct_noheader(T,fwiOutName);

