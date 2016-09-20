function create_summed_FWB_coverage_file(fwbfiles,target_list,fwbOutName,P)
% create_summed_FWB_coverage_file(fwbfiles,target_list,fwbOutName)
%
% fwbfiles = list of fwb's to sum: must have width=1
% target_list = filename of target list (gene chr st en) that they were generated on, i.e. their index(2,3,4)
% fwbOutName = filename of fwb to write: will have width=16
%
% Note: Only works with input files that all share the same region list and have width=1
%       To sum FWB an arbitrary collection of FWB files, use create_summed_FWB_coverage_file_2

if ~exist('P','var'), P=[]; end

if exist(fwbOutName)
  fprintf('create_summed_FWB_coverage_file: output file already exists\n');
  return
end

inWidth = 1;
outWidth = 16; % can handle up to ~32,000 patients
outFormat = ['uint' num2str(outWidth)];
outEndianness = 'b';

if ~iscell(fwbfiles), error('fwbfiles should be a cell array of filenames'); end
demand_files(fwbfiles);
fwiOutName = regexprep(fwbOutName,'.fwb$','.fwi');
if strcmpi(fwbOutName,fwiOutName), error('fwbOutName should end with .fwb'); end
fprintf('Creating summed FWB coverage file: %s\n',fwbOutName);

demand_file(target_list);
T = load_target_file(target_list,P);
I = keep_fields(T,{'chr','start','end'});
save_struct_noheader(I,fwiOutName);

bp = sum(T.len);
in_len = ceil(bp/8);
out_len = bp;

fprintf('Opening files: ');
nf = length(fwbfiles);
f = cell(nf,1);
for fi=1:nf, if ~mod(fi,10), fprintf('%d/%d ',fi,nf); end
  f{fi} = fopen(fwbfiles{fi},'rb');
  fseek(f{fi}, 0, 'eof');
  filelen = ftell(f{fi});
  if filelen~=in_len
    error('fwbfile(s) have incorrect size!  index mismatch?');
  end
  fseek(f{fi}, 0, 'bof');
end,fprintf('\n');

fprintf('Summing and writing output file:\n');
out = fopen(fwbOutName,'wb');
chunksize = 100000;
aa=tic;
data = zeros(chunksize, nf, 'uint8');
for pos=1:chunksize:in_len, fprintf('%d/%d ',ceil(pos/chunksize),ceil(filelen/chunksize));toc(aa);
  if pos+chunksize > in_len
    last_chunk = true;
    bytes_to_read = in_len - pos + 1;
    data = zeros(bytes_to_read,nf,'uint8');
  else
    last_chunk = false;
    bytes_to_read = chunksize;
  end
  for fi=1:nf, data(:,fi) = fread(f{fi},bytes_to_read,'uint8=>uint8'); end
  tot = sum_FWB_1bit_to_16bit(data);   % mex C function
  if last_chunk
    % trim last chunk if necessary (to remove the padding bits from the width=1 representation)
    bytes_to_trim = mod(8-mod(bp,8),8);
    if bytes_to_trim>0
      tot = tot(1:end-bytes_to_trim);
    end
  end
  fwrite(out,tot,outFormat,outEndianness);
end, fprintf('\n');
fclose(out);

fprintf('Closing input files\n');
for fi=1:nf, fclose(f{fi}); end

fprintf('Done\n');

