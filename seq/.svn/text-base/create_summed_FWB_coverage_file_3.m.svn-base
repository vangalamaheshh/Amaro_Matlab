function create_summed_FWB_coverage_file_3(infiles,inwidth,outfile,outwidth)
% create_summed_FWB_coverage_file_3(infiles,target_list,outfile)
%
% infiles = list of fwb's to sum: must all be on the same index and have the same width
% inwidth = width of input infiles: can be 8, 16, or 32.
% outfile = filename of fwb to write
% outwidth = width of fwb to write: can be 8, 16, or 32.
%
% create_summed_FWB_coverage_file_3 is good for summing *many* width>=8 FWBs that are all on the same index

if exist(outfile,'file')
  fprintf('output file already exists.\n');
  return
end

if ~ismember(inwidth,[8 16 32]), error('unsupported inwidth'); end
if ~ismember(outwidth,[8 16 32]), error('unsupported outwidth'); end
infmt = ['uint' num2str(inwidth)];
outfmt = ['uint' num2str(outwidth)];
endianness = 'b';

% open input files
nf = length(infiles);
F = cell(nf,1);
fprintf('Opening files: ');
for i=1:nf, if ~mod(i,100), fprintf('%d/%d ',i,nf); end
  F{i} = fopen(infiles{i},'rb');
end, fprintf('\n');

% get file size
sz = filesize(infiles{1}) / (inwidth/8);

% open output file
out = fopen(outfile,'wb');
outmax = (2.^(outwidth-1))-1;

% sum it
chunksize = 1e6;
fprintf('Summing by %dbp chunks:\n',chunksize);
tic
for st=1:chunksize:sz-1, fprintf('%d/%d ',ceil(st/chunksize),ceil(sz/chunksize));
  en = min(sz,st+chunksize-1);
  len = en-st+1;
  for f=1:nf, if ~mod(f,1000), if f==1000, fprintf('(')); end; fprintf('%d/%d ',f,nf); end
    x = double(fread(F{f},len,[infmt '=>' infmt],endianness));
    if f==1, tot=x; else tot=tot+x; end
  end, if f>=1000, fprintf(') '); end
  over = (tot>outmax);
  if any(over)
    fprintf('\t:Note %d/%d entries were over width %d max = %d\n',sum(over),len,outwidth,outmax);
    tot(over)=outmax;
  end
  fwrite(out,tot,outfmt,endianness);
  toc
end

% finish
fclose(out);
for i=1:nf, fclose(F{i}); end


