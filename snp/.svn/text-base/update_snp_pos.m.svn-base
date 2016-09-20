function C=update_snp_pos(C)

f=fopen('~/projects/snp/data/SNP/affy_id.snp.pos.txt');
s=textscan(f,'%d %s %s %d');
fclose(f);

if range(range(strvcat(C.marker)-strvcat(s{2})))>0
  error('files don''t match');
else
  C.pos=double(s{4})/1e6;
end
