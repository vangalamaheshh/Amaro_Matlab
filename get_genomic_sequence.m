function s=get_genomic_sequence(path,chr,st,en)
%s=get_genomic_sequence('/xchip/data01/gadgetz/vogelstein/hg17/',10,1001000:1002000);
if ~exist('en','var')
  en=max(st);
  st=min(st);
end

m=memmapfile([ path 'chr' num2chromosome(chr) '.bin']);
s=char(m.data(st:en))';
