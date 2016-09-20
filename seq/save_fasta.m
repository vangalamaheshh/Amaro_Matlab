function save_fasta(S,filename)
%
% saves the given FASTA struct to the specified filename
%
% struct should have fields 'header' and 'seq'
%
% Mike Lawrence 2008-06-24
%
require_fields(S,{'header';'seq'});
if length(S.header) ~= length(S.seq), error('header and seq must be same length'); end
ns = length(S.header);

O=cell(5,1);
O{1} = repmat({'>'},ns,1);
O{2} = S.header;
O{3} = repmat({char(10)},ns,1);
O{4} = S.seq;
O{5} = O{3};
O = strcat(O{:});   % collapse to lines
O = [O{:}];         % collapse to file

out=fopen(filename,'wt');
fwrite(out,O);
fclose(out);
