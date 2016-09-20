function [tbl,D]=read_qs_file(fname,D);

tbl=read_table(fname,['%s' repmat('%f',1,9)],char(9),1);

if exist('D','var') && ~isempty(D)
  [Mt,m1,m2]=match_string_sets_hash(tbl.dat{1},D.sdesc);
  error('Finish me');
end
