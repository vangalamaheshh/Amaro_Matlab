function write_D_ctwc(D,fname,dbname)

write_eisen_dat([fname '.dat'],strvcat(fill_empty(D.gacc,'EMPTY')),strvcat(fill_empty(D.gdesc,'EMPTY')),D.sdesc,dbname, ...
                D.dat,[],[],0);

if isfield(D,'supdat')
  write_eisen_dat([fname '.supdat'],strvcat(D.supacc),strvcat(D.supdesc),D.sdesc,'LABELS', ...
                  D.supdat,[],[],1);
end

