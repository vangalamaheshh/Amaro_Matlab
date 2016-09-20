function D=append_D_sup(D,supid,addtoacc,addtodesc,rc)

if ~exist('rc','var')
  rc='cols';
end


if is_col(rc)
  sacc=cellstr(D.supacc);
  sacc{supid}=[ sacc{supid} addtoacc ];
  D.supacc=strvcat(sacc);
  sdesc=cellstr(D.supdesc);
  sdesc{supid}=[ sdesc{supid} addtoacc ];
  D.supdesc=strvcat(sdesc);
else
  sacc=cellstr(D.gsupacc);
  sacc{supid}=[ sacc{supid} addtoacc ];
  D.gsupacc=strvcat(sacc);
  sdesc=cellstr(D.gsupdesc);
  sdesc{supid}=[ sdesc{supid} addtoacc ];
  D.gsupdesc=strvcat(sdesc);  
end

