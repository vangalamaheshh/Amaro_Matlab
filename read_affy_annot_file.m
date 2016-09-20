function annot=read_affy_annot_file(fname)

unixstr=['cat ' fname ' | tr "|" "\\" | sed ''s/","/|/g'' | tr -d ''"'' > ' fname '.txt' ];
disp([ 'UNIX: ' unixstr ]);
unix(unixstr);

d=read_dlm_file([fname '.txt'],'|');

if (0)
  f=fopen([fname '.txt']);

  l=fgetl(f);
  d=dlmsep(l,'|');
  
  form=repmat('%s',1,length(d));
  %form=form(1:(end-1));
  c=textscan(f,form,'delimiter','|','bufsize',1000000);
  
  fclose(f);
end

annot.headers=d{1};
annot.data=cat(1,d{2:end});

annot=parse_affy_annot(annot);
