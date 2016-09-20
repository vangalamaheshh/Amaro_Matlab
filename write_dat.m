function write_dat(fname,D,file_type)

switch(file_type)
 case 'ctwc' 
  write_eisen_dat(fname,strvcat(D.gacc),strvcat(D.gdesc), ...
                  strvcat(getfields(D.scans,'name')),'AFFY_U133',D.dat,...
                  [],[],0);
  
 otherwise
  disp('NO SUCH TYPE');
end


