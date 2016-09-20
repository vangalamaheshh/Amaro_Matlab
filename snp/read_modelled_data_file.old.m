function M=read_modelled_data_file(fname)

[f,fid]=read_dlm_file(fname,char(9),1);
form=[ '%s' repmat('%f%*s',1,floor((length(f{1})-2)/2)) '%*s'];
keep_reading=1;
M_dat_full={};
while (keep_reading)
  M_dat=textscan(fid,form,10000,'bufSize',1000000);
  if length(M_dat{1})<10000
    keep_reading=0;
  end
  disp(['read ' num2str(length(M_dat{1}))]);
  if isempty(M_dat_full)
    M_dat_full=M_dat;
  else
    for j=1:length(M_dat)
      M_dat_full{j}=[M_dat_full{j}; M_dat{j}] ;
    end
  end
end
M.marker=M_dat_full{1};
M.dat=cat(2,M_dat_full{2:end}); % we have %*s on the calls 
M.sdesc=f{1}(2:2:(end-1));
fclose(fid);

