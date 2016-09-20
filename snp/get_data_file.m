function S=get_data_file(fname,birdseed_format,floor_val,logdata,overwrite_mat)
if ~exist('logdata','var') || isempty(logdata)
  logdata=1;
end
if ~exist('birdseed_format','var') || isempty(birdseed_format)
  birdseed_format=0;
end
if ~exist('overwrite_mat','var') || isempty(overwrite_mat)
  overwrite_mat=0;
end

[tmp,fname_noext]=file_ext(fname);
mat_fname=[ fname_noext '.mat'];
if exist(mat_fname,'file') && ~overwrite_mat
  verbose(['Reading ' fname_noext '.mat'],10);
  S=load_D(mat_fname);
else
  verbose(['Reading ' fname],10);
  if birdseed_format
    verbose('Assuming birdseed format',30);
    S=read_modelled_data_file(fname,-1,-1,0,0,0,0,0,5e5);     
  else
    S=read_modelled_data_file(fname,-1,-1,1,0,0,0,0,5e5); 
    if ~exist('floor_val','var')
      floor_val=0.001;
    end
    verbose(['Setting floor_val to ' num2str(floor_val)],30);
    S.dat(S.dat<floor_val)=floor_val;
    if  logdata
      verbose('Taking log2(x)-1',30);
      S.dat=log2(S.dat)-1;
    end
  end
  verbose([ 'Saving ' mat_fname],10);
  try
    save_D(mat_fname,S);
  catch
    disp('Could not write mat file');
  end
end

