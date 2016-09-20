function check_file_exists_and_delete(fname)
cf_='[lsf::wait]';
if exist(fname,'file')
  verbose([cf_ 'deleting ' fname]);
  delete(fname);
end
