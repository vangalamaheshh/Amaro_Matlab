function ensure_writeable(fname)

try
  [path name ext] = fileparts(fname);
  if ~isempty(path), ensure_dir_exists(path); end
  testfile = [fname '.' rand() '.test'];
  save_textfile('test',testfile);
  delete(testfile);
catch me
  error('%s is not writeable',fname);
end
