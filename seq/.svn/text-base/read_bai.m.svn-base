function B = read_bai(fname)
B = [];
if ~exist(fname,'file'), error('File not found!'); end
f = fopen(fname,'rb');
magic_num = fread(f,4,'uint8=>char');
if ~strcmp(magic_num,['BAI' 1]'), error('Not a baifile!'); end

B.n_ref = double(fread(f,1,'uint8=>int32'));

for r=1:B.n_ref
  B.n_bin(r) = double(fread(f,1,'uint8=>int32'));
  for b=1:B.n_bin(r)
    B.bin(r,b) = double(fread(f,1,'uint8=>uint32'));
    B.n_chunk(r,b) = double(fread(f,1,'uint8=>int32'));
    for c=1:B.n_chunk(r,b)
      B.chunk_beg(r,b,c) = double(fread(f,1,'uint8=>uint64'));
      B.chunk_end(r,b,c) = double(fread(f,1,'uint8=>uint64'));
    end
  end
  B.n_intv(r) = double(fread(f,1,'uint8=>int32'));
  for i=1:B.n_intv(r)
    B.ioffset(r,i) = double(fread(f,1,'uint8=>uint64'));
  end
end

fclose(f);










