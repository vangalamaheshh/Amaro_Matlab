function result = get_block(filename,type,st,en,outtype)
% get_block(filename,type,st[,en])
%
% "st" and "en" may be vectors of the same length,
%  in which case results will be concatenated in order requested
%
% filename = file to pull from
%
% type = 'byte'          = signed int8
%        'short'         = signed int16 big-endian
%        'int'           = signed int32 big-endian
%        'long'          = signed int64 big-endian
%        'short-little'  = signed int16 little-endian
%        'int-little'    = signed int32 little-endian
%        'long-little''  = signed int64 little-endian
%
% st = 0-based, in units of "type"
% en = if missing, en=st
%
% returns results as "double", regardless of "type"
%
% Mike Lawrence 2009-08-25

if isempty(st), result=[]; return; end

if ~exist('en','var'), en=st; end
if any(size(st)>1) | any(size(en)>1)   % vector mode
  if length(st)~=length(st(:)) | length(en)~=length(en(:)), error('"st" and "en" can be vectors but not matrices'); end
  if length(st)~=length(en), error('"st" and "en" must be same length'); end
end

if ~ischar(type), error('"type" should be a string'); end
if strcmp(type,'byte'), format='int8'; width=1; endian='l';
elseif strcmp(type,'short'), format='int16'; width=2; endian='b';
elseif strcmp(type,'int'), format='int32'; width=4; endian='b';
elseif strcmp(type,'long'), format='int64'; width=8; endian='b';
elseif strcmp(type,'short-little'), format='int16'; width=2; endian='l';
elseif strcmp(type,'int-little'), format='int32'; width=4; endian='l';
elseif strcmp(type,'long-little'), format='int64'; width=8; endian='l';
else error('Unknown type %s',type);
end

if ~exist(filename,'file'), error('%s not found', filename); end
d = dir(filename);
len = d.bytes/width;
if any(st<0), error('get_block: st<0'); end
if any(st==1), fprintf('WARNING:  get_block expects zero-based coordinates\n'); end
if any(st>=len), error('get_block: st>=filelength'); end
if any(en>=len), error('get_block: en>=filelength'); end

f = fopen(filename);
if f==-1, error('Can''t open %s',filename); end

if ~exist('outtype','var'), outtype = 'double'; end

% allocate space
sz = en-st+1;
if any(sz<1), fprintf('get_block: sz<1\n'); keyboard; end
cumsz = cumsum(sz);

if strcmp(outtype,'double')
  result = zeros(cumsz(end),1);
  fstring = [format '=>double'];
else
  result = zeros(cumsz(end),'int64');
  fstring = '*int64';
end

% read from disk
pos = 1;
for i=1:length(st)
%  if strcmp(type,'byte') & strcmp(outtype,'

  if fseek(f,st(i)*width,'bof')==-1, error('Error seeking to position %d',st(i)); end
  result(pos:cumsz(i)) = fread(f,sz(i),fstring,endian);
  pos = cumsz(i)+1;
end

fclose(f);

% typecast if necessary

if ~strcmp(outtype,'double')
  result = typecast(result,outtype);
end
