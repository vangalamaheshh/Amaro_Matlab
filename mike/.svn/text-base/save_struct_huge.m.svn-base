function save_struct_huge(S,filename,noheader)
% save_struct_huge(S, filename)
%
% slow version of save_struct that writes directly to disk instead of using a memory buffer.
%
% writes a tab-delimited table from the given struct.
% writes a header line based on the field names.
%
% fields are written in the order in which they appear in the struct.
%
% Mike Lawrence 2008-08-04

omit_headers = false;
if exist('noheader','var')
  if strcmpi(noheader,'no_headers'), omit_headers = true;
  else error('third parameter should be "no_headers" or nothing.');
  end
end

% pre-process fields

fld = fieldnames(S);
nf = length(fld);
ftype = nan(nf,1);
nr = -1;
C = cell(nf,1);
for f=1:nf
  C{f} = getfield(S,fld{f});    % get column
  % check for legal type
  if isempty(C{f})
    ftype(f) = 0; % blank
  elseif isnumeric(C{f}) || islogical(C{f})
    ftype(f) = 1; % number
  else
    elem = C{f}{1};
    if sum(size(elem)>1)>1, error('Field %s is multidimensional', fld{f}); end
    if iscell(elem), error('Field %s is a cell', fld{f}); end
    ftype(f) = 2; % string
  end
  % check for consistent column length
  if nr==-1, nr=length(C{f}); end
  if nr~=length(C{f}), error('Field %s is a different length', fld{f}); end
end    

% write file

out = fopen(filename,'wt');

% header
if ~omit_headers
  for f=1:nf
    fprintf(out,'%s',fld{f});
    if f<nf, fprintf(out,'\t'); else fprintf(out,'\n'); end
  end
end

% body
for r=1:nr
  if ~mod(r,10000), fprintf('%d/%d ',r,nr); end
  for f=1:nf
    switch ftype(f)
      case 0, % empty
      case 1, fprintf(out,'%d',C{f}(r));    % number
      case 2, fprintf(out,'%s',C{f}{r});    % string
      otherwise, error('unknown column type');
    end
    if f<nf, fprintf(out,'\t'); else fprintf(out,'\n'); end
  end
end

fclose(out);
