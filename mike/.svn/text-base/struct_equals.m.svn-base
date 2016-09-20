function result = struct_equals(A,B)
% struct_equals(A,B)
%
% given two structs, A and B
%
% returns false if:
%    they don't have the same number of fields
%    they have a different set of fieldnames
%          (different orderings OK)
%    for any field
%        the two structs have fields are of different type
%        the two structs have fields of different length
%        the two structs have fields of nonidentical contents
%
% Mike Lawrence 2009-10-16

if ~isstruct(A) | ~isstruct(B), error('A and B must be structs'); end

result = false;

fA = sort(fields(A));
fB = sort(fields(B));

if length(fA)~=length(fB)
  fprintf('structs have different fields\n');
  if ~isempty(setdiff(fA,fB)), fprintf('\n  only in first struct:\n'); disp(setdiff(fA,fB)), end
  if ~isempty(setdiff(fB,fA)), fprintf('\n  only in second struct:\n'); disp(setdiff(fB,fA)), end
  return
end

if slength(A)~=slength(B)
  fprintf('structs have different lengths:  %d vs. %d\n',slength(A),slength(B));
  return
end

result = true;

for f=1:length(fA)

  if ~strcmp(fA{f},fB{f})
    fprintf('nonidentical fieldnames:  %s vs. %s\n',fA{f},fB{f});
    fname = [fA{f} '/' fB{f}];
  else
    fname = fA{f};
  end
    
  a = getfield(A,fA{f});
  b = getfield(B,fB{f});

  if length(a)~=length(b)
    fprintf('field %s has nonidentical lengths:  %d vs. %d\n',fname,length(a),length(b));
    result = false;
  else
    if (isnumeric(a) & isnumeric(b)) | (islogical(a) & islogical(b))
      if any(~nanequals(a,b))
        fprintf('field %s has nonidentical contents\n',fname);
        idx = find(~nanequals(a,b),1);
        fprintf('\te.g. %d vs. %d\n',a(idx),b(idx));
        result = false;
      end
    elseif iscell(a) & iscell(b)
      if ~all(strcmp(a,b))
        fprintf('field %s has nonidentical contents\n',fname);
        idx = find(~strcmp(a,b),1);
        fprintf('\te.g. %s vs. %s\n',a{idx},b{idx});
        result = false;
      end
    else
      fprintf('field %s has nonidentical types\n',fname);
      whos a b
      result = false;
    end
  end
end

