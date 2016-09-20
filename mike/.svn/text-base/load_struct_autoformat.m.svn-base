function X = load_struct_autoformat(infile,header_lines,lowercase_fieldnames)
% load_struct_autoformat(infile,header_lines,lowercase_fieldnames)
%
% same as load_struct, but automatically infers format:
% assumes all columns are numeric until proven otherwise.
%
% header_lines is 1 by default.
% if header_lines is set to 0, field names are col1, col2, etc.
%
% if lowercase_fieldnames is 1, fieldnames are converted to lowercase
%
% Mike Lawrence 2009-06-11

if ~exist('header_lines', 'var'), header_lines = 1; end
if ~exist('lowercase_fieldnames', 'var'), lowercase_fieldnames = false; end
if ~exist(infile,'file'), error('%s not found',infile); end

try

fprintf('Loading %s\n',infile);
L = load_lines(infile);
nl = length(L);

fprintf('  Determining column count\n');
for i=1:nl, nt(i,1) = sum(L{i}==char(9)); end
nf = max(nt)+1;

if header_lines>0
  H = split(L{header_lines},char(9));
  L = L(header_lines+1:end);
  nl = length(L);
else
  for i=1:nf, H{i,1} = sprintf('col%d',i); end
end

% remove illegal characters from column headings
% and convert to list of unique field names

H_orig = H;
if lowercase_fieldnames, H = lower(H); end
H = regexprep(H, '\W','');   % remove any characters except A-Z, a-z, 0-9, underscore
H = genvarname(H);

% preserve "end", because it's only going to be a field name, not a variable name

for f=1:nf
  if strcmp(H_orig{f}, 'end')
     H{f} = 'end';
     break
  end
  if strcmp(H_orig{f}, 'End')
     if lowercase_fieldnames, H{f} = 'end';
     else H{f} = 'End'; end
     break
  end
end

% extract data

fprintf('Extracting data\n');
F = repmat({''},nl,nf);
for l=1:nl
  I = split(L{l},char(9));
  F(l,1:length(I)) = I;
end

% build struct (trying to make everything numeric)

fprintf('  Building struct\n');
X = [];
for f=1:nf, fprintf('%d/%d ',f,nf);
  c = F(:,f);
  tmp = str2double(c);
  if ~any(isnan(tmp) & ~strcmpi(c,'nan')), c = tmp; end
  X = setfield(X,H{f},c);
end,fprintf('\n');

fprintf('  Done!\n');

catch me, excuse(me); end
