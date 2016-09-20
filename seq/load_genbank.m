function G = load_genbank(filename)
%
% Mike Lawrence 2008/06/12
%
% currently retrieves: start, end, strand
%                      class (CDS/rRNA/tRNA)
%                      gene, locus_tag, product
%

G.start = zeros(2,1);
G.end = zeros(2,1);
G.strand = cell(2,1);
G.class = cell(2,1);      % CDS, rRNA, tRNA, etc.
G.gene = cell(2,1);
G.locus_tag = cell(2,1);
G.product = cell(2,1);

% open file and seek ahead to FEATURES

in = fopen(filename, 'rt');
while(1)
  L = fgetl(in);
  if L==-1, error('File contains no FEATURES'); end
  if strncmp(L, 'FEATURES', 8), break; end
end

f=1;
while (1)  
  if strncmp(L,'ORIGIN',6), break; end  % end of features
  % process feature
  while(1)
    if length(L)>=21, L1 = L(1:21);
    else, L1 = L(1:end); end
    tmp = regexp(L1, '\s*(\S*)', 'tokens');
    class = tmp{1}{1};
    if strcmp(class, 'CDS'), break; end   % found a feature
    if strcmp(class, 'tRNA'), break; end   % found a feature
    if strcmp(class, 'rRNA'), break; end   % found a feature
    % skip "gene" (redundant)
    L = fgetl(in);   % try next line
    if L==-1, error('Premature file end'); end
  end
  G.class{f} = class;

  % extract coordinates
  L2 = L(22:end);
  L2 = regexprep(L2, 'join\((\d*)\.\.\d*,\d*\.\.(\d*)\)','$1\.\.$2'); % simplify "join"s
  if strncmp(L2, 'complement',10)
    G.strand{f} = '-';
    tmp = regexp(L2,'complement\((\d*)\.\.(\d*)\)','tokens');
    G.start(f) = str2double(tmp{1}{1});
    G.end(f) = str2double(tmp{1}{2});
  else
    G.strand{f} = '+';
    tmp = regexp(L2,'(\d*)\.\.(\d*)','tokens');
    G.start(f) = str2double(tmp{1}{1});
    G.end(f) = str2double(tmp{1}{2});
  end

  G.locus_tag{f} = '';
  G.product{f} = '';
  G.gene{f} = '';
  % read process info until see something at column 6 place
  while(1)
    L = fgetl(in);
    % see if this should be processed elsewhere
    if strncmp(L,'ORIGIN',6), break; end  % end of features
    if length(L)>=21, L1 = L(1:21);
    else L1 = L(1:end); end
    tmp = regexp(L1, '\s*(\S*)', 'tokens');
    class = tmp{1}{1};
    if ~strcmp(class,''), break; end

    % process info
    L2 = L(22:end);
      while(L2(end)~='"')
        L=fgetl(in);
        L2 = [L2 ' ' L(22:end)];
      end
    tmp = regexp(L2,'/([^=]*)="([^"]*)"','tokens');
    info = tmp{1}{1};
    content = tmp{1}{2};
    if strcmp(info,'gene'), G.gene{f} = content;
    elseif strcmp(info,'locus_tag'), G.locus_tag{f} = content;
    elseif strcmp(info,'product'), G.product{f} = content;
    end
  end

  f=f+1;
  % process next feature
end


