function C2 = convert_chr_backXXX(C1,P)
%
% convert_chr_back(chromosome_list)
%
% converts numbers 1-24 to chromosome identifiers (1-->chr1, 23-->chrX, 24-->chrY)
%
% Mike Lawrence 2009-05-19

if ~exist('P','var'), P=[]; end

if ischar(P)
  build=P;
  P=[];
  P.build=build;
end

if ~exist('build','var')
  build = [];
end

P = impose_default_value(P,'build',build);
P = impose_default_value(P,'refdir',[]);

if ~isempty(P.build)
    
  if ~isempty(P.refdir)
      filename=[P.refdir '/' P.build '_info.txt'];
      if ~exist(filename,'file')
        error(['missing ' filename])  
      end
      CC = load_table(filename);
      %CC = trimStruct(CC,strfindk(CC.use,'D'));
      M = containers.Map;
      for j=1:length(CC.num);
         M(CC.num(j))=CC.chromosome{j});
      end
      C2=repmat({''},(size(C1));
      UC=unique(sort(C1));
      m=ismember(UC,M.keys);
      UC=UC(m);
      for u=1:length(UC)
          k=find(UC(u),C1);
          C2(k)=M(UC{u})+0*k;
      end
      return;
  end
  
    
  ct = get_chrcount(P.build);
else
  % assume human
  ct = 24;
end


if isempty(C1)
  C2 = C1;
else
  if isnumeric(C1)
    C1num=C1;
    for i=1:length(C1), C2{i,1} = num2str(C1(i)); end
  else
    C1num=str2double(C1);
    C2 = C1;
  end

  if any(C1num>ct)
    fprintf('WARNING: chromosome number(s) out of range for this genome\n');
    disp(unique(C1num(C1num>ct)));
  end

  C2 = regexprep(C2,'(.*)','chr$1');
  
  C2 = regexprep(C2,'chr0','chrM');
  C2 = regexprep(C2,['chr' num2str(ct-1)],'chrX');
  C2 = regexprep(C2,['chr' num2str(ct)],'chrY');

  C2 = regexprep(C2,'^chrchr','chr');
end
