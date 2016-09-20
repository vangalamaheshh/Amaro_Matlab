function C2 = convert_chrXXX(C1,P)
%
% convert_chr(chromosome_list)
%
% converts text chromosome identifiers to numbers 1-24 (X=23, Y=24)
%
% Mike Lawrence 2008-05-01

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

if isnumeric(C1), C2=C1; return; end
if ~iscell(C1), C1={C1}; end      

if ~isempty(P.build)

  if ~isempty(P.refdir)

      % look up number to gene mapping from <build>_info.txt file under refdir
      filename=[P.refdir '/' P.build '_info.txt'];
      if ~exist(filename,'file')
        error(['missing ' filename])  
      end
      CC = load_table(filename);
      %CC = trimStruct(CC,strfindk(CC.use,'D'));
      CC.chromosome=upper(CC.chromosome); 
      M = containers.Map;
      for j=1:length(CC.num);
         M(CC.chromosome{j})=CC.num(j);
         q=split(upper(CC.alias{j}),',');
         for s=1:length(q)
             M(q{s})=CC.num(j);
         end
      end
      C2=NaN*zeros(size(C1));
      UC=upper(unique(sort(C1)));
      m=ismember(UC,M.keys);
      UC=UC(m);
      for u=1:length(UC)
          k=strmatch(UC{u}, upper(C1), 'exact');          
          C2(k)=M(UC{u})+0*k;
      end
      return;
      
  end
  assumed_human = false;
  ct = get_chrcount(P.build);

else
  assumed_human = true;
  ct = 24;
end

if isnumeric(C1), C2=C1; return; end

if ~iscell(C1), C1={C1}; end

if size(C1,1)==1, flag=true; C1=C1'; else flag=false; end

  [ccs cci ccj] = unique(C1);

%  if size(ccs,1)==1    % to fix a weird bug
%    cci=cci(1);        %
%    ccj=ccj(1);        %
%  end                  %

  ccs = regexprep(ccs, '^([Cc][Hh][Rr])', '');
  ccs = regexprep(ccs, '^(Mt|MT)$','0');
  ccs = regexprep(ccs, '^[Mm]$', '0');

  idx = grep('^[XxYy]$',ccs,1);
  if ~isempty(idx)
    if assumed_human
      fprintf('convert_chr: assuming human for chrX/chrY\n');
    end
    ccs = regexprep(ccs, '^[Xx]$', num2str(ct-1));
    ccs = regexprep(ccs, '^[Yy]$', num2str(ct));
  end

  ccs = str2double(ccs);
  C2 = ccs(ccj);

if flag, C2=C2'; end

if any(C2>ct)
  fprintf('WARNING: chromosome number(s) out of range for this genome\n');
  disp(unique(C2(C2>ct)));
end
