function [G,M,u_nsc]=read_pharma_data(fname,cells,just_read)

dlm=',';

% 89 ,M ,-4.00,Non-Small Cell Lung ,NCI-H23 ,1 ,1 ,4.209,3 ,3 ,0.361
%  1  2    3       4                 5       6  7   8    9  10  11
%[nsc,unit,lconc,panel,cell_line,panelnbr,cell_line_nbr,nloggi50,indn,totn,stddev]=...
%    textread(fname,'%d,%s,%f,%s,%s,%d,%d,%f,%d,%d,%f');

if exist([ fname '_small.mat'],'file')
  load([ fname '_small.mat']);
else
  fid=fopen([ fname '.txt' ]);
  C=textscan(fid,'%d%s%f%s%s%d%d%f%d%d%f','delimiter', ',','headerLines',1);
  fclose(fid);
  C{4}=[];
  C{5}=[];
  save([ fname '_small.mat'],'C');
end

if exist('just_read','var') && just_read
  return
end

if (0)
  clt=double(C{6})+double(C{7})/100;
  [u,i,j]=unique(clt);

  cell_lines=find(q>25000);
% checked that these are indeed the NCI60 lines;
% cell_lines=[     1     2     3     5     6     7     8     9 ...
%                 11    23    24    25    26    27 ...
%                 29    30    50    51    54    55    59    60  ...
%                 61    62    70    71    72    73    74    75  ...
%                 79    80    81    82    83    84    94    95  ...
%                 97    98    99   101   104   105   109   110  ...
%                112   113   114   115   117   118   121   122  ...
%                126   127   129   130   131   132 ];


   maxnsc=max(C{1});
   S=sparse(maxnsc,length(cell_lines));
end


clt=double(C{6})+double(C{7})/100;
cells2=cells(:,1)+cells(:,2)/100;

idx=find(ismember(clt,cells2));
clt2=clt;
for i=1:length(cells2)
  clt2(clt==cells2(i))=i;
end

maxnsc=max(C{1});
maxtotn=max(C{9});

nsc=sparse(double(maxnsc),1);
u_nsc=unique(double(C{1}));
nsc(u_nsc)=1:length(u_nsc);

M=cell(length(u_nsc),length(cells));

for i=1:length(idx)
  ci=idx(i);
  clid=clt2(ci);
  cmp=full(nsc(C{1}(ci)));
%  if C{9}(ci)==C{10}(ci) % it is the max
%    G(cmp,clid)=C{8}(ci);
%  end   
  q=M{cmp,clid};
  if ~isempty(q)
    q{end+1}=[ ci double(C{9}(ci)) double(C{10}(ci)) double(C{8}(ci)) double(C{11}(ci)) ];
%    disp([ cmp clid]);
  else
    q={[ ci double(C{9}(ci)) double(C{10}(ci)) double(C{8}(ci)) double(C{11}(ci)) ]};
  end
  M{cmp,clid}=q;
  if mod(i,10000)==0
    disp(i);
  end
end

G=NaN*ones(size(M,1),size(M,2),4);
for i=1:size(M,1)
    for j=1:size(M,2)
        if ~isempty(M{i,j})
            x=cat(1,M{i,j}{:});
            [xm,xi]=sortrows(x(:,[2 4])); % largest mean from the ones with the same number of experiments
            xi=xi(end);
            G(i,j,1)=x(xi,4); %mean
            G(i,j,2)=x(xi,5); %std
            G(i,j,3)=x(xi,2); %nexp
            G(i,j,4)=xi; %row in M
        end 
    end 
end 

% add nsc data
%G.dat=G;
%G.gacc=
