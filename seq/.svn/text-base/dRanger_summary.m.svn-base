function dRanger_summary(sample,P)
% dRanger_summary(sample,P)
%
% old parameter style: dRanger_summary(sample,outfile)

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P=impose_default_value(P,'results_name','dRanger_results');
P=impose_default_value(P,'summary_name',[]);
P=impose_default_value(P,'filter_by_tumreads',false);
P=impose_default_value(P,'tumreads_cutoff',[]);
P=impose_default_value(P,'filter_by_S',false);
P=impose_default_value(P,'maximum_S',4);

direc = ['/xchip/tcga_scratch/lawrence/' sample];
name2 = upper(regexprep(sample,'/','-'));
fname = [direc '/' P.results_name '.txt'];
if ~exist(fname,'file'), error('Can''t find file %s',fname);end
Xall = load_struct(fname);
Xall = make_numeric(Xall,{'chr1','chr2','normreads','filter'});

O = [];
  function fprintf_output(message,varargin)
    O = [O sprintf(message,varargin{:})];
  end

try

fprintf_output('BEFORE FILTERING\n');

Xall.inter = (Xall.chr1~=Xall.chr2);
Xall.intra = (Xall.chr1==Xall.chr2);
Xall.som = (Xall.normreads==0);
Xall.germ = (Xall.normreads>0);

fprintf_output('               Total = %d\t(%d inter, %d intra)\n', slength(Xall), sum(Xall.inter), sum(Xall.intra));
fprintf_output('            Germline = %d\t(%d inter, %d intra)\n', sum(Xall.germ), sum(Xall.germ & Xall.inter), ...
        sum(Xall.germ & Xall.intra));
fprintf_output('             Somatic = %d\t(%d inter, %d intra)\n', sum(Xall.som), sum(Xall.som & Xall.inter), ...
        sum(Xall.som & Xall.intra));

fprintf_output('\n');

fprintf_output('AFTER FILTERING\n');

X = reorder_struct(Xall,Xall.filter==0);

if P.filter_by_S
  if isempty(P.maximum_S)
    error('Must specify P.maximum_S');
  else
    fprintf_output('Filtered by S <= %d\n',P.maximum_S);
    X = make_numeric(X,{'S1','S2'});
    X = reorder_struct(X,X.S1<=P.maximum_S & X.S2<=P.maximum_S);
  end
end


if P.filter_by_tumreads
  if isempty(P.tumreads_cutoff)
    error('Must specify P.tumreads_cutoff');
  else
    fprintf_output('Filtered by tumreads >= %d\n',P.tumreads_cutoff);
    X = make_numeric(X,{'tumreads'});
    X = reorder_struct(X,X.tumreads>=P.tumreads_cutoff);
  end
end

fprintf_output('               Total = %d\t(%d inter, %d intra)\n', slength(X), sum(X.inter), sum(X.intra));
fprintf_output('            Germline = %d\t(%d inter, %d intra)\n', sum(X.germ), sum(X.germ & X.inter), sum(X.germ & X.intra));
fprintf_output('             Somatic = %d\t(%d inter, %d intra)\n', sum(X.som), sum(X.som & X.inter), sum(X.som & X.intra));

fprintf_output('\n');

% 4x4 histograms

fprintf_output('Breakpoint locations (after filtering):\n');

X.class1 = regexprep(X.site1,'[: ].*$','');
X.class2 = regexprep(X.site2,'[: ].*$','');
X.class1 = regexprep(X.class1,'\d''-','');
X.class2 = regexprep(X.class2,'\d''-','');
X.class1 = regexprep(X.class1,'Exons','Exon');
X.class2 = regexprep(X.class2,'Exons','Exon');

siteclass = {'Exon','Intron','UTR','IGR'};

X.sc1 = listmap(X.class1,siteclass);
X.sc2 = listmap(X.class2,siteclass);

idx = find(X.sc1 < X.sc2);
tmp = X.sc1(idx); X.sc1(idx)=X.sc2(idx); X.sc2(idx)=tmp;

for i=1:3
  fprintf_output('\n');
  switch i, case 1, Y=X; fprintf_output('TOTAL');
            case 2, Y=reorder_struct(X,X.germ); fprintf_output('GERMLINE');
            case 3, Y=reorder_struct(X,X.som);  fprintf_output('SOMATIC');
  end
  fprintf_output('\n\t');
  h = zeros(length(siteclass),length(siteclass));
  for j=1:slength(Y)
    h(Y.sc1(j),Y.sc2(j)) = h(Y.sc1(j),Y.sc2(j)) + 1;
  end
  for j=1:length(siteclass), fprintf_output('\t%s',siteclass{j}); end
  for k=length(siteclass):-1:1
    fprintf_output('\n\t%s',siteclass{k});
    for j=1:length(siteclass)
      fprintf_output('\t%d',h(k,j));
      if j==k, break; end
    end
  end
  fprintf_output('\n');
end

fprintf('%s',O);

if isempty(P.summary_name), summary_name = [P.results_name '_summary.txt'];
else summary_name = P.summary_name; end
save_textfile(O,[direc '/' summary_name]);

catch me, excuse(me); end

end % main function
