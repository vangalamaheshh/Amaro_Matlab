function [ratios, cis] = analyze_strand_bias(C,file_pointer)

if ~exist('file_pointer','var'), file_pointer = 1; end  % default = print to screen

require_fields(C,{'cat','n','N','gene'});
require_fields(C.n,{'gc'});   % gc = gene x category
require_fields(C.N,{'gc'});
require_fields(C.gene,{'strand'});
require_fields(C.cat,{'base','left','right'});

plus = find(strcmp(C.gene.strand,'+'));
minus = find(strcmp(C.gene.strand,'-'));

ratios = nan(9, 1); 
cis = nan(9, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mutation rate in CpG dinucleotides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cpgidx = find((strcmp(C.cat.base,'C') & strcmp(C.cat.right,'G')) |...
              (strcmp(C.cat.base,'G') & strcmp(C.cat.left,'C')));
n = sum(sum(C.n.gc([plus;minus],cpgidx)));
N = sum(sum(C.N.gc([plus;minus],cpgidx)));
[ratio sd] = ratio_and_sd(n,N); ci = 1.84*sd;
fprintf(file_pointer,'%% CpG mutations:   total                  n = %5d  N = %12d  rate = %0.2d +- %0.2d\n',n,N,ratio,sd);
ratios(1) = ratio; 
cis(1, :) = sd;

% strand bias
c_cpg = find((strcmp(C.cat.base,'C') & strcmp(C.cat.right,'G')));
g_cpg = find((strcmp(C.cat.base,'G') & strcmp(C.cat.left,'C')));

% CpG mutations at C on coding strand
n = sum(sum(C.n.gc(plus,c_cpg))) + sum(sum(C.n.gc(minus,g_cpg)));
N = sum(sum(C.N.gc(plus,c_cpg))) + sum(sum(C.N.gc(minus,g_cpg)));
[ratio sd] = ratio_and_sd(n,N); ci = 1.84*sd;
fprintf(file_pointer,'%% CpG mutations:   C on coding strand     n = %5d  N = %12d  rate = %0.2d +- %0.2d\n',n,N,ratio,sd);
ratios(2) = ratio;
cis(2, :) = sd;

n1=n;N1=N;

% CpG mutations at C on noncoding strand
n = sum(sum(C.n.gc(minus,c_cpg))) + sum(sum(C.n.gc(plus,g_cpg)));
N = sum(sum(C.N.gc(minus,c_cpg))) + sum(sum(C.N.gc(plus,g_cpg)));
[ratio sd] = ratio_and_sd(n,N); ci = 1.84*sd;
fprintf(file_pointer,'%% CpG mutations:   C on noncoding strand  n = %5d  N = %12d  rate = %0.2d +- %0.2d\n',n,N,ratio,sd);
ratios(3) = ratio;
cis(3, :) = sd;

n2=n;N2=N;

dir = (n1/N1)>(n2/N2);
p = compare_ratios(n1,N1,n2,N2);
fprintf(file_pointer,'%%        p = %0.2d    dir = %d\n', p,dir);

fprintf(file_pointer,'\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mutation rate at A:T basepairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

atidx = find((strcmp(C.cat.base,'A') | strcmp(C.cat.base,'T')));
n = sum(sum(C.n.gc([plus;minus],atidx)));
N = sum(sum(C.N.gc([plus;minus],atidx)));
[ratio sd] = ratio_and_sd(n,N); ci = 1.84*sd;
fprintf(file_pointer,'%% A:T mutations:   total                  n = %5d  N = %12d  rate = %0.2d +- %0.2d\n',n,N,ratio,sd);
ratios(4) = ratio;
cis(4, :) = sd;

% strand bias
aidx = find(strcmp(C.cat.base,'A'));
tidx = find(strcmp(C.cat.base,'T'));

% A:T mutations at A on coding strand
n = sum(sum(C.n.gc(plus,aidx))) + sum(sum(C.n.gc(minus,tidx)));
N = sum(sum(C.N.gc(plus,aidx))) + sum(sum(C.N.gc(minus,tidx)));
[ratio sd] = ratio_and_sd(n,N); ci = 1.84*sd;
ratios(5) = ratio;
cis(5, :) = sd;

fprintf(file_pointer,'%% A:T mutations:   A on coding strand     n = %5d  N = %12d  rate = %0.2d +- %0.2d\n',n,N,ratio,sd);

n1=n;N1=N;

% A:T mutations at A on noncoding strand
n = sum(sum(C.n.gc(minus,aidx))) + sum(sum(C.n.gc(plus,tidx)));
N = sum(sum(C.N.gc(minus,aidx))) + sum(sum(C.N.gc(plus,tidx)));
[ratio sd] = ratio_and_sd(n,N); ci = 1.84*sd;
ratios(6) = ratio;
cis(6, :) = sd;
fprintf(file_pointer,'%% A:T mutations:   A on noncoding strand  n = %5d  N = %12d  rate = %0.2d +- %0.2d\n',n,N,ratio,sd);

n2=n;N2=N;

dir = (n1/N1)>(n2/N2);
p = compare_ratios(n1,N1,n2,N2);
fprintf(file_pointer,'%%        p = %0.2d    dir = %d\n', p,dir);

fprintf(file_pointer,'\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mutation rate at C:G basepairs not in CpG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cpgidx = find((strcmp(C.cat.base,'C') & strcmp(C.cat.right,'G')) |...
              (strcmp(C.cat.base,'G') & strcmp(C.cat.left,'C')));
cgidx = find((strcmp(C.cat.base,'C') | strcmp(C.cat.base,'G')));
cgidx = setdiff(cgidx,cpgidx);
n = sum(sum(C.n.gc([plus;minus],cgidx)));
N = sum(sum(C.N.gc([plus;minus],cgidx)));
[ratio sd] = ratio_and_sd(n,N); ci = 1.84*sd;
ratios(7) = ratio;
cis(7, :) = sd;
fprintf(file_pointer,'%% C:G nonCpG muts: total                  n = %5d  N = %12d  rate = %0.2d +- %0.2d\n',n,N,ratio,sd);

% strand bias
cidx = setdiff(find(strcmp(C.cat.base,'C')),cpgidx);
gidx = setdiff(find(strcmp(C.cat.base,'G')),cpgidx);

% C:G mutations at C on coding strand
n = sum(sum(C.n.gc(plus,cidx))) + sum(sum(C.n.gc(minus,gidx)));
N = sum(sum(C.N.gc(plus,cidx))) + sum(sum(C.N.gc(minus,gidx)));
[ratio sd] = ratio_and_sd(n,N); ci = 1.84*sd;
ratios(8) = ratio;
cis(8, :) = sd;
fprintf(file_pointer,'%% C:G nonCpG muts: C on coding strand     n = %5d  N = %12d  rate = %0.2d +- %0.2d\n',n,N,ratio,sd);

n1=n;N1=N;

% C:G mutations at C on noncoding strand
n = sum(sum(C.n.gc(minus,cidx))) + sum(sum(C.n.gc(plus,gidx)));
N = sum(sum(C.N.gc(minus,cidx))) + sum(sum(C.N.gc(plus,gidx)));
[ratio sd] = ratio_and_sd(n,N); ci = 1.84*sd;
ratios(9) = ratio;
cis(9, :) = sd;
fprintf(file_pointer,'%% C:G nonCpG muts: C on noncoding strand  n = %5d  N = %12d  rate = %0.2d +- %0.2d\n',n,N,ratio,sd);

n2=n;N2=N;

dir = (n1/N1)>(n2/N2);
p = compare_ratios(n1,N1,n2,N2);
fprintf(file_pointer,'%%        p = %0.2d    dir = %d\n', p,dir);


