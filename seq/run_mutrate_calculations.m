function run_mutrate_calculations(C)

aidx = grep('A in',C.cat.name,1);
cidx = grep('C in',C.cat.name,1);
gidx = grep('G in',C.cat.name,1);
tidx = grep('T in',C.cat.name,1);
toA = 1;
toC = 2;
toG = 3;
toT = 4;
toN = 1:4;
plus = grep('^\+$',C.gene.strand,1);
minus = grep('^\-$',C.gene.strand,1);
allgenes = [plus;minus];
igr = C.ng;
total = 1:C.ng;
highexpr = find(C.gene.expr>10);  % (3173 genes)
lowexpr = find(C.gene.expr<5);    % (3727 genes)
intron = grepi('intron',C.cat.name,1);
exon = grepi('exon',C.cat.name,1);
good = grepi('good',C.cat.name,1);


fprintf('\nMUTATIONS AT A:T BASEPAIRS:\n');

fprintf('total A:T\t\t\t\t');
mutrate(C,total,[aidx;tidx],toN)

fprintf('IGR A:T\t\t\t\t\t');
mutrate(C,igr,[aidx;tidx],toN)

fprintf('non-IGR A:T\t\t\t\t');
mutrate(C,allgenes,[aidx;tidx],toN)

fprintf('non-IGR A:T, A on coding strand\t\t');
mutrate(C,plus,aidx,toN,minus,tidx,toN)

fprintf('non-IGR A:T, A on non-coding strand\t');
mutrate(C,minus,aidx,toN,plus,tidx,toN)

fprintf('non-IGR A:T, A on plus strand\t\t');
mutrate(C,allgenes,aidx,toN)

fprintf('non-IGR A:T, A on minus strand\t\t');
mutrate(C,allgenes,tidx,toN)

fprintf('\nMUTATIONS AT C:G BASEPAIRS:\n');

fprintf('total C:G\t\t\t\t');
mutrate(C,total,[cidx;gidx],toN)

fprintf('IGR C:G\t\t\t\t\t');
mutrate(C,igr,[cidx;gidx],toN)

fprintf('non-IGR C:G\t\t\t\t');
mutrate(C,allgenes,[cidx;gidx],toN)

fprintf('non-IGR C:G, C on coding strand\t\t');
mutrate(C,plus,cidx,toN,minus,gidx,toN)

fprintf('non-IGR C:G, C on non-coding strand\t');
mutrate(C,minus,cidx,toN,plus,gidx,toN)

fprintf('non-IGR C:G, C on plus strand\t\t');
mutrate(C,allgenes,cidx,toN)

fprintf('non-IGR C:G, C on minus strand\t\t');
mutrate(C,allgenes,gidx,toN)
