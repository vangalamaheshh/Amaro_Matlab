function X = ibp(individual,chr,pos)
% try BreakPointer on a putative somatic insertion
% individual = e.g. 'MM-0309'

X = [];
X.chr1 = chr;
X.pos1 = pos-20;
X.str1 = 0;
X.chr2 = chr;
X.pos2 = pos+20;
X.str2 = 1;
X.tumreads = 10;
X.normreads = 0;
X.num = 1;

X = tryBP(individual,X,1);
