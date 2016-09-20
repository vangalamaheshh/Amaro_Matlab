% from Mike Chapman 5/31/2010 12:33am
%
% Comparison of intronic mutation rates in expressed
% vs. non-expressed genes (Mike L. or Gaddy - I've attached the two gene
% lists to this email; I did this based on Affy A/P calls at 2
% thresholds: >5% (expressed) vs. <5% (not expressed) P calls and >0%
% vs. 0% P calls. Either list could be used). 

build = 'hg18'

dr = '/xchip/cga1/lawrence/db/mmexpr';
e5 = load_struct([dr '/expr5pct.txt']);
ne5 = load_struct([dr '/nonexpr5pct.txt']);
e1 = load_struct([dr '/expr1.txt']);
ne1 = load_struct([dr '/nonexpr1.txt']);
e50 = load_struct([dr '/expr50pct.txt']);
ne50 = load_struct([dr '/nonexpr50pct.txt']);
e75 = load_struct([dr '/expr75pct.txt']);
ne25 = load_struct([dr '/nonexpr25pct.txt']);

all=[];
all.name = union([e5.gene;ne5.gene;e1.gene;ne1.gene;e50.gene;ne50.gene;e75.gene],ne25.gene);

e5 = listmap(e5.gene,all.name);
ne5 = listmap(ne5.gene,all.name);
e1 = listmap(e1.gene,all.name);
ne1 = listmap(ne1.gene,all.name);
e50 = listmap(e50.gene,all.name);
ne50 = listmap(ne50.gene,all.name);
e75 = listmap(e75.gene,all.name);
ne25 = listmap(ne25.gene,all.name);

all.expr = 5*ones(slength(all),1);
all.expr(e75)=6;
all.expr(ne50)=4;
all.expr(ne25)=3;
all.expr(ne5)=2;
all.expr(ne1)=1;

all.name = regexprep(all.name,'^(\d*)-Mar','MARCH$1');
all.name = regexprep(all.name,'^(\d*)-Sep','SEPT$1');
all.name = regexprep(all.name,'^(\d*)-Dec','DEC$1');
all = sort_struct(all,'name');

save_struct(all,[dr '/mm_expr_genes.txt']);

% match to Refseq

R = load_refseq(build)

aidx = match_genelists(upper(R.gene),upper(all.name));

