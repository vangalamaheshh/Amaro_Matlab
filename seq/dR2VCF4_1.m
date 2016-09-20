function V=dR2VCF1(X,PAR)
% make/print VCF structure  from dRanger matlab struct
% V=dR2VCF(X,PAR)
% inputs: 
%       X dRanger matlab struct
%     PAR parameter struct
%         PAR.ref_build (default hg18)
%         PAR.ref_area  (default none)
%         PAR.mixture (vector of tumor and normal tissue fractions in tumor, default [0.5 0.5])
%         PAR.vcf_filename  (output vcf filename, default print to screen)
% output:  zipped VCF file with name P.vcf_filename
%          V struct with VCF fields
%  

if (nargin<2)
    PAR.ref_build='HG18';
    PAR.filename='';
end

if (~isfield(PAR,'ref_build'))
    PAR.ref_build='HG18';
end
if (~isfield(PAR,'ref_area'))
    PAR.ref_area='~/Projects/NCBI/build_36.3';
end

if (~isfield(PAR,'mixture'))
    PAR.mixture=[0.5 0.5];
end

if (~isfield(PAR,'min_somaticscore'))
    PAR.min_somaticscore=3;
end

fid=1;
if (isfield(PAR,'vcf_filename'))
    f=PAR.vcf_filename;
    if (length(f)>1)
        fid=fopen(f,'wt');
    end
end

fprintf('fh_dRanger2VCF\n');
fprintf(['  vcf_filename = ' PAR.vcf_filename '\n']);
fprintf('  min_somaticscore = %d\n', PAR.min_somaticscore);
fprintf(['  ref_area = ' PAR.ref_area '\n']);
fprintf(['  ref_build = ' PAR.ref_build '\n']);

V.H.fileformat='VCFv4.1';
V.H.fileDate=datestr(now, 'yyyymmdd');

sam=unique(X.individual);
if (length(sam)>1)
    fprintf(1,' dR2VCF1 is limited to one individual ')
    for i=1:length(sam)
        fprintf(1,'%d\t%s\n ',i,sam{i})   
    end
    return
end

fprintf('  X elements input  = %d\n', length(X.pos1));

X=trimStruct(X,find(X.somatic_score>=PAR.min_somaticscore));

fprintf('  X elements selected  = %d\n', length(X.pos1));


cls=unique(X.class);

tabulate(X.class)

V.H.reference=PAR.ref_build;
V.H.INFO={};
V.H.INFO=[V.H.INFO;{'##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">'}];
V.H.INFO=[V.H.INFO;{'##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">'}];
V.H.INFO=[V.H.INFO;{'##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'}];
V.H.INFO=[V.H.INFO;{'##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of mate breakends">'}];
V.H.INFO=[V.H.INFO;{'##INFO=<ID=SVCLASS,Number=1,Type=String,Description="dRanger class of structural variant">'}];
V.H.INFO=[V.H.INFO;{'##INFO=<ID=SOMATICSCORE,Number=1,Type=Integer,Description="Somatic Score">'}];
V.H.INFO=[V.H.INFO;{'##INFO=<ID=MATECHROM,Number=1,Type=String,Description="Mate Breakend chromosome">'}];
V.H.INFO=[V.H.INFO;{'##INFO=<ID=MATEPOS,Number=1,Type=Integer,Description="Mate Breakend coordinate">'}];
V.H.INFO=[V.H.INFO;{'##INFO=<ID=MATESTRAND,Number=1,Type=String,Description="Mate Breakend strand (+/-)">'}];
V.H.INFO=[V.H.INFO;{'##INFO=<ID=STRAND,Number=1,Type=String,Description="Breakend strand (+/-)">'}];

V.H.FMT={};
V.H.FMT=[V.H.FMT;{'##FORMAT=<ID=NALT,Number=1,Type=Integer,Description="number of ALT support Reads">'}];
V.H.FMT=[V.H.FMT;{'##FORMAT=<ID=NREF,Number=1,Type=Integer,Description="number of REF support Reads">'}];

V.H.FILT={};
V.H.FILT=[V.H.FILT;{sprintf('##FILTER=<ID=LOWSCORE,Description="Somatic Score < %.1f">',PAR.min_somaticscore)}];

V.H.SAM={};
V.H.sam={};
for s=1:length(sam)
    mx=PAR.mixture;
    V.H.sam=[V.H.sam; {[sam{s} ':N']}];
    V.H.sam=[V.H.sam; {[sam{s} ':T']}];
    V.H.SAM=[V.H.SAM;{['##SAMPLE=<SampleName=<' sam{s} ':N>,Genomes=<Germline>,Description=<"Patient ' sam{s} ' germline genome">>']}];
    V.H.SAM=[V.H.SAM;{['##SAMPLE=<SampleName=<' sam{s} ':T>,Genomes=<Germline,Tumor>,Mixture=<' num2str(mx,'%.1f,%.1f') '>,Description=<"Patient germline genome","Patient tumor genome">>']}];
end

% record column header line
head={'#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'};
V.H.H=head{1};
for s=2:length(head)
    V.H.H=[V.H.H sprintf('\t%s',head{s})];
end   

for s=1:length(V.H.sam)
    V.H.H=[ V.H.H sprintf('\t%s',V.H.sam{s})];
end       



N=length(X.pos1);

% required fields
V.CHROM=num2chrom(X.chr1);
V.POS=X.pos1;
V.ID=cellstr([num2str(X.num) repmat(':',N,1) char(X.individual) repmat(':1',N,1)]);
V.ID=regexprep(V.ID,' ','');
V.REF=repmat({'.'},N,1);    
V.ALT=repmat({'.'},N,1);    
V.QUAL=round(-10*log10(1-X.quality));
V.QUAL(V.QUAL>99)=99;
V.FILTER=repmat({'PASS'},N,1);    
k=find(X.somatic_score<PAR.min_somaticscore);
V.FILTER(k)={'LOWSCORE'};    

% INFO fields
V.SVTYPE=repmat({'BND'},N,1);
V.STRAND=repmat({'-'},N,1);
V.STRAND(find(~X.str1))={'+'};
V.MATECHROM=num2chrom(X.chr2);
V.MATEPOS=X.pos2;
V.MATESTRAND=repmat({'-'},N,1);
V.MATESTRAND(find(~X.str2))={'+'};
V.MATEID=cellstr([num2str(X.num) repmat(':',N,1) char(X.individual) repmat(':2',N,1)]);
V.MATEID=regexprep(V.MATEID,' ','');
V.SVCLASS=X.class;
V.SOMATICSCORE=X.somatic_score;
% set postion CI to +/- 100 bp (~Mike L)
V.CIPOS=repmat([-100 100],N,1);
V.CIPOS(find(X.BPresult),:)=0*V.CIPOS(find(X.BPresult),:);
V.IMPRECISE=true(N,1);
V.IMPRECISE(find(X.BPresult))=false;

% FMT fields
V.GT=repmat({'.'},N,1);
V.NALT(:,1)=X.normreads;
V.NALT(:,2)=X.tumreads;

% flip for MATE breakend records
VM=V;
VM.CHROM=num2chrom(X.chr2);
VM.POS=X.pos2;
VM.ID=V.MATEID;
VM.REF=repmat({'.'},N,1);    
VM.ALT=repmat({'.'},N,1);    
VM.QUAL=V.QUAL;
VM.FILTER=V.FILTER;

% INFO fields
VM.SVTYPE=repmat({'BND'},N,1);
VM.STRAND=V.MATESTRAND;
VM.MATECHROM=V.CHROM;
VM.MATEPOS=X.pos1;
VM.MATESTRAND=V.STRAND;
VM.MATEID=V.ID;
VM.SVCLASS=X.class;
VM.SOMATICSCORE=V.SOMATICSCORE;
VM.CIPOS=V.CIPOS;
VM.IMPRECISE=V.IMPRECISE;

% FMT fields
VM.GT=V.GT;
VM.NALT=V.NALT;
         
% catenate V V2
V=mergeStruct(V,VM);

if (V.N>0)
    % sort on chromosome pos
    x=chrom2num(V.CHROM)*1e10+V.POS;
    [x,k]=sort(x);
    V=trimStruct(V,k);
end

% init reference genome lookup...
DOREF=length(PAR.ref_area)>1;
a0=0;  

% loop over events
for i=1:V.N
    
    % chromosome number
    a=chrom2num(V.CHROM(i));

    % get reference base
    if (DOREF)
        if (a~=a0)
            a0=a;
            %f=[PAR.ref_area '/hs_ref_chr' char(num2chrom(a0)) '.mat']
            %HG=load(f);
            f=[PAR.ref_area '/chr' char(num2chrom(a0)) '.txt']
            HG.s=loadRefText(f);
        end
        V.REF{i,1}=HG.s(V.POS(i));
    else
        V.REF{i,1}='.';
    end
    
    ptag=sprintf('chr%s:%d',V.MATECHROM{i},V.MATEPOS(i));
    ttag=V.REF{i};
    
    % ILLUMINA RP +- PROPER ORIENTATION
    
    if ( (V.STRAND{i}=='+')&    (V.MATESTRAND{i}=='+') )     
           V.ALT{i,1}=sprintf('%s]%s]',ttag,ptag);  
    end
    
    if ( (V.STRAND{i}=='+')&    (V.MATESTRAND{i}=='-') )
           V.ALT{i,1}=sprintf('%s[%s[',ttag,ptag);
    end
    
    if ( (V.STRAND{i}=='-')&    (V.MATESTRAND{i}=='+') )
           V.ALT{i,1}=sprintf('[%s[%s',ptag,ttag);
    end
    
    if ( (V.STRAND{i}=='-')&    (V.MATESTRAND{i}=='-') )
           V.ALT{i,1}=sprintf(']%s]%s',ptag,ttag);
    end
                   
end

% write vcf
fprintf(fid,'##fileformat=%s\n',V.H.fileformat);
fprintf(fid,'##fileDate=%s\n',V.H.fileDate);
fprintf(fid,'##reference=%s\n',V.H.reference);
for i=1:length(V.H.INFO)    
    fprintf(fid,'%s\n',V.H.INFO{i});
end
for i=1:length(V.H.FMT)    
    fprintf(fid,'%s\n',V.H.FMT{i});
end
for i=1:length(V.H.FILT)    
    fprintf(fid,'%s\n',V.H.FILT{i});
end
for i=1:length(V.H.SAM)    
    fprintf(fid,'%s\n',V.H.SAM{i});
end

fprintf(fid,'%s\n',V.H.H);

for i=1:V.N
    
    fprintf(fid,'%s\t%d\t%s\t%s\t%s\t%d\t%s',V.CHROM{i},V.POS(i),V.ID{i},V.REF{i},V.ALT{i},V.QUAL(i),V.FILTER{i});
    fprintf(fid,'\tMATEID=%s',V.MATEID{i});

    if (V.IMPRECISE(i))
        fprintf(fid,';IMPRECISE;CIPOS=%d,%d',V.CIPOS(i,1),V.CIPOS(i,2));
    end

    fprintf(fid,';SVTYPE=%s;SVCLASS=%s',V.SVTYPE{i},V.SVCLASS{i});
    fprintf(fid,';STRAND=%s;MATESTRAND=%s;MATEPOS=%d;MATEID=%s',V.STRAND{i},V.MATESTRAND{i},V.MATEPOS(i),V.MATEID{i});
    if (~strcmp(V.MATECHROM{i},V.CHROM{i}));
        fprintf(fid,';MATECHROM=%s',V.MATECHROM{i});
    end
    fprintf(fid,';SOMATICSCORE=%d',V.SOMATICSCORE(i));
 
    fprintf(fid,'\tNALT');
    
    for s=1:length(V.H.sam);        
        fprintf(fid,'\t%d',V.NALT(i,s));
    end
        
    fprintf(fid,'\n');
    

end




function s=loadRefText(f)

fid = fopen(f,'rt');
s = textscan(fid, '%s','bufsize',2^30);
s=char(s{1});
fclose(fid);

 
function test0

clear
cd /Users/stewart/Projects/dRanger/tools/matlab
%rmpath('/Users/stewart/Projects/Spanner/tools/matlab')
%rmpath('/Users/stewart/Projects/matlab')
%rmpath('/Users/stewart/Documents/MATLAB')
A='~/Projects/dRanger/data/Prostate/'
%A='~/Projects/dRanger/data/GBM/'

d=dir([A '*.mat']);
X=[]; i=1;
for i=1:length(d)
    f=d(i).name
    load([A f]);        

    P.ref_build='hg18';
    P.vcf_filename=[A regexprep(f,'.mat$','.vcf')];
    P.ref_area='/Volumes/cga1/annotation/db/ucsc/hg18_v2';
    %P.ref_area='~/Projects/NCBI/build_36.3';
    P.min_somaticscore=4;    
    
    V=dR2VCF4_1(X,P);
    
    cmd=sprintf('/usr/local/bin/bgzip  %s  ',P.vcf_filename)
    system(cmd)
end

function test1
    A='~/Projects/Cancer/dRanger/An_TCGA_Benchmark4/'
    system(['mkdir -p ' A])
    f0='/local/cga-fh/cga/An_TCGA_Benchmark4/Individual/HCC1143.n20t80/jobs/wgs/ra/CC1143.n20t80.dRanger_results.detail.all.mat'    
    system(['rsync -av ' f0 ' ' A])
    [a f]=fileparts(f0); f=[f '.mat']
    load([A '/' f]);
    P.ref_build='hg19';
    P.vcf_filename=[A regexprep(f,'.mat$','.vcf')];
    P.ref_area='/xchip/cga1/annotation/db/ucsc/hg19';
    P.min_somaticscore=4;        
    V=dR2VCF4_1(X,P);
    
    cmd=sprintf('bgzip  %s  ',P.vcf_filename)
    system(cmd)
end



