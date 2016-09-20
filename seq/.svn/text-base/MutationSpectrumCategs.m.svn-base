function [XY,X,C1]=MutationSpectrumCategs(fmaf,covfile,catfile, custom_coverage_WGS)
%function [XY,X,C1]=MutationSpectrumCategs(fmaf,covfile,catfile)
%
% builds XY struct for mutation contecxt lego plots
% inputs:   fmaf = name of MutSig MAF file with context65 field
%               can also be a MutSig Structure with the same fields as a MAF
%           covfile = name of MutSig coverage file
%               can also be a C1 coverage structure
%               can also be the string "exome" or "genome", in which case a default coverage distribution will be used,
%                       based on a precomputed average from many previous well-covered samples.
%                       The number of patients will be estimated from length(unique(maf.patient))--
%                           this only matters for the absolute rate estimate, if it is designated to be plotted on the z-axis.
%           category file (optional - defaults to /xchip/cga1/lawrence/db/hg19/c65e/categs.txt
%
% outputs:  XY: structure for lego plot
%               XY.n:  8x12 matrix of mutation counts per lego bin
%               XY.ncov: 8x12 matrix of base coverage per lego bin
%               XY.cat:  8x12 cell matrix of context+mutation categories per lego bin
%                   eg. G-C>A-T is a C>A mutation after a G before a T
%               XY.col:  8x12x3 matrix of RGB colors per lego bin
%           X: MAF struction with Mutsig fields & added CT field with mutation category cellstr
%           C1:coverge structure stripped down to a few fields including orig_cov
%
%
% CS:  04Mar2012, based on Mike Lawrence's draw_mutation_spectrum_3d_barplot
% ML:  15Mar2012, added option to use precomputed average coverage ("exome" or "genome")
% CS:  20Jul2012, also allow custom coverage vector (65x1) in covfile
%

if (nargin<1)
    fmaf='/xchip/cga1/firehose_output/Sigma_Pediatric/Individual_Set/AN_SIGMA_Rhabdoid_16Feb2012_Diagnostic/mutsig_1.5_WES/*coverage.mat/AN_SIGMA_Rhabdoid_16Feb2012_TX/mutsig2.0/AN_SIGMA_Rhabdoid_16Feb2012_Diagnostic.final_analysis_set.maf'
end

%if (nargin==4)
%  custom_coverage_WGS = varargin{1};
%end 

if ischar(fmaf)
    X=load_table(fmaf);
    k=strfindk(X.Hugo_Symbol,'Hugo_Symbol');
    if (length(k)>0)
        fprintf(' multiple headers in maf... trying to fix...');
        k1=1:length(X.Hugo_Symbol); k1(k)=[];
        X=trimStruct(X,k1);
        X.context65=strnum(X.context65);
        X.t_alt_count=str2num(X.t_alt_count);
        X.t_ref_count=strnum(X.t_ref_count);
        if isfield(X,'i_tumor_f')
            X.i_tumor_f=strnum(X.i_tumor_f);
        end
        if ~isfield(X,'i_t_lod_fstar')
            X.i_t_lod_fstar=strnum(X.i_t_lod_fstar);
        end
    end
    
else
    X=fmaf;
end

if ~isfield(X,'Variant_Type') && isfield(X,'classification')
    X = rename_field(X,'classification','Variant_Type');
end
if ~isfield(X,'Variant_Type')
    error(' need "Variant_Type" field in maf')
else
    X=trimStruct(X,strfindk(X.Variant_Type,'SNP'));
    if (X.N<1)
        error(['No SNPs in ' fmaf])
    end

end
if ~isfield(X,'patient')
    X.patient=0;
end

if isnumeric(X.patient)
    % older mutsigs have a patient number rather than name.
    % try to replace patient number with the name from Tumor_Sample_Barcode
    % or convert the number to cellstr
    if isfield(X,'Tumor_Sample_Barcode')
        X.patient=regexprep(X.Tumor_Sample_Barcode,'-Tumor$','');
        X.patient=regexprep(X.patient,'-1$','');
    else
        X.patient=cellstr(num2str(X.patient));
    end
end

npat = length(unique(X.patient));

if (nargin<2)
    %covfile='/xchip/cga1/firehose_output/Sigma_Pediatric/Individual_Set/AN_SIGMA_Rhabdoid_16Feb2012_Diagnostic/mutsig_1.5_WES/*coverage.mat/AN_SIGMA_Rhabdoid_16Feb2012_TX/mutsig2.0/AN_SIGMA_Rhabdoid_16Feb2012_Diagnostic.coverage.mat'
    covfile=[];
end
if isempty(covfile)
    %covfile='~/Projects/Cancer/tools/matlab/exome.coverage.mat';
    covfile='unit';
end

if ischar(covfile)
    if strcmpi(covfile,'exome')
        fprintf('Using precomputed average exome coverage\n');
        N = average_exome_coverage();
        C1=[]; C1.orig_cov = repmat(N',npat,1);
    elseif strcmpi(covfile,'genome')
        fprintf('Using precomputed average genome coverage\n');
        N = average_genome_coverage();
        C1=[]; C1.orig_cov = repmat(N',npat,1);
    elseif strcmpi(covfile,'unit')
        fprintf('Not using coverage for mutation rates\n');
        N = 0.5*ones(64,1)/npat;
        C1=[]; C1.orig_cov = repmat(N(:)',npat,1);
        
    else
        load(covfile);
        ff=fieldnames(C1); ff=ff(strfindk(ff,'gene'));
        C1=rmfield(C1,ff);
        ff=fieldnames(C1); ff=ff(strfindk(ff,'fcov'));
        C1=rmfield(C1,ff);
        ff=fieldnames(C1); ff=ff(strfindk(ff,'targ'));
        C1=rmfield(C1,ff);
        catfile='~/Projects/Cancer/tools/matlab/categs.txt';
        if (~exist(catfile))            
          %catfile='/xchip/cga1/lawrence/db/hg19/c65e/categs.txt';
          catfile='/xchip/cga/reference/hg19/annotation/db/tracks/hg19/c65e/categs.txt';
        end   
    end
elseif isnumeric(covfile)
    C1.orig_cov = repmat(covfile(:)',npat,1);
    if ((length(covfile)>65)|(length(covfile)<64))
        error(sprintf('coverage vector should be length 64, find %d',length(covfile)));
    end
else
    C1=covfile;
end


if ((nargin<3) & ~exist('catfile','var'))
   
    catfile=false;
    % conxtext struct
    base='ACGT';
    n=0; C=[];
    for i=1:4; for j=1:4;  for k=1:4
        n=n+1;
        C.name(n,1)={[base(i) ' in ' base(j) '_' base(k)]};
        C.id(n,1)=n;
    end;    end ;end
    n=n+1;
    C.name(n,1)={'any N'};
    C.id(n,1)=n;
    C.N=length(C.id);
    

    
end
if isstruct(catfile)
    C=catfile;
elseif ischar(catfile)
    categs=load_table(catfile);
    categs.longname=categs.name;
    q=regexp(categs.name,':','split');
    categs.name=cellfun(@(x) x(1),q);
    % reduce 845 3-base categories to 65 3-base context categories
    k = map_categories_to_65(catfile);
    [C.name i]=unique(categs.name);
    C.id=k(i);
    C=trimStruct(C,sort(C.id));
    C.N=length(C.id);
end

if ~exist('C','var'), error('need mutation categ info'); end


% context65
if ~isfield(X,'context65')
    if ~isfield(X,'ref_context')
        error(' need "context65" or "ref_context" field in maf')
    end
    % calculate context65 from ref_context if needed
    [X.context65]=ref_context_to_context65(X.ref_context,C)
end

% coverage
if (size(C1.orig_cov,2)>65)
    C.Ncov = round_collapse_categories_to_65(sum(C1.orig_cov,1)');
else  % already in 65 categories
    C.Ncov=sum(C1.orig_cov,1);
end
% SNP from->to convention (flip  such that from is always C or T)
FT=cellfun(@(x) x(1), C.name);
kFLIP1=find( (FT(:,1)=='G')|(FT(:,1)=='T'));
FT(kFLIP1,:)=tr(FT(kFLIP1,:),'ACGTacgt','TGCAtgca');
CT=C.name;
CT0=char(CT); CT0=CT0(:,[6 8]);
CT0F=char(CT); CT0F=CT0F(:,1);
CT1=CT0; CT1F=CT0F;
% if SNP is on R strand,context should also be reverse complemented
CT1(kFLIP1,:)=tr(fliplr(CT1(kFLIP1,:)),'ACGTacgt','TGCAtgca');
CT1F(kFLIP1,:)=tr(fliplr(CT1F(kFLIP1,:)),'ACGTacgt','TGCAtgca');

C.CT=cellstr([CT1(:,1) repmat('-',C.N,1) FT(:,1) repmat('>?',C.N,1) repmat('-',C.N,1) CT1(:,2)] );
C.T=tabulate(C.CT);

% SNPs in MAF
% SNP from->to convention (flip  such that from is always C or T)
if ~isfield(X,'Reference_Allele') && isfield(X,'ref_allele'), X = rename_field(X,'ref_allele','Reference_Allele'); end
if ~isfield(X,'Tumor_Seq_Allele1') && isfield(X,'newbase'), X = rename_field(X,'newbase','Tumor_Seq_Allele1'); end
if (all(strcmp(X.Reference_Allele,X.Tumor_Seq_Allele1)))
    display('Tumor_Seq_Allele1 is set to Reference_Allele, try Tumor_Seq_Allele2')
    X.Tumor_Seq_Allele1=X.Tumor_Seq_Allele2;
    if (all(strcmp(X.Reference_Allele,X.Tumor_Seq_Allele1)))
        error('no alternate alleles in Tumor_Seq_Allele1 or Tumor_Seq_Allele2')
    end
end

FT=[char(X.Reference_Allele) char(X.Tumor_Seq_Allele1)];
kFLIP=find( (FT(:,1)=='G')|(FT(:,1)=='T'));
FT(kFLIP,:)=tr(FT(kFLIP,:),'ACGTacgt','TGCAtgca');

CT=C.name(X.context65);
CT0=char(CT); CT0=CT0(:,[6 8]);
CT0F=char(CT); CT0F=CT0F(:,1);
CT1=CT0;
CT1F=CT0F;

CT1(kFLIP,:)=tr(fliplr(CT1(kFLIP,:)),'ACGTacgt','TGCAtgca');
CT1F(kFLIP,:)=tr(fliplr(CT1F(kFLIP,:)),'ACGTacgt','TGCAtgca');


X.CT=cellstr([CT1(:,1) repmat('-',X.N,1) FT(:,1) repmat('>',X.N,1) FT(:,2) repmat('-',X.N,1) CT1(:,2)] );
T=tabulate(X.CT);

colors = [1 1 0;0 0.7 0.7;1 0 0;0.1 0.8 0.1;0 0.2 0.8;0.5 0.3 0.7;];

nn=cell2mat(T(:,2));
cc=T(:,1);
XYZ=[];
XY=[];
base={'T';'C';'A';'G'};
snp={'C>T';'C>A';'C>G';'A>G';'A>C';'A>T'};
%baseMap = containers.Map(base,1:4) ;
%snpMap = containers.Map(snp,1:6) ;
XYZ.n=zeros(4,4,6);
XY.n=zeros(8,12);

% maf
for b1=1:4, for b2=1:4, for s=1:6
            c=cellstr([char(base(b1)) '-' char(snp(s)) '-' char(base(b2))]);
            XYZ.cat(b1,b2,s)=c;
            XYZ.col(b1,b2,s,:)=colors(s,:);
            %k=strfindk(cc,c{1});
            k=find(ismember(cc,c{1}));
            y=b2+4*mod(7-s-1,3);
            x=b1+4*floor(s/4);
            if (~isempty(k))
                XYZ.n(b1,b2,s)=nn(k);
                XY.n(x,y)=nn(k);
            end
            XY.col(x,y,:)=colors(s,:);
            XY.cat(x,y)=c;
        end,end,end

% coverage
for b1=1:4, for b2=1:4, for s=1:6
            c=cellstr([char(base(b1)) '-' char(cellfun(@(x) x(1),snp(s))) '>?-' char(base(b2))]);
            %k=strfindk(C.CT,c{1});
            k=find(ismember(C.CT,c{1}));
            y=b2+4*mod(7-s-1,3);
            x=b1+4*floor(s/4);
            if (~isempty(k))
                XY.ncov(x,y)=sum(C.Ncov(k));
            end
        end,end,end

if (nargout<1)
    bar3_with_colors(XY.n,XY.col)
end

function [context65]=ref_context_to_context65(ref_context,C65)
        if (nargin~=2)
            fprintf('*** ref_context_to_context65 needs 2 arguments')
            context65=-1
            return;
        end
        if length(C65.name)~=65
            fprintf('*** ref_context_to_context65 C65 need 65 categories')
            context65=-2
            return;
        end
        
        
        
        XABC=cellfun(@(x) x(10:12),upper(ref_context),'UniformOutput',false);
        
        abc=char(C65.name);
        ABC=abc(:,6:8); ABC(:,2)=abc(:,1);
        context65=NaN*zeros(size(ref_context));
        for i=1:65
            k=find(ismember(XABC,ABC(find(C65.id==i),:)));
            context65(k)=i;
        end

function [Y]=round_collapse_categories_to_65(X)
    if ndims(X)>2, error('doesn''t work with 3D matrices'); end
    
    k = 1+floor( ((1:845)'-1)/13);
    
    if any(isnan(k)) || any(k<1) || any(k>65) || any(k~=round(k))
        error('problem with round_collapse_categories_to_65');
    end
    if length(k)~=size(X,1), error('X has %d rows, whereas %s has %d categories',size(X,1),categs_txt,length(k)); end
    
    Y = zeros(65,size(X,2));
    for i=1:length(k)
        Y(k(i),:) = Y(k(i),:) + X(i,:);
    end



function test0
SET='AN_TCGA_THCA_22Feb2012_TB'
[XY,X,C1]=MutationSpectrumCategs
X0=X;

function test1
clf
[XY]=MutationSpectrumCategs(X,C1)
bar3_with_colors(XY.n,XY.col)
bar3_with_colors(XY.n./XY.ncov,XY.col)

bar3_with_colors(XY.ncov,XY.col)


function makeContext65;

% /xchip/cga1/lawrence/db/hg19/c65e/categs.txt
catfile='~/Projects/Cancer/tools/matlab/categs.txt'
k = map_categories_to_65(catfile);
categs = load_struct(catfile);
categs.longname=categs.name;
q=regexp(categs.name,':','split');
categs.name=cellfun(@(x) x(1),q);
[context.name i]=unique(categs.name);
context.id=k(i);
context=trimStruct(context,sort(context.id))
f='categs65.txt';
printStruct(context,[],f)




function N = average_exome_coverage
N=[...
    516645      336113      545347      365625      462172      488457      797766      489731
    549735      357085      540864      397040      230606      269162      226201      252346
    469526      429764      160802      389408      664422      477669      235217      616696
    509693      512648      186053      541099      537050      568802      169936      589050
    577447      538629      601391      389620      168530      186767      234416      162102
    556167      509431      469219      426880      536918      514744      661776      473758
    252140      397574      486994      367907      226652      552864      803001      556516
    267878      359998      486805      338852      231118      560177      467349      524231
    ];
N=N';
N=N(:);

function N = average_genome_coverage
N=[...
    94564228    36158942    49494329    61919448    46789183    35927089    48482195    45521663
    48882158    22927567    40130275    33037787    51852385    28182087    32110875    50214797
    49545058    27998332     5899505    39750998    44139802    30001011     5654254    42486275
    34775000    27270757     4806689    33579428    48349332    37043310     5098333    54501232
    54403725    33559069    42441923    39792917     5088120     4802745     5651927     5914345
    37066788    27263586    30032551    28031174    48352258    34807746    44225804    49728134
    50164185    33013168    45519675    61986409    32073940    40143640    48528214    49570546
    28197095    22968750    36009505    36288035    51937650    48943412    46991718    94906144
    ];
N=N';
N=N(:);


