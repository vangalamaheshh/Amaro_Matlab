function [X]=alleleFreqPurity(f)

if ischar(f)
    M=load_table(f,char(9));
else
    M=f;
end
% chromosomes 1-22
if (iscell(M.Chromosome))
    a=chromosome2num(M.Chromosome);
else
    a=M.Chromosome;
end
M=trimStruct(M,find(a<23));

M=trimStruct(M,strfindk(M.Variant_Type,'SNP'));
x=M.t_alt_count;
if (iscell(x))
    M.t_alt_count=str2num(char(M.t_alt_count));
    M.t_ref_count=str2num(char(M.t_ref_count));
end

% samples
sample=unique(M.Tumor_Sample_Barcode);


% allele fraction bins
afb=0:0.005:1;

% output struct
X.sample=sample;
X.afb=repmat(afb,length(sample),1);
X.af=0*X.afb;
X.afhist=0*X.afb;
X.afbeta=0*X.afb;
X.afp=NaN*zeros(length(sample),1);
X.pur=NaN*zeros(length(sample),1);
X.Nmut=zeros(length(sample),1);
X.afCI=NaN*zeros(length(sample),2);
M=trimStruct(M,find(M.i_t_lod_fstar>8));
for i=1:length(X.sample)
    s1=X.sample{i};
    % mutations in this sample
    ks=strfindk(M.Tumor_Sample_Barcode,s1);
    
    if (isempty(ks)),continue; end
   
    % allele fraction histogram (no smoothing)
    X.afhist(i,:)=hist(M.t_alt_count(ks)./(M.t_alt_count(ks)+M.t_ref_count(ks)),afb);

    % binomial smoothing
    for k1=ks'
        %binomial pdf for each mutation
        x1=binopdf(M.t_alt_count(k1),M.t_alt_count(k1)+M.t_ref_count(k1),afb);
        % add each mutation unit normalized (binopdf is normalized to k not p) 
        X.af(i,:)=X.af(i,:)+(x1/sum(x1));
        
        %beta pdf for each mutation
        x2=betapdf(afb,M.t_alt_count(k1)+1,M.t_ref_count(k1)+1);
        % add each mutation
        X.afbeta(i,:)=X.afbeta(i,:)+(x2/sum(x2));
        
    end
    % this sample's af dist
    x1=X.af(i,:);
        
    % find het peak below 0.5 
    j=2:(length(x1)-1);
    kp=find((x1(j-1)<=x1(j))&(x1(j+1)<x1(j)));    
    xp=x1(kp+1);
    ap=afb(kp+1);
    ap=ap(ap<0.55);
    xp=xp(ap<0.55);
    ap=ap(xp>0.000025);
    xp=xp(xp>0.000025);
    [ap,k]=max(ap);
    ap=min([0.5 ap]);
    if (~isempty(k))
        X.afp(i)=ap;
        X.pur(i)=2*X.afp(i);    
        X.Nmut(i)=length(ks);
        % CI on purity (assuming it's the het peak)
        laf=log(X.af(i,:)+1e-10);
        b0=find(afb==ap);
        lafT=laf(b0)-1;
        for b=b0:length(afb)
            if (laf(b)<=(lafT))
                X.purCI(i,2)=min([1 2*(afb(b))]);
                break;
            end
        end
        for b=b0:-1:1
            if (laf(b)<=(lafT))
                X.purCI(i,1)=2*(afb(b));
                break;
            end
        end
    end
    
end


function test
%%
close all
clear
%path(path,'/Users/stewart/CancerGenomeAnalysis/trunk/matlab')
%path(path,'/Users/stewart/CancerGenomeAnalysis/trunk/matlab/seq')
% -------------------------------------------------------------------------
% maf file name goes here
f='~/Projects/Cancer/CLL/TRIOS_II/CLL_Trios_II.maf.txt'
f='~/Projects/Cancer/Pediatric/EwingsSarcoma/EwingsSarcoma_30Oct2011.maf.txt'
f='/Users/stewart/Projects/Cancer/Pediatric/EwingsSarcoma/WGS/muTect/EwngSRC-SJDES006-Tumor.maf'
f='/Users/stewart/Projects/Cancer/Pediatric/EwingsSarcoma/WGS/muTect/EwingsSarcoma_WGS.maf'
f='/Users/stewart/Projects/Cancer/BLCA/PR_TCGA_BLCA_Capture-59.final_analysis_set.maf'
f='/Users/stewart/Projects/Cancer/BLCA/PR_TCGA_BLCA_Capture-59.final_analysis_set.maf'
f='/local/cga-fh/cga/Sigma_Pediatric/Individual_Set/AN_Sigma_EwingSarcoma_29Jan2013_TN26/jobs/mutsig1.5_Sigma_Pediatric/AN_Sigma_EwingSarcoma_29Jan2013_TN26.final_analysis_set.maf'
f='/Users/stewart/Projects/Cancer/Pediatric/EwingsSarcoma/AN_Sigma_Ewings_Sarcoma_27Feb2013_32TN.combined.capture.maf'
% -------------------------------------------------------------------------

[X]=alleleFreqPurity(f)

clf
i=0
for s=1:length(X.sample)
    i=i+1;
    if (i>6)
        figure
        i=1; 
    end
    ksub=find(X.afb(s,:)<(X.afp(s)*0.75));
    X.Nsubcl(s,1)=sum(X.af(s,ksub))
    subplot(3,2,i)
    semilogy(X.afb(s,:),X.af(s,:),'b-')
    ylim([.001,10]);
    line(X.afp(s)+[0 0],[.001 max(X.af(s,:))*1.1],'color','r','linestyle','--')
    %title(X.sample{s});
    sam=regexprep(X.sample{s},'-Tumor',''); 
    text(0.5, 0.9, sam,'units','normalized','interpreter','none');
    text(0.5, 0.8, sprintf('n=%d',X.Nmut(s)),'units','normalized','interpreter','none');
    xlabel('allele fraction ')
    grid on;
end


printStruct(X)


sum(X.Nsubcl)/sum(X.Nmut)

%%


