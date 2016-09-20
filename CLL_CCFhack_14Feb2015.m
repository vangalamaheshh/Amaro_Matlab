% clear
% % X chrom
% XX=load_table('~/Documents/FullCll8chrX_pre_abs.maf')
% % ccfs for each mutation
% X=load_table('~/Documents/ABSOLUTE_maf.all.Stilgenbauer_Tps.DFCI.ICGC.maf')
% % gender
% G=load_table('/Users/amaro/Documents/GenderAssignmentFullset.txt')
% % ABSOLUTE purity/ploidy
% A=load_table('~/Documents/FullSetAllCLLSamplesForCLL8.txt')
% % ABSOLUTE seg tab
% T=load_table('~/Documents/MergedABS_DFCI_Stilgenbauer_ICGC.seg')
% X.ABSOLUTE=ones(size(X.Start_position));
% X.Tumor_Sample_Barcode=X.sample;
% X.t_alt_count=X.alt;
% X.t_ref_count=X.ref;
% XX.ABSOLUTE=zeros(size(XX.Start_position));
% X=mergeStruct(X,XX);
% save(['/Users/amaro/Documents/CLL_Stilgenbauer.mat'],'X','XX','G','A','T')
% 
% 
% 
% %%
% % clear
% % cd /Users/stewart/Projects/Cancer/CLL/CLL8
% load(['/Users/stewart/Projects/Cancer/CLL/CLL8/CLL_Stilgenbauer.mat'])
case_control_table=load_struct('/Volumes/xchip_cga_home/amaro/CLL/StilgenbauerMafFreeze2.0/case_control_table.tsv');
case_control_table_stil=load_struct('/Users/amaro/Documents/StilgenbauerABSresults/case_control_tableStil.tsv');
case_control_table=mergeStruct(case_control_table,case_control_table_stil);
case_control_table=rmfield(case_control_table,'N');
[i m]=ismember(case_control_table.pair_id,T.sample);
T.sample(m(m>0))=case_control_table.case_sample(i);

sample=sort(A.sample)
%X.sample=X.Tumor_Sample_Barcode;
%X.tid=cellfun(@(x) x{1},regexp(X.Tumor_Sample_Barcode,'-Tumor','split'),'UniformOutput',false)
X.tid=upper(X.sample);
A.tid=upper(A.sample);
T.tid=upper(T.sample);
% A.tid=upper(regexprep(A.sample,'-TP-','-TUMOR-'));
% A.tid=upper(regexprep(A.tid,'-NT',''));
% A.tid=upper(regexprep(A.tid,'-NUNK',''));
% A.tid=upper(regexprep(A.tid,'-SM-.....$',''));
tabulate(A.tid)

% T.tid=upper(regexprep(T.sample,'-TP-','-TUMOR-'));
% T.tid=upper(regexprep(T.tid,'-NT',''));
% T.tid=upper(regexprep(T.tid,'-NUNK',''));
% T.tid=upper(regexprep(T.tid,'-SM-.....$',''));
% tabulate(T.tid)
% tid=A.tid;

% check matching tid's (missing tid's should be empty) 
[i m]=ismember(X.tid,A.tid);
tabulate(sort(X.tid(~i)))
[i m]=ismember(A.tid,X.tid);
tabulate(sort(A.tid(~i)))
[i m]=ismember(T.tid,A.tid);
tabulate(sort(T.tid(~i)))
[i m]=ismember(A.tid,T.tid);
tabulate(sort(A.tid(~i)))

% ok

X.x=xhg19(X.Chromosome,X.Start_position);
T.x1=xhg19(T.Chromosome,T.Startbp);
T.x2=xhg19(T.Chromosome,T.Endbp);
N=length(A.tid)
X.corrected_total_cn=NaN*X.x;
X.modal_total_cn=NaN*X.x;
X.LOH=NaN*X.x;
X.modal_a1=NaN*X.x;
X.modal_a2=NaN*X.x;
X.cancer_cell_frac_a1=NaN*X.x;
X.cancer_cell_frac_a2=NaN*X.x;
%X.cancer_cell_frac=NaN*X.x;
X.ccf_hat_hack=NaN*X.x;
X.ccf_CI95_low_hack=NaN*X.x;
X.ccf_CI95_high_hack=NaN*X.x;
% af binning
daf=0.01;
% af range
aft=0:daf:1

i=1
X1=[]

for i=1:N
    s1=sample{i}
    tid1=A.tid{i};
    x1=trimStruct(X,ismember(X.tid,tid1));
    t1=trimStruct(T,ismember(T.tid,tid1));
    % CN at each mutation site
    [ma,mb,dab,fa,nb]=matchI(x1.x,t1,100,0);
    %[chrom2num(x1.Chromosome(ma>0)) t1.Chromosome(ma(ma>0))]
    i1=find(ma>0);
    m1=ma(i1);
    % seg_tab info 
    x1.corrected_total_cn(i1)=t1.corrected_total_cn(m1);
    x1.modal_total_cn(i1)=t1.modal_total_cn(m1);
    x1.LOH(i1)=(t1.LOH(m1));
    x1.modal_a1(i1)=(t1.modala1(m1));
    x1.modal_a2(i1)=(t1.modala2(m1));
    x1.cancer_cell_frac_a1(i1)=(t1.cancercellfraca1(m1));
    x1.cancer_cell_frac_a2(i1)=(t1.cancercellfraca2(m1));

        
        
    % purity
    x1.pur=A.purity(ismember(A.tid,tid1))+0*x1.Start_position;
    %x1.dna_fraction_in_tumor=(x1.pur.*x1.corrected_total_cn)./(x1.pur.*x1.corrected_total_cn+2*(1-x1.pur));
    x1.dna_fraction_in_tumor=(x1.pur.*x1.modal_total_cn)./(x1.pur.*x1.modal_total_cn+2*(1-x1.pur));
    x1.dna_fraction_in_tumor(x1.corrected_total_cn==0)=NaN;
    x1.ccf_CI95_low_hack=NaN*x1.Start_position;
    x1.ccf_CI95_high_hack=NaN*x1.Start_position;
    %x1.cancer_cell_frac=NaN*x1.Start_position;
    % initial ccf guess 
    x1.ccf_guess=x1.i_tumor_f./(x1.pur.*x1.modal_total_cn);
    x1.ccf_hat_hack=NaN*x1.Start_position;
    x1.ccf_mode_hack=NaN*x1.Start_position;
    x1.i_tumor_f_corrected=NaN*x1.Start_position;

    j=1;
    % loop over each mutation 
    for j=1:length(x1.pur)
        x1.paf(j,:)=daf*betapdf(aft,x1.t_alt_count(j)+1,x1.t_ref_count(j)+1);
        % x1.paf(j,:) not always normalized (t_alt_count = 0 problem) 
        x1.paf(j,:)=x1.paf(j,:)/sum(x1.paf(j,:));
        x1.pccf(j,:)=NaN*x1.paf(j,:);
        x1.paft(j,:)=NaN*x1.paf(j,:);
        % skip events that ABSOLUTE skipped
        if isnan(x1.dna_fraction_in_tumor(j)), continue; end
        if x1.modal_total_cn(j)<1, continue; end
        % possible ref count in tumor
        nreft=0:x1.t_ref_count(j);
        % likelihood for each tumor ref count
        wreft=binopdf(nreft+x1.t_alt_count(j),x1.t_ref_count(j)+x1.t_alt_count(j),x1.dna_fraction_in_tumor(j));
        % init tumor specific allele fraction dist
        paft=0*aft; paftv=[];
        % loop over each tumor ref count
        for n=1:length(nreft)
            paft1=daf*wreft(n)*betapdf(aft,x1.t_alt_count(j)+1,nreft(n)+1);
            paftv=[paftv; paft1];
            paft=paft+paft1;
        end
        % store tumor allele fracton pdf to maf
        x1.paft(j,:)=paft/sum(paft);
        % possible multiplicities depending on when the mutation happened relative to the SCNA 
        qm=unique(sort([1 x1.modal_a1(j) x1.modal_a2(j)]));
        qm=qm(qm>0);
        % equal weight for each qm - 
        wqm=(1/length(qm))+0*qm;
        
        %ccf=0:daf:1;
        ccf=0:daf:1;
        pccf=0*ccf;
        pccfv=[];
        m=1;
        for m=1:length(qm)
            % scale ccf down by n
            ccfb=x1.modal_total_cn(j)*(ccf/qm(m));
            % project aft into ccfb regions with ccf bins
            p1=histw(x1.paft(j,:),ccfb,ccf);
            % weight by wqm
            pccf1=wqm(m)*p1;
            pccf=pccf+pccf1;
            pccfv=[pccfv; pccf1];
        end
        %stairs(ccf,[pccfv; pccf]');
        %line(ccf,x1.paft(j,:),'linestyle','--','color',0.5*[1 1 1]);
        pccfx=1-sum(pccf);
        %pccf=pccf/sum(pccf);
        x1.pccf(j,:)=pccf;
        x1.pccf(j,end)=x1.pccf(j,end)+1-sum(pccf);
        if (x1.pccf(j,end)<0), x1.pccf(j,end)=0;end
        if isnan(sum(pccf)), continue; end
        %pause()
        cccf=cumsum(pccf);
        ilow=sum(cccf<0.025);
        ihigh=sum(cccf<0.975);
        [pmax,imax]=max(pccf);
        if (ilow<1), ilow=1;end
        if (ihigh<=ilow), ihigh=ilow+1;end
        if (ihigh>x1.N), ihigh=x1.N;end
        x1.ccf_CI95_low_hack(j,1)=ccf(ilow);
        x1.ccf_CI95_high_hack(j,1)=ccf(ihigh);
        %x1.cancer_cell_frac(j,1)=ccf(imax);
        x1.ccf_mode_hack(j,1)=ccf(imax(1));
        x1.ccf_hat_hack(j,1)=sum(ccf.*pccf)./sum(pccf);
        [q k]=max(x1.paft(j,:));
        x1.i_tumor_f_corrected(j,1)=aft(k);
    end
    if isempty(X1)
        X1=x1;
    else
        X1=mergeStruct(X1,x1);
    end
    
end

%save(['/xchip/cga_home/adunford/DLBCL/CCF_MAF/DLBCL_Pairs/DLBCL_Pairs.CCF.' TODAY 'a.maf.mat'],'-struct','X1')
save(['CLL_Stilgenbauer.a.maf.mat'],'-struct','X1')
printStruct(X1,-1,['CLL_Stilgenbauer.a.maf'])

P.XLAB='ABS CCF'; P.YLAB='HACK CCF';
hist2(X1.ccf_hat,X1.ccf_hat_hack,100,100,[],P)

%%



%%
% catenate segtabs
a='/xchip/cga_home/amaro/CLL/ABSOLUTE/ABSOLUTE_results/Redo_2.11.CLL_Stilgenbauer/reviewed/SEG_MAF/'   
f1=rdir([a '*.segtab.txt'])
T=[]
for i=1:length(f1)
    x1=load_table(f1(i).name);
    [ i sum(x1.End_bp-x1.Start_bp)]
    if isempty(T)
        T=x1;
    else
        T=mergeStruct(T,x1);
    end
end
printStruct(T,-1,['/xchip/cga_home/adunford/CLL/CCF_MAF/CLL_Stilgenbauer.' TODAY '.segtab.txt'])

% catenate segtabs
a='/xchip/cga_home/amaro/CLL/ABSOLUTE/ABSOLUTE_results/Redo_2.11.CLL_Stilgenbauer/reviewed/SEG_MAF/'   
f1=rdir([a '*ABS_MAF.txt'])
T=[]
for i=1:length(f1)
    x1=load_table(f1(i).name);
    i
    if isempty(T)
        T=x1;
    else
        T=mergeStruct(T,x1);
    end
end
printStruct(T,-1,['/xchip/cga_home/adunford/CLL/CCF_MAF/CLL_Stilgenbauer.' TODAY '.ABS_MAF.txt'])
