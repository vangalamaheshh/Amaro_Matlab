clear
% X chrom
XX=load_tsv('/Users/amaro/Documents/StilgenbauerABSresults/FullStilgenbauerMafchrX.maf')
% ccfs for each mutation
X=load_tsv('/Users/amaro/Documents/StilgenbauerABSresults/')
% gender
G=load_table('/Users/stewart/Projects/Cancer/CLL/CLL8/XchrData_gender_assignment.txt')
% ABSOLUTE purity/ploidy
A=load_table('/Users/stewart/Projects/Cancer/CLL/CLL8/Final_Stilgenbauer_2.4.ATW.ABSOLUTE.table.txt')
% ABSOLUTE seg tab
T=load_table('/Users/stewart/Projects/Cancer/CLL/CLL8/CLL_Stilgenbauer.14Feb2015.segtab.txt')
X.ABSOLUTE=ones(size(X.Start_position));
X.Tumor_Sample_Barcode=X.sample;
X.t_alt_count=X.alt;
X.t_ref_count=X.ref;
XX.ABSOLUTE=zeros(size(XX.Start_position));
X=mergeStruct(X,XX);
save(['/Users/stewart/Projects/Cancer/CLL/CLL8/CLL_Stilgenbauer.' TODAY '.mat'],'X','XX','G','A','T')



%%
clear
cd /Users/stewart/Projects/Cancer/CLL/CLL8
load(['/Users/stewart/Projects/Cancer/CLL/CLL8/CLL_Stilgenbauer.14Feb2015.mat'])

sample=sort(A.sample)
%X.Tumor_Sample_Barcode=X.sample;
%X.tid=cellfun(@(x) x{1},regexp(X.Tumor_Sample_Barcode,'-Tumor','split'),'UniformOutput',false)
X.tid=upper(X.Tumor_Sample_Barcode)
tabulate(X.tid)
X.tid(ismember(X.Chromosome,'X'))

A.tid=upper(regexprep(A.sample,'-TP-','-TUMOR-'));
A.tid=upper(regexprep(A.tid,'-NT',''));
A.tid=upper(regexprep(A.tid,'-NUNK',''));
A.tid=upper(regexprep(A.tid,'-SM-.....$',''));
tabulate(A.tid)

T.tid=upper(regexprep(T.sample,'-TP-','-TUMOR-'));
T.tid=upper(regexprep(T.tid,'-NT',''));
T.tid=upper(regexprep(T.tid,'-NUNK',''));
T.tid=upper(regexprep(T.tid,'-SM-.....$',''));
tabulate(T.tid)

G.tid=upper(regexprep(G.pair_id,'-TP-','-TUMOR-'));
G.tid=upper(regexprep(G.tid,'-NT',''));
G.tid=upper(regexprep(G.tid,'-NUNK',''));
G.tid=upper(regexprep(G.tid,'-SM-.....$',''));
tabulate(G.tid)


tid=A.tid;

% check matching tid's (missing tid's should be empty) 
[i m]=ismember(X.tid,A.tid);
tabulate(sort(X.tid(~i)))
[i m]=ismember(A.tid,X.tid);
tabulate(sort(A.tid(~i)))
[i m]=ismember(A.tid,T.tid);
tabulate(sort(A.tid(~i)))
[i m]=ismember(T.tid,A.tid);
tabulate(sort(T.tid(~i)))

% ok
X.x=xhg19(X.Chromosome,X.Start_position);
T.x1=xhg19(T.Chromosome,T.Start_bp);
T.x2=xhg19(T.Chromosome,T.End_bp);
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
daf=0.001;
% af range
aft=(daf/2):daf:1

i=1
X1=[]

for i=1:N
    s1=sample{i}
    tid1=A.tid{i};
    x1=trimStruct(X,ismember(X.tid,tid1));
    t1=trimStruct(T,ismember(T.tid,tid1));
    g1=trimStruct(G,ismember(G.tid,tid1));
    % CN at each mutation site
    [ma,mb,dab,fa,nb]=matchI(x1.x,t1,100,0);
    %[chrom2num(x1.Chromosome(ma>0)) t1.Chromosome(ma(ma>0))]
    i1=find(ma>0);
    m1=ma(i1);
    % seg_tab info 
    x1.corrected_total_cn(i1)=t1.corrected_total_cn(m1);
    x1.modal_total_cn(i1)=t1.modal_total_cn(m1);
    x1.LOH(i1)=(t1.LOH(m1));
    x1.modal_a1(i1)=(t1.modal_a1(m1));
    x1.modal_a2(i1)=(t1.modal_a2(m1));
    x1.cancer_cell_frac_a1(i1)=(t1.cancer_cell_frac_a1(m1));
    x1.cancer_cell_frac_a2(i1)=(t1.cancer_cell_frac_a2(m1));

        
    if any(x1.q_hat==0)
        x1.q_hat(x1.q_hat==0)=1;
        x1.HS_q_hat_2(x1.q_hat==0)=1;
    end
        
    if any(ismember(x1.Chromosome,'X'))
        %keyboard
        k=find(ismember(x1.Chromosome,'X'));
        if ismember(g1.gender,'Male')
            x1.q_hat(k)=1;
            x1.HS_q_hat_1(k)=0;
            x1.HS_q_hat_2(k)=1;
        else
            x1.q_hat(k)=2;
            x1.HS_q_hat_1(k)=1;
            x1.HS_q_hat_2(k)=1;
        end
    end
    
    % purity
    x1.pur=A.purity(ismember(A.tid,tid1))+0*x1.Start_position;
    %x1.dna_fraction_in_tumor=(x1.pur.*x1.corrected_total_cn)./(x1.pur.*x1.corrected_total_cn+2*(1-x1.pur));
    %x1.dna_fraction_in_tumor=(x1.pur.*x1.modal_total_cn)./(x1.pur.*x1.modal_total_cn+2*(1-x1.pur));
    x1.dna_fraction_in_tumor=(x1.pur.*x1.q_hat)./(x1.pur.*x1.q_hat+2*(1-x1.pur));
    %x1.dna_fraction_in_tumor(x1.corrected_total_cn==0)=NaN;
    x1.dna_fraction_in_tumor(x1.q_hat==0)=NaN;
    x1.ccf_CI95_low_hack=NaN*x1.Start_position;
    x1.ccf_CI95_high_hack=NaN*x1.Start_position;
    %x1.cancer_cell_frac=NaN*x1.Start_position;
    % initial ccf guess 
    x1.ccf_guess=x1.i_tumor_f./(x1.pur.*x1.modal_total_cn);
    x1.ccf_guess=x1.i_tumor_f./(x1.pur.*x1.q_hat);
    x1.ccf_mode_hack=NaN*x1.Start_position;
    x1.ccf_mean_hack=NaN*x1.Start_position;
    x1.ccf_median_hack=NaN*x1.Start_position;
    x1.i_tumor_f_corrected=NaN*x1.Start_position;
    x1.clonal_hack = NaN*x1.Start_position;
    x1.subclonal_hack = NaN*x1.Start_position;

    j=1;
    % loop over each mutation 
    for j=1:length(x1.pur)
        % allel fraction pdf
        x1.paf(j,:)=daf*betapdf(aft,x1.t_alt_count(j)+1,x1.t_ref_count(j)+1);
        % x1.paf(j,:) not always normalized (t_alt_count = 0 problem) 
        if abs(sum(x1.paf(j,:))-1)>0.1,'paf',sum(x1.paf(j,:)), keyboard;        end
        x1.paf(j,:)=x1.paf(j,:)/sum(x1.paf(j,:));
        x1.pccf(j,:)=NaN*x1.paf(j,:);
        x1.paft(j,:)=NaN*x1.paf(j,:);
        % skip events that ABSOLUTE skipped
        if isnan(x1.dna_fraction_in_tumor(j)), continue; end
        % skip regions with no tumor DNA
        %if x1.modal_total_cn(j)<1, continue; end
        if x1.q_hat(j)<1
            'bad q_hat'
            printStruct(x1,j)
            continue; 
        end
        % all possible ref count in tumor from 0 to t_ref_count
        nreft=0:x1.t_ref_count(j);
        % likelihood for each possible tumor ref count
        wreft=binopdf(nreft+x1.t_alt_count(j),x1.t_ref_count(j)+x1.t_alt_count(j),x1.dna_fraction_in_tumor(j));
        % if abs(sum(wreft)-1)>0.1,'wreft',sum(wreft), keyboard;        end
        % normalize wreft within possible nreft counts
        wreft=wreft/sum(wreft);
        % init tumor specific allele fraction dist
        paft=0*aft; paftv=[];
        % loop over each tumor ref count
        for n=1:length(nreft)
            paft1=daf*wreft(n)*betapdf(aft,x1.t_alt_count(j)+1,nreft(n)+1);
            paftv=[paftv; paft1];
            paft=paft+paft1;
        end
        % store tumor allele fracton pdf to maf
        if abs(sum(paft)-1)>0.1,'paft',sum(paft), keyboard;        end
        x1.paft(j,:)=paft/sum(paft);
        % possible multiplicities depending on when the mutation happened relative to the SCNA 
        %qm=unique(sort([1 x1.modal_a1(j) x1.modal_a2(j)]));
        qm=unique(sort([1 x1.HS_q_hat_1(j) x1.HS_q_hat_2(j)]));           
        qm=qm(qm>0);
        % equal weight for each qm - 
        wqm=(1/length(qm))+0*qm;
        
        %ccf=0:daf:1;
        ccf=(daf/2):daf:1;
        pccf=0*ccf;
        pccfv=[];
        m=1;
        for m=1:length(qm)
            % scale ccf down by n
            %ccfb=x1.modal_total_cn(j)*(ccf/qm(m));
            ccfb=x1.q_hat(j)*(ccf/qm(m));
            % project aft into ccfb regions with ccf bins
            p1=histw(x1.paft(j,:),ccfb,ccf);
            %keyboard;  
            % weight by wqm
            pccf1=wqm(m)*p1;
            pccf=pccf+pccf1;
            pccfv=[pccfv; pccf1];
        end
        stairs(ccf,[pccfv; pccf]');
        line(ccf,x1.paft(j,:),'linestyle','--','color',0.5*[1 1 1]);
        pccfx=1-sum(pccf);
        %pccf=pccf/sum(pccf);
        x1.pccf(j,:)=pccf;
        x1.pccf(j,end)=x1.pccf(j,end)+1-sum(pccf);
        if (x1.pccf(j,end)<0), x1.pccf(j,end)=0;end
        if isnan(sum(pccf)), continue; end
        %pause()
        cccf=cumsum(x1.pccf(j,:));
        ilow=sum(cccf<0.025);
        ihigh=sum(cccf<0.975);
        imedian=sum(cccf<0.5);
        if (imedian<1), imedian=1; end
        [pmax,imax]=max(pccf);
        if (ilow<1), ilow=1;end
        if (ihigh<=ilow), ihigh=ilow+1;end
        if (ihigh>x1.N), ihigh=x1.N;end
        x1.ccf_CI95_low_hack(j,1)=ccf(ilow);
        x1.ccf_CI95_high_hack(j,1)=ccf(ihigh);
        %x1.cancer_cell_frac(j,1)=ccf(imax);
        x1.ccf_mode_hack(j,1)=ccf(imax(1));
        x1.ccf_mean_hack(j,1)=sum(ccf.*x1.pccf(j,:))./sum(x1.pccf(j,:));
        x1.ccf_median_hack(j,1)=ccf(imedian(1));
        [q k]=max(x1.paft(j,:));
        x1.i_tumor_f_corrected(j,1)=aft(k);
        
        % set clonal flags
        x1.clonal_hack(j) = x1.ccf_mode_hack(j)>0.99;
        x1.subclonal_hack(j) = x1.ccf_mode_hack(j)<0.99;

    end
    if isempty(X1)
        X1=x1;
    else
        X1=mergeStruct(X1,x1);
    end
    
end

ff=fieldnames(X1)
ff=ff(strfindk(ff,'v_'));n0=length(ff)
X1.ccf_ABS=0:(1/(n0-1)):1;
X1.pccf_ABS=NaN(X1.N,n0)
for j=1:n0
    X1.pccf_ABS(:,j)=X1.(ff{j});
end


%save(['/xchip/cga_home/adunford/DLBCL/CCF_MAF/DLBCL_Pairs/DLBCL_Pairs.CCF.' TODAY 'a.maf.mat'],'-struct','X1')
save(['CLL_Stilgenbauer.' TODAY 'a.maf.mat'],'-struct','X1')
printStruct(X1,-1,['CLL_Stilgenbauer.' TODAY 'a.maf'])


P.XLAB='ABS CCF'; P.YLAB='HACK CCF mode';
hist2(X1.ccf_hat,X1.ccf_mode_hack,100,100,[],P)

feps=['CLL_Stilgenbauer.ccf_hat_hack_mode_scatterplot.' TODAY '.eps']
saveas(gcf,feps,'epsc')


% paper recipe for clonal_ix
qA = sum(X1.pccf_ABS(:,95:end),2)>0.5;
tabulate(qA)

% hack recipe for clonal_ix
q=zeros(size(X1.clonal_hack));
q(X1.subclonal_hack>0)=-1;
q(X1.clonal_hack>0)=1;
tabulate(q)

%%

P.XLAB='ABS CCF'; P.YLAB='HACK CCF median';
hist2(X1.ccf_hat,X1.ccf_median_hack,100,100,[],P)

P.XLAB='ABS CCF'; P.YLAB='HACK CCF mean';
hist2(X1.ccf_hat,X1.ccf_mean_hack,100,100,[],P)

n1=length(ccf); n2=n1/10; x=mean(reshape(ccf,10,n2),1)';
for i=1:X1.N
  X1.pccf_HACK(i,:)=sum(reshape(X1.pccf(i,:),10,n2),1)'; 
end

k=1:300
figure(1)
cmap=jet(1023); cmap(1,:)=[1 1 1];    
imagesc(-log(X1.pccf_ABS(k,:))); colormap(cmap)
figure(2)
imagesc(X1.pccf_HACK(k,:));colormap(cmap)

clf
t=regexprep(strcat(X1.tid,'--',X1.Chromosome,':',cellstr(num2str(X1.Start_position)),'--',cellstr(num2str(X1.i_tumor_f))), ' ','')
for i=1:X1.N  
  plot(X1.ccf_ABS, X1.pccf_ABS(i,:),x,  X1.pccf_HACK(i,:) )
  title(t{i})
  xlabel('CCF')
  legend({'ABS','HACK'})
  text(0.1,0.9,sprintf('ABS  ccf_hat=%.2f',X1.ccf_hat(i)),'units','normalized')
  text(0.1,0.85,sprintf('HACK ccf_mode=%.2f',X1.ccf_mode_hack(i)),'units','normalized')
  pause()
end


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
