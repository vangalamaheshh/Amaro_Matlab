% % % clear
% % % % X chrom
% %  XX=load_struct('/Volumes/xchip_cga_home/amaro/CLL/ABSOLUTE/FullSetchrX.maf');
% % % % ccfs for each mutation
% %  X=load_struct('/Users/amaro/FullSetCLL8.CCF_Clustered.ForceCall.maf');
% % % % gender
% %  G=load_table('/Volumes/xchip_cga_home/amaro/CLL/HetsOnX.txt');
% % % % ABSOLUTE purity/ploidy
% %  A=load_table('/Volumes/xchip_cga_home/amaro/CLL/ABSOLUTE/FullSetABSOLUTE.Table.CLL8.txt');
% % % % ABSOLUTE seg tab
% %  T=load_table('/Volumes/xchip_cga_home/amaro/CLL/ABSOLUTE/FullSet.CLL8.ABSOLUTECopy.seg');
% %  X.ABSOLUTE=ones(size(X.Start_position));
% %  X.Tumor_Sample_Barcode=X.sample;
% %  X.t_alt_count=X.alt;
% %  X.t_ref_count=X.ref;
% %  XX.ABSOLUTE=zeros(size(XX.Start_position));
% %  X=mergeStruct(X,XX);
% %  save(['/Users/amaro/Documents/CLL_Stilgenbauer.3.5.mat'],'X','XX','G','A','T')
% % 
% 
% 
% %%
% clear
% %cd /Users/stewart/Projects/Cancer/CLL/CLL8
% %load(['/Users/stewart/Projects/Cancer/CLL/CLL8/CLL_Stilgenbauer.14Feb2015.mat'])
% load(['/Users/amaro/Documents/CLL_Stilgenbauer.3.5.mat'])
% P=load_struct('/Volumes/xchip_cga_home/amaro/CLL/ABSOLUTE/case_control_individual_table.fullset.txt')
% 
% %X.Tumor_Sample_Barcode=X.sample;
% %X.tid=upper(X.Tumor_Sample_Barcode)
% %X.tid=cellfun(@(x) x{1},regexp(X.Tumor_Sample_Barcode,'-Tumor','split'),'UniformOutput',false)
% % tabulate(X.tid)
% % X.tid(ismember(X.Chromosome,'X'))
% % tabulate(A.sample)
% % 
% % A.tid=upper(A.sample)
% % tabulate(A.tid)
% % 
% % %T.sample has many pair_ids
% % [i m]=ismember(T.sample,P.pair_id);
% % T.tid=upper(T.sample);
% % T.tid(i)=upper(P.case_sample(m(i)));
% % tabulate(T.tid)
% % 
% % [i m]=ismember(G.pair_id,P.pair_id);
% % G.tid=upper(G.pair_id);
% % G.tid(i)=upper(P.case_sample(m(i)));
% % tabulate(G.tid)
% % 
% % 
% % 
% % tid=A.tid;
% 
% % check matching tid's (missing tid's should be empty) 
% [i m]=ismember(X.tid,A.tid);
% tabulate(sort(X.tid(~i)))  
% [i m]=ismember(A.tid,X.tid);
% tabulate(sort(A.tid(~i)))
% [i m]=ismember(A.tid,T.tid);
% tabulate(sort(A.tid(~i)))
% [i m]=ismember(T.tid,A.tid);
% tabulate(sort(T.tid(~i)))   
% 
% % ok
% % X.Start_position=str2double(X.Start_position);
% %X.corrected_total_cn=str2double(X.corrected_total_cn);
% X.q_hat=str2double(X.q_hat);
% X.HS_q_hat_1=str2double(X.HS_q_hat_1);
% X.HS_q_hat_2=str2double(X.HS_q_hat_2);
% X.t_alt_count=str2double(X.t_alt_count);
% X.t_ref_count=str2double(X.t_ref_count);
% X.i_tumor_f=X.t_alt_count./(X.t_alt_count+X.t_ref_count);
% %X.LOH=str2double(X.LOH);
% % X.x=xhg19(X.Chromosome,X.Start_position);
% % T.x1=xhg19(T.Chromosome,T.Startbp);
% % T.x2=xhg19(T.Chromosome,T.Endbp);
% N=length(A.tid)
% X.corrected_total_cn=NaN*X.x;
% X.modal_total_cn=NaN*X.x;
% X.LOH=NaN*X.x;
% X.modal_a1=NaN*X.x;
% X.modal_a2=NaN*X.x;
% X.cancer_cell_frac_a1=NaN*X.x;
% X.cancer_cell_frac_a2=NaN*X.x;
% %X.cancer_cell_frac=NaN*X.x;
% X.ccf_hat_hack=NaN*X.x;
% X.ccf_CI95_low_hack=NaN*X.x;
% X.ccf_CI95_high_hack=NaN*X.x;
% % af binning
% daf=0.001;
% % af range
% aft=(daf/2):daf:1;

%% CS 16 Mar 2015
% W is the SCNA CCF struct 
W=load_tsv('/Users/Stewart/Projects/Cancer/UCS/ABSOLUTE/UCS.SCNA_ccf.16Mar2015.txt')
% do whatever it takes to make W.Tid match X.TID
% W.tid=regexprep(W.sample,'CLL-','')
X.SCNA_NA=NaN*X.Start_position;
X.SCNA_NB=NaN*X.Start_position;
X.SCNA_CCF_hat=NaN*X.Start_position;
X.SCNA_CCF_Low=NaN*X.Start_position;
X.SCNA_CCF_High=NaN*X.Start_position;
X.SCNA_CR=NaN*X.Start_position;
X.IS_SCNA=NaN*X.Start_position;
Tid=sort(unique(X.Tid)); N=length(Tid)
for i=1:N
    kX=find(ismember(X.Tid,Tid{i}));
    kW=find(ismember(W.Tid,Tid{i}));
    xX=xhg19(X.Chromosome(kX),X.Start_position(kX));
    w1=[];
    w1.gender=W.gender(kW);
    w1.SCNA_NA=W.NA(kW);
    w1.SCNA_NB=W.NB(kW);
    w1.SCNA_CCF_hat=W.CCF_hat(kW);
    w1.SCNA_CCF_Low=W.CCF_Low(kW);
    w1.SCNA_CCF_High=W.CCF_High(kW);
    w1.SCNA_CR=W.tau(kW);
    w1.SCNA_purity=W.purity(kW);
    w1.SCNA_ploidy=W.ploidy(kW);
    w1.IS_SCNA=W.IS_SCNA(kW);
    w1.x1=xhg19(W.Chromosome(kW),W.Start_bp(kW));
    w1.x2=xhg19(W.Chromosome(kW),W.End_bp(kW));
    [ma,mb,dab,fa,nb]=matchI(xX,w1,100,0);
    i1=find(ma>0);
    m1=ma(i1);
    X.SCNA_NA(kX(i1),1)=w1.SCNA_NA(m1);
    X.SCNA_NB(kX(i1),1)=w1.SCNA_NB(m1);
    X.SCNA_CCF_hat(kX(i1),1)=w1.SCNA_CCF_hat(m1);
    X.SCNA_CCF_Low(kX(i1),1)=w1.SCNA_CCF_Low(m1);
    X.SCNA_CCF_High(kX(i1),1)=w1.SCNA_CCF_High(m1);
    X.SCNA_CR(kX(i1),1)=w1.SCNA_CR(m1);   
    X.SCNA_purity(kX(i1),1)=w1.SCNA_purity(m1);    
    X.SCNA_ploidy(kX(i1),1)=w1.SCNA_ploidy(m1);    
    X.IS_SCNA(kX(i1),1)=w1.IS_SCNA(m1);    
end
% now X has SCNA info 
%% CS 16 Mar 2015

i=1
X1=[]
ccfA=0:0.01:1;
for i=1:N
    s1=A.sample{i}
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
    x1.modal_a1(i1)=(t1.modala1(m1));
    x1.modal_a2(i1)=(t1.modala2(m1));
    x1.cancer_cell_frac_a1(i1)=(t1.cancercellfraca1(m1));
    x1.cancer_cell_frac_a2(i1)=(t1.cancercellfraca2(m1));

    %% CS 16 Mar 2015
    % this shouldn't be needed with the SCNA CCF     
    if any(x1.q_hat==0)
        x1.q_hat(x1.q_hat==0)=1;
        x1.HS_q_hat_2(x1.q_hat==0)=1;
    end
        
    x1.CNN=2+0*x1.Start_position;
    x1.SCNA_q_hat=x1.SCNA_NA+x1.SCNA_NB;
    if any(ismember(x1.Chromosome,'X'))
        %keyboard
        k=find(ismember(x1.Chromosome,'X'));
        if ismember(g1.gender,'Male')
            x1.q_hat(k)=1;
            x1.HS_q_hat_1(k)=0;
            x1.HS_q_hat_2(k)=1;
            x1.SCNA_q_hat(k)=1;
            x1.SCNA_NA(k)=0;
            x1.SCNA_NB(k)=1;            
            x1.CNN(k)=1;
        else
            x1.q_hat(k)=2;
            x1.HS_q_hat_1(k)=1;
            x1.HS_q_hat_2(k)=1;
            x1.SCNA_q_hat(k)=2;
            x1.SCNA_NA(k)=1;
            x1.SCNA_NB(k)=1;      
        end
    end
    
    % purity
    x1.pur=A.purity(ismember(A.tid,tid1))+0*x1.Start_position;
    %x1.dna_fraction_in_tumor=(x1.pur.*x1.corrected_total_cn)./(x1.pur.*x1.corrected_total_cn+2*(1-x1.pur));
    %x1.dna_fraction_in_tumor=(x1.pur.*x1.modal_total_cn)./(x1.pur.*x1.modal_total_cn+2*(1-x1.pur));
    % x1.dna_fraction_in_tumor=(x1.pur.*x1.q_hat)./(x1.pur.*x1.q_hat+x1.CNN.*(1-x1.pur));
    %% CS 16 Mar 2015
    x1.dna_fraction_in_tumor=(x1.pur.*x1.SCNA_q_hat)./(x1.pur.*x1.SCNA_q_hat+x1.CNN.*(1-x1.pur));

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
        
        % CS 16 Mar 2015
        %qm=unique(sort([1 x1.modal_a1(j) x1.modal_a2(j)]));
        %qm=unique(sort([1 x1.HS_q_hat_1(j) x1.HS_q_hat_2(j)]));           
        %qm=qm(qm>0);
        % equal weight for each qm - 
        %wqm=(1/length(qm))+0*qm;
        qm=unique(sort([1 x1.SCNA_NA(j) x1.SCNA_NB(j)]));           
        qm=qm(qm>0);
        % equal weight for each qm - 
        wqm=(1/length(qm))+0*qm;
        % correct for SCNA ccf
        % mutation can appear on SCNA_ccf_hat or copy 1 wiht 1-SCNA_ccf_hat
        wqm(1)=1-x1.SCNA_CCF_hat(j)*wqm(1);
        wqm(2:end)=x1.SCNA_CCF_hat(j)*wqm(2:end);
        wqm/sum(wqm);
        
        %ccf=0:daf:1;
        ccf=(daf/2):daf:1;
        pccf=0*ccf;
        pccfv=[];
        m=1;
        for m=1:length(qm)
            % scale ccf down by n
            %ccfb=x1.modal_total_cn(j)*(ccf/qm(m));
            ccfb=x1.SCNA_q_hat(j)*(ccf/qm(m));
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

        y=sum(reshape(x1.pccf(j,:),10,100),1);
        xb=mean(reshape(ccf,10,100),1);
        y1 = interp1(xb,y,ccfA,'linear','extrap'); y1(y1<0)=0; y1=y1/sum(y1);
        x1.pccf_hack_100(j,:)=y1;    
        [ymax,imode]=max(y1);
        x1.ccf_mode_hack_100(j,1)=ccfA(imode);
    
        % set clonal flags
        x1.clonal_hack(j) = x1.ccf_mode_hack(j)>0.99;
        x1.subclonal_hack(j) = x1.ccf_mode_hack(j)<0.99;
        x1.clonal_hack_100(j) = x1.ccf_mode_hack_100(j)>=0.95;
        x1.subclonal_hack_100(j) = x1.ccf_mode_hack_100(j)<0.95;

    end
    if isempty(X1)
        X1=x1;
    else
        X1=mergeStruct(X1,x1);
    end
    
end

X1.clonal_hack_100 = X1.ccf_mode_hack_100>=0.95;
X1.subclonal_hack_100 = X1.ccf_mode_hack_100<0.95;

        
ff=fieldnames(X1)
k=find(cellfun(@length,regexp(ff,'^x.+','match'))>0)
ff=ff(k);

n0=length(ff)
X1.ccf_ABS=0:(1/(n0-1)):1;
X1.pccf_ABS=NaN(X1.N,n0)
for j=1:n0
    X1.pccf_ABS(:,j)=X1.(ff{j});
    vlab=sprintf('pccf_%d',round(100*X1.ccf_ABS(j)));
    X1=RenameField(X1,ff{j},vlab);    
end
for j=1:length(X1.ccf_ABS)
    vlab=sprintf('pccf_hack_%d',round(100*X1.ccf_ABS(j)))
    X1.(vlab)=X1.pccf_hack_100(:,j);
end


save(['CLL_Stilgenbauer.' TODAY '.maf.mat'],'-struct','X1')
printStruct(X1,-1,['CLL_Stilgenbauer.' TODAY '.maf'])


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

%% test
clear
cd /Users/stewart/Projects/Cancer/CLL/CLL8
M=load_tsv('/Users/stewart/Projects/Cancer/CLL/CLL8/mismatches.txt')
X=load('/Users/stewart/Projects/Cancer/CLL/CLL8/CLL_Stilgenbauer.16Feb2015a.maf.mat')
M.id=regexprep(strcat(M.sample,'@',cellstr(num2str(M.Chromosome)),':',cellstr(num2str(M.Start_position))),' ','')
X.id=regexprep(strcat(X.sample,'@',X.Chromosome,':',cellstr(num2str(X.Start_position))),' ','')
X1=trimStruct(X,find(ismember(X.id,M.id)))
%% TP53, ATM, BRAF
X1=trimStruct(X,find(ismember(X.Hugo_Symbol ,{'TP53','ATM','BRAF'})))
X1=trimStruct(X1,find(~ismember(X1.Variant_Classification ,{'Silent','Intron','IGR'})))

[i m]=ismember(X1.id,M.id)
M=trimStruct(M,m)

ff=fieldnames(M)
ff=regexp(ff,'^x\d+','match');ff=[ff{:}]; ff=ff(~ismember(ff,'x1000'));n0=length(ff)
M.ccf_ABS=0:(1/(n0-1)):1;
M.pccf_ABS=NaN(length(M.x0),n0)
for j=1:n0
    M.pccf_ABS(:,j)=M.(ff{j});
end

%%
ccf=((1:1000)-1)/1000;
X1.af=X1.t_alt_count./(X1.t_alt_count+X1.t_ref_count);
for i=1:length(M.x0)
    q=mod(i-1,4)+1;
    subplot(2,2,q);
    [M.id{i} ' ' M.Hugo_Symbol{i}]
    plot(M.ccf_ABS,M.pccf_ABS(i,:))
    title([M.id{i} ' ' M.Hugo_Symbol{i}])

    text(0.05,0.9,sprintf('pur=%.2f af=%.2f',X1.pur(i),X1.af(i)),'units','normalized')
    text(0.05,0.85,sprintf('t_alt=%d t_ref=%d',X1.t_alt_count(i),X1.t_ref_count(i)),'units','normalized')
    text(0.05,0.8,sprintf('q1=%d q2=%d',X1.HS_q_hat_1(i),X1.HS_q_hat_2(i)),'units','normalized')
    k=find(X1.pccf(i,:)>0);
    x=ccf(k)';
    y=10*X1.pccf(i,k)';
    x1 = M.ccf_ABS';
    y1 = interp1(x,y,x1,'linear','extrap'); y1=y1/sum(y1);
    line(x1,y1,'color','r')
    line(M.modal(i)*[1 1 ],[0 max([y1' M.pccf_ABS(i,:)])],'color','r')
    line(X1.ccf_mode_hack(i)*[1 1 ],[0 max([y1' X1.pccf_ABS(i,:)])],'color','m')
    line(M.ccf_hat(i)*[1 1 ],[0 max([y1' M.pccf_ABS(i,:)])],'color','b')
    xa=str2double(M.i_tumor_f(i))*2/M.purity(i);
    %line(xa*[1 1 ],[0 max([y M.pccf_ABS(i,:)])],'color',0.5*[1 1 1])
    xlim([0 1.01])
    xlabel('CCF')
    ylabel('pCCF')
    legend('ABS','HACK','location','N')
    %keyboard ; % dbcont
    if (q==4)
        %feps=['/Users/stewart/Projects/Cancer/CLL/plots/CLL_Stilgenbauer.' M.id{i} '.' M.Hugo_Symbol{i} '.' TODAY '.eps']
        feps=['/Users/stewart/Projects/Cancer/CLL/plots/CLL_Stilgenbauer.mismatched.ccfs.' num2str(round(i/4)) '.' TODAY '.eps']
        saveas(gcf,feps,'epsc')
        clf
    end
end
feps=['/Users/stewart/Projects/Cancer/CLL/plots/CLL_Stilgenbauer.mismatched.ccfs.' num2str(round(i/4)) '.' TODAY '.eps']
saveas(gcf,feps,'epsc')


%%
clear
cd /Users/stewart/Projects/Cancer/CLL/CLL8
X=load('/Users/stewart/Projects/Cancer/CLL/CLL8/CLL_Stilgenbauer.16Feb2015a.maf.mat')
X.id=regexprep(strcat(X.sample,'@',X.Chromosome,':',cellstr(num2str(X.Start_position))),' ','')
%% TP53, ATM, BRAF
X1=trimStruct(X,find(ismember(X.Hugo_Symbol ,{'TP53','ATM','BRAF'})))
X1=trimStruct(X1,find(~ismember(X1.Variant_Classification ,{'Silent','Intron','IGR'})))
[x, k]=sort(X1.Tumor_Sample_Barcode)
X1=trimStruct(X1,k)
[x, k]=sort(X1.Hugo_Symbol)
X1=trimStruct(X1,k)

tabulate(X.ccf_median_hack>=0.95)
tabulate(X.ccf_median_hack>=0.85)

tabulate(X1.ccf_median_hack>=0.95)
tabulate(X1.ccf_median_hack>=0.85)

MSG={'SF3B1', 'ATM', 'TP53', 'POT1', 'NOTCH1', 'XPO1', 'BIRC3', 'RPS15', 'BRAF', 'EGR2', 'MYD88', 'KRAS', 'MAP2K1', 'DDX3X', 'NRAS', 'CHD2', 'SAMHD1', 'FUBP1', 'FBXW7', 'DYRK1A', 'HIST1H1E', 'NXF1', 'IRF4', 'MGA', 'PTPN11', 'ZMYM3', 'FAM50A', 'IKZF3', 'MED12', 'EWSR1', 'IGLL5', 'TRAF3', 'BAZ2A', 'TTK', 'ELF4', 'GNB1', 'CARD11', 'TRAF2', 'BRCC3', 'CHEK2', 'HIST1H1B', 'ADAM30', 'XPO4', 'ASXL1', 'PIM1', 'NHS'}

X1=trimStruct(X,find(ismember(X.Hugo_Symbol ,MSG)))
X1=trimStruct(X1,find(~ismember(X1.Variant_Classification ,{'Silent','Intron','IGR'})))
tabulate(X1.ccf_median_hack>=0.95)
tabulate(X1.ccf_median_hack>=0.85)

%%
ccf=((1:1000)-1)/1000;
x=mean(reshape(ccf,10,100),1)
ccfA=X.ccf_ABS
X.af=X.t_alt_count./(X.t_alt_count+X.t_ref_count);
for i=1:length(X.x)
    y=sum(reshape(X.pccf(i,:),10,100),1);
    y1 = interp1(x,y,ccfA,'linear','extrap'); y1=y1/sum(y1);
    X.pccf_hack_100(i,:)=y1;    
    [ymax,imode]=max(y1);
    X.ccf_mode_hack_100(i,1)=X.ccf_ABS(imode);
end
P=[]
P.XLAB='mode CCF 1000 bins'
P.YLAB='mode CCF 100  bins'
hist2(X.ccf_mode_hack,X.ccf_mode_hack_100,100,100,[],P)

MSG={'SF3B1', 'ATM', 'TP53', 'POT1', 'NOTCH1', 'XPO1', 'BIRC3', 'RPS15', 'BRAF', 'EGR2', 'MYD88', 'KRAS', 'MAP2K1', 'DDX3X', 'NRAS', 'CHD2', 'SAMHD1', 'FUBP1', 'FBXW7', 'DYRK1A', 'HIST1H1E', 'NXF1', 'IRF4', 'MGA', 'PTPN11', 'ZMYM3', 'FAM50A', 'IKZF3', 'MED12', 'EWSR1', 'IGLL5', 'TRAF3', 'BAZ2A', 'TTK', 'ELF4', 'GNB1', 'CARD11', 'TRAF2', 'BRCC3', 'CHEK2', 'HIST1H1B', 'ADAM30', 'XPO4', 'ASXL1', 'PIM1', 'NHS'}
X1=trimStruct(X,find(ismember(X.Hugo_Symbol ,MSG)))
X1=trimStruct(X1,find(~ismember(X1.Variant_Classification ,{'Silent','Intron','IGR'})))

k=find((X1.ccf_mode_hack>0.99)&(X1.ccf_mode_hack_100<0.95))

X1=trimStruct(X1,k)

%%
clf
ccf=((1:1000)-1)/1000;
ccfA=X1.ccf_ABS
X1.af=X1.t_alt_count./(X1.t_alt_count+X1.t_ref_count);
for i=1:length(X1.x)
    q=mod(i-1,4)+1;
    subplot(2,2,q);
    [X1.id{i} ' ' X1.Hugo_Symbol{i}]
    plot(X1.ccf_ABS,X1.pccf_ABS(i,:),ccfA,X1.pccf_hack_100(i,:))
    title([X1.id{i} ' ' X1.Hugo_Symbol{i}])

    text(0.05,0.9,sprintf('pur=%.2f af=%.2f',X1.pur(i),X1.af(i)),'units','normalized')
    text(0.05,0.85,sprintf('t_alt=%d t_ref=%d',X1.t_alt_count(i),X1.t_ref_count(i)),'units','normalized')
    text(0.05,0.8,sprintf('q1=%d q2=%d',X1.HS_q_hat_1(i),X1.HS_q_hat_2(i)),'units','normalized')
    k=find(X1.pccf(i,:)>0);
    x=ccf(k)';
    y=10*X1.pccf(i,k)';
    y1 = interp1(x,y,ccfA,'linear','extrap'); y1=y1/sum(y1);
    line(ccfA,y1,'color','m')
    line(X1.ccf_mode_hack(i)*[1 1 ],[0 max([y1 X1.pccf_ABS(i,:)])],'color','m')
    xlim([0 1.01])
    xlabel('CCF')
    ylabel('pCCF')
    legend('ABS','HACK','location','N')

    
    keyboard ; % dbcont
    if (q==4)
        %feps=['/Users/stewart/Projects/Cancer/CLL/plots/CLL_Stilgenbauer.' M.id{i} '.' M.Hugo_Symbol{i} '.' TODAY '.eps']
        feps=['/Users/stewart/Projects/Cancer/CLL/plots/CLL_Stilgenbauer.mismatched.ccfs.diff.' X1.Hugo_Symbol{i} '.'   num2str(round(i/4)) '.' TODAY '.eps']
        saveas(gcf,feps,'epsc')
        clf
    end
    
end
feps=['/Users/stewart/Projects/Cancer/CLL/plots/CLL_Stilgenbauer.mismatched.ccfs.diff.'  X1.Hugo_Symbol{i} '.'  num2str(round(i/4)) '.' TODAY '.eps']
saveas(gcf,feps,'epsc')

hist2(X1.ccf_mode_hack,X1.ccf_mode_hack_100,100)


%% test these 
x1=X1;
% af binning
daf=0.001;
% af range
aft=(daf/2):daf:1
j=find(ismember(x1.id,'CLL-GCLL-0049-Tumor-SM-4MHYW@9:139390649'))

for j=1:length(x1.pur)
    [x1.id{j} ' ' x1.Hugo_Symbol{j}]
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
    %line(ccf,x1.paft(j,:),'linestyle','--','color',0.5*[1 1 1]);
    pccfx=1-sum(pccf);
    %pccf=pccf/sum(pccf);
    x1.pccf(j,:)=pccf;
    
    
    
    keyboard
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
