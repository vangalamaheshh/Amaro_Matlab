function [TiN, pTiN, CI_h, CI_l, b]=deTiN_allele_shift(germline_snps,firehose,bed,pair_id,germline_bed,target_interval_region,noisy_sites_blacklist,Exac)

%    figure('visible','off') 
%set(0, 'DefaultFigureVisible', 'off');
overlap = @(x1,x2,y1,y2) max([min(x2,y2)-max(x1,y1),0]);


disp('using annotated callstats for LOH TiN')
call=load_struct(germline_snps); % PoN annotated call_stats file required
disp('loaded call stats file preprocessing file')


ExAC_database=load_struct(Exac);
ExAC_database.key=strcat(ExAC_database.CHROM,ExAC_database.POS);
call.key=strcat(call.contig,call.position);
ExAC_database=reorder_struct(ExAC_database,ismember(ExAC_database.key,call.key));
call=reorder_struct(call,ismember(call.key,ExAC_database.key));
ExAC_database.AF=str2double(ExAC_database.AF);
[i m]=ismember(ExAC_database.key,call.key);
call.ExAF(m(m>0),1)=ExAC_database.AF(i);

call.contig=chromosome2num_legacy(call.contig);
call.position=str2double(call.position);
call.t_alt_count=str2double(call.t_alt_count);
call.t_ref_count=str2double(call.t_ref_count);
call.n_alt_count=str2double(call.n_alt_count);
call.n_ref_count=str2double(call.n_ref_count);
call.alt_count_greater10_af_greater_20percent=str2double(call.alt_count_greater10_af_greater_20percent);
call.normal_f=str2double(call.normal_f);
call.tumor_f=str2double(call.tumor_f);
germline_bed=load_struct_noheader(germline_bed);
germline_bed.contig=chromosome2num_legacy(germline_bed.col1);
germline_bed.position=str2double(germline_bed.col2);
call.x=xhg19(call.contig,call.position);
germline.x=xhg19(germline_bed.contig,germline_bed.position);
call=reorder_struct(call,ismember(call.x,germline.x));
targets=load_struct_noheader(target_interval_region);
targets.contig=chromosome2num_legacy(targets.col1);
targets.start=str2double(targets.col2);
targets.end=str2double(targets.col3);
targets.xs=xhg19(targets.contig,targets.start);
targets.xe=xhg19(targets.contig,targets.end);
blacklist=load_struct(noisy_sites_blacklist);
blacklist.chrom=chromosome2num_legacy(blacklist.chrom);
blacklist.start=str2double(blacklist.start);
blacklist.end=str2double(blacklist.end);
blacklist.xs=xhg19(blacklist.chrom,blacklist.start);
blacklist.xe=xhg19(blacklist.chrom,blacklist.end);

blacklist_calls=ones(slength(call),1);
for i=1:slength(call)
    if sum(call.x(i)>blacklist.xs & call.x(i)<blacklist.xe)>0
        blacklist_calls(i)=0;
    end
end

call=reorder_struct(call,blacklist_calls==1);

if isfield(call,'total_reads')
call.total_reads=str2double(call.total_reads);
else
    call.total_reads=call.t_alt_count+call.t_ref_count+call.n_alt_count+call.n_ref_count;
end
call.observed_in_normals_count=str2double(call.observed_in_normals_count);


af_range=[0:.01:1];

xCoV=call.total_reads>30;
xAF=call.normal_f<.93&call.normal_f>.07;
xEXaf=call.ExAF<.4&call.ExAF>.1;
xShift=abs(.5-call.normal_f)<=abs(.5-call.tumor_f); 
% call=reorder_struct(call,xPoN&xAF&xART&xCoV);
call=reorder_struct(call,xCoV&xAF&xEXaf);
call_targeted=ones(slength(call),1);
for i=1:slength(call)
    k=isempty(targets.xs<call.x(i) & targets.xe> call.x(i));
    if k==1
        call_targeted(i)=0;
    end
end

call=reorder_struct(call,call_targeted==1);


call=reorder_struct(call,~isnan(call.contig));
call=reorder_struct(call,((call.t_alt_count+call.t_ref_count)>20)&((call.n_alt_count+call.n_ref_count)>20));
call=reorder_struct(call,~(call.contig==23|call.contig==24));
disp('Done preprocessing file')

disp('Using Bed File for LOH detection')

    b=load_struct(bed);
    
      ff=fieldnames(b)
    if ~ismember('Num_Probes',ff)
        
        if ismember('width',ff)& ismember('cn',ff) &  ismember('start',ff)    &  ismember('end',ff)      
            b.Num_Probes=cellstr(num2str(1+round(str2double(b.width)/1e3))); % hack 1 probe per 1kbp
            b.Chromosome=b.chrom;   
            b.Start=b.start;   
            b.End=b.end;   
        else
            fprintf('require  Chromosome, Start,End, Num_Probes fields in seg file')
        end
    end
    
    af_low=[0:.01:.5];
    af_high=fliplr([.5:.01:1]);
    
    b.Num_Probes=str2double(b.Num_Probes);
    b=reorder_struct(b,b.Num_Probes>20);
    b.Chromosome=chromosome2num_legacy(b.Chromosome);
    b.Start=str2double(b.Start);
    b.End=str2double(b.End);
    b.xs=xhg19(b.Chromosome,b.Start);
    b.xe=xhg19(b.Chromosome,b.End);
    t=1;
    bb=1;
    b.overlap_len=zeros(slength(b),1);
    
    while t<=slength(targets) && bb<=slength(b)
        if targets.xe(t)<b.xs(bb)
            t=t+1;
            continue
        end
        if b.xe(bb)<targets.xs(t)
            bb=bb+1;
            continue
        end
        b.overlap_len(bb)=overlap(b.xs(bb),b.xe(bb),targets.xs(t),targets.xe(t))+b.overlap_len(bb);
        t=t+1;
    end
    
    for i=1:slength(b)
     i
    win=reorder_struct(call,call.x>=b.xs(i)&call.x<=b.xe(i));
    if sum(win.tumor_f>.5)>0 && sum(win.normal_f>.5)>0
    p_up=1-poisscdf(sum(win.tumor_f>.5),sum(win.tumor_f<.5));
    p_down=1-poisscdf(sum(win.tumor_f<.5),sum(win.tumor_f>.5));
    counter=1;
    while p_down<.5
    	[v l]=min(win.tumor_f);
    	x=ones(slength(win),1);
    	x(l)=0;
    	win=reorder_struct(win,x==1);
    	counter=counter+1;
    	p_down=1-poisscdf(sum(win.tumor_f<.5),sum(win.tumor_f>.5));
    end
    p_up=1-poisscdf(sum(win.normal_f>.5),sum(win.normal_f<.5));
    p_down=1-poisscdf(sum(win.normal_f<.5),sum(win.normal_f>.5));
   if abs(.5-p_up)>abs(.5-p_down) && p_up >0 && p_down>0
        while p_up<.5 && p_up>0
            [v l]=max(win.normal_f);
            x=ones(slength(win),1);
            x(l)=0;
            win=reorder_struct(win,x==1);
            counter=counter+1;
            p_up=1-poisscdf(sum(win.normal_f>.5),sum(win.normal_f<.5));
        end
   elseif p_up >0 && p_down>0
        while p_down<.5 && p_down>0
        [v l]=min(win.normal_f);
        x=ones(slength(win),1);
        x(l)=0;
        win=reorder_struct(win,x==1);
        counter=counter+1;
        p_down=1-poisscdf(sum(win.normal_f<.5),sum(win.normal_f>.5));
        end
    end




    call=reorder_struct(call,~((call.x>=b.xs(i)&call.x<=b.xe(i))&~ismember(call.x,win.x)));
    
    
    
    b.nMutRemovedByPois(i,1)=counter;
    
    win.direction=win.normal_f>.5;
    counter=1;
    if sum(win.direction==0)>1 && sum(win.direction==1)>1
    [tt p]=ttest2(abs(.5-win.normal_f(win.direction==0)),abs(.5-win.normal_f(win.direction==1)));
    while p < .1 && slength(win)>2
        if mean(abs(.5-win.normal_f(win.direction==0))) > mean(abs(.5-win.normal_f(win.direction==1)))
        [v l]=min(win.normal_f);
        x=ones(slength(win),1);
        x(l)=0;
        win=reorder_struct(win,x==1);
        counter=counter+1;
        [tt p]=ttest2(abs(.5-win.normal_f(win.direction==0)),abs(.5-win.normal_f(win.direction==1)));
        elseif mean(abs(.5-win.normal_f(win.direction==0))) < mean(abs(.5-win.normal_f(win.direction==1)))
        [v l]=max(win.normal_f);
        x=ones(slength(win),1);
        x(l)=0;
        win=reorder_struct(win,x==1);
        counter=counter+1;
        [tt p]=ttest2(abs(.5-win.normal_f(win.direction==0)),abs(.5-win.normal_f(win.direction==1)));
        end
    end
    b.nMutRemovedbyTtest(i,1)=counter;
    call=reorder_struct(call,~((call.x>=b.xs(i)&call.x<=b.xe(i))&~ismember(call.x,win.x)));
    end
    
    counter=1;
    if sum(win.direction==0)>1 && sum(win.direction==1)>1
        p=1-binocdf(sum(win.direction==0),slength(win),.5);
        while p < .1 && sum(win.direction==1)>1 && sum(win.direction==0)>1
            [v l]=min(win.normal_f);
            x=ones(slength(win),1);
            x(l)=0;
            win=reorder_struct(win,x==1);
            counter=counter+1;
            p=1-binocdf(sum(win.direction==0),slength(win),.5);
        end
        while p > .9 && sum(win.direction==1)>1 && sum(win.direction==0)>1
            [v l]=max(win.normal_f);
            x=ones(slength(win),1);
            x(l)=0;
            win=reorder_struct(win,x==1);
            counter=counter+1;
            p=1-binocdf(sum(win.direction==0),slength(win),.5);
        end
    
    b.nMutRemovedByBinomail(i,1)=counter;
    call=reorder_struct(call,~((call.x>=b.xs(i)&call.x<=b.xe(i))&~ismember(call.x,win.x)));
    end
    
    p_low=zeros(slength(win),51);
    p_high=zeros(slength(win),51);
    for h=1:slength(win)
        p_low(h,:)=betapdf(af_low,win.t_alt_count(h)+1,win.t_ref_count(h)+1);
        p_high(h,:)=betapdf(af_high,win.t_alt_count(h)+1,win.t_ref_count(h)+1);
    end
    b_sum_low=sum(p_low);
    b_sum_high=sum(p_high);
    beta_total=(b_sum_low.*b_sum_high);
    beta_total=beta_total./sum(beta_total);

    [~, peak_af]=max(beta_total);
    int_beta_total=cumsum(beta_total);
    m_p_up=min([1-int_beta_total(peak_af) .16]);
    m_p_down=min([int_beta_total(peak_af) .16]);
    ci_low=af_low(max([find(int_beta_total<(.32-m_p_up),1,'last') 1]));
    ci_high=af_low(max([find(int_beta_total>.68+m_p_down,1,'first') peak_af]));
    
    b.sigma_af(i,1)=(ci_high-ci_low);
    b.peak_af(i,1)=af_low(peak_af);
    
    
    llikelihoods=(binocdf(win.t_alt_count,(win.t_alt_count+win.t_ref_count),.5));
    llikelihoods_artifacts=(binocdf(win.t_alt_count,(win.t_alt_count+win.t_ref_count),.1));
    b.p_artifacts(i,1)=sum(llikelihoods_artifacts<.01)/length(llikelihoods_artifacts);
    b.p_loh(i,1)=sum(llikelihoods>.98|llikelihoods<.02)/length(llikelihoods);
    b.nMutbelow(i,1)=sum(win.tumor_f<.5);
    b.nMutabove(i,1)=sum(win.tumor_f>.5);
    b.nMoreExtreme(i,1)=sum(abs(win.tumor_f-.5)<abs(win.normal_f-.5));
    b.p_down(i,1)=1-poisscdf(b.nMutbelow(i),b.nMutabove(i));
    
    else
        b.sigma_af(i,1)=.5;
        b.peak_af(i,1)=.5;
        b.p_artifacts(i,1)=1;
        b.nMutbelow(i,1)=sum(win.tumor_f<.5);
        b.nMutabove(i,1)=sum(win.tumor_f>.5);
        b.nMoreExtreme(i,1)=sum(abs(win.tumor_f-.5)<abs(win.normal_f-.5));
        b.p_down(i,1)=1-poisscdf(b.nMutbelow(i),b.nMutabove(i));
        b.nMutRemovedByPois(i,1)=0;

    end
    
    end
    
disp('Done with LOH detection')
b.length=b.End-b.Start;
b.per_removed=b.nMutRemovedByPois./(b.nMutabove+b.nMutbelow+b.nMutRemovedByPois);
b.density=(b.nMutabove+b.nMutbelow+b.nMutRemovedByPois)./b.overlap_len;
b=reorder_struct(b,(b.peak_af+b.sigma_af<.5) & b.nMutbelow>1 & b.nMutabove>1 );
b.nMut=b.nMutbelow+b.nMutabove;

if slength(b)==0
          disp('No LOH regions detected')
              pTiN=ones(101,1);

    TiN=-1;
    CI_l=-1;
    CI_h=-1;
    save(sprintf('%s.TinN_txt',pair_id),'TiN','-ascii','-double')
    
    log_odds_curve=zeros(101,1);
    save('log_odds_curve.mat','log_odds_curve');
    if isequal(firehose,'1')
    
    quit()
    end
    call_filtered.empty=[];
else



    
    call_filtered=reorder_struct(call,call.x>=b.xs(1)&call.x<=b.xe(1));
    if slength(call_filtered)<10
        call_filtered=reorder_struct(call_filtered,zeros(slength(call_filtered),1)==1);
    end
if slength(b)>1
    for i=2:slength(b)
        cfilt=reorder_struct(call,call.x>=b.xs(i)&call.x<=b.xe(i));
        if slength(cfilt)>5 && (cfilt.position(end)-cfilt.position(1))>10000
        call_filtered=mergeStruct(cfilt,call_filtered);
        call_filtered=rmfield(call_filtered,'N');
        end
    end
end

end  

if  slength(call_filtered)>0
    


lohd=unique(call_filtered.contig);
disp('Using the following regions to estimate TiN:')
strcat(num2str(lohd))
TiN=0:.01:1;
 
call_filtered.direction=call_filtered.tumor_f>.5;
%Segment by segment
pTiN=ones(1,101);
b.pTiN_sum=ones(slength(b),101);
b.wSeg=zeros(slength(b),1);
b.pTiN_conv=ones(slength(b),101);
b.modal_TiN=NaN(slength(b),1);
for seg=1:slength(b)
    cfilt=reorder_struct(call_filtered,call_filtered.x>=b.xs(seg)&call_filtered.x<=b.xe(seg));
    if slength(cfilt)>1
    b.isused(seg,1)=1;
    AFt(1)=median(cfilt.tumor_f(cfilt.direction==1))-.5;
     AFt(2)=.5-median(cfilt.tumor_f(cfilt.direction==0));
%     
    % CHANGE THIS LINE TO BE STDEV ON MEDIAN*2 =
    % 2*1.253*sigma/np.sqrt(n)
    b.stT(seg,1)=std(abs(.5-cfilt.tumor_f));
    b.stTM(seg,1)=1.253*b.stT(seg,1)./sqrt(slength(cfilt));
    b.wSeg(seg,1)=(min(AFt)/.5)^2;

    for i=1:slength(cfilt)
       if cfilt.direction(i)==1
          pTiN(i,:)=betapdf(.5+(AFt(1)*TiN),cfilt.n_alt_count(i)+1,cfilt.n_ref_count(i)+1);
       else
          pTiN(i,:)=betapdf(.5-(AFt(2)*TiN),cfilt.n_alt_count(i)+1,cfilt.n_ref_count(i)+1);
       end
        
        pTiN(i,:)=pTiN(i,:)./sum(pTiN(i,:));
              
        
    end
    b.pTiN_down(seg,:)=zeros(1,101);
    b.pTiN_up(seg,:)=zeros(1,101);
    b.pTiN_down(seg,:)=sum(log(pTiN(cfilt.direction==0,:)));
    b.pTiN_up(seg,:)=sum(log(pTiN(cfilt.direction==1,:)));
    
    b.pTiN_down(seg,:)=b.pTiN_down(seg,:)+(1-max(b.pTiN_down(seg,:)));
    b.pTiN_up(seg,:)=b.pTiN_up(seg,:)+(1-max(b.pTiN_up(seg,:)));
    
    b.pTiN_up(seg,:)=exp(b.pTiN_up(seg,:));
    b.pTiN_down(seg,:)=exp(b.pTiN_down(seg,:));
    
    b.pTiN_up(seg,:)=b.pTiN_up(seg,:)./sum(b.pTiN_up(seg,:));
    b.pTiN_down(seg,:)=b.pTiN_down(seg,:)./sum(b.pTiN_down(seg,:));
    
    [l peak_af_down]=max(b.pTiN_down(seg,:));
    [l peak_af_up]=max(b.pTiN_up(seg,:));
    [v l]=max(pTiN');
    disp(sprintf('Fitting Segment %d',seg));


    b.pTiN_sum(seg,:)=sum(log(pTiN),1);
    
    pTiN_joint=b.pTiN_sum(seg,:)+(1-max(b.pTiN_sum(seg,:)));
    pTiN_joint=exp(pTiN_joint);

    
    
    TumorSigma(seg,:)=normpdf([-5:1:5]/100,0,b.stTM(seg));
    TumorSigma(seg,:)=TumorSigma(seg,:)/sum(TumorSigma(seg,:));
    b.pTiN_conv(seg,:)=conv(pTiN_joint,TumorSigma(seg,:),'same');
    [v l]=max(b.pTiN_conv(seg,:));
    b.modal_TiN(seg)=TiN(l);
    b.pTiN_conv(seg,:)=b.pTiN_conv(seg,:)./sum(b.pTiN_conv(seg,:));
    b.TiNCIL(seg,:)=TiN(max([find(cumsum(b.pTiN_conv(seg,:))<0.025,1,'last') 1 ]));
    b.TiNCIH(seg,:)=TiN(min([find(cumsum(b.pTiN_conv(seg,:))>0.975,1,'first') 101 ]));
    else
         b.isused(seg,1)=0;
    end
    
    pTiN=ones(1,101);
end
b=reorder_struct(b,~isnan(b.modal_TiN));

b=reorder_struct(b,~isnan(b.wSeg));
for i=1:slength(b)
[v l]=max(b.pTiN_down(i,:)); 
b.p_consistent_up_down(i,1)=b.pTiN_up(i,l);
end
b=reorder_struct(b,b.p_consistent_up_down>.005)

b.wSeg=b.wSeg./max(b.wSeg);
pTiN=sum((b.wSeg.^2)'*b.pTiN_conv,1)/sum(b.wSeg.^2);

pTiN=pTiN-min(pTiN);
pTiN=(pTiN./sum(pTiN));
[v l]=max(pTiN);



if sum(pTiN(1:l))>.16
CI_l=TiN(max([find(cumsum(pTiN)<0.16,1,'last') 1 ]));
else
    CI_l=TiN(1);
end

if sum(pTiN(l:101))>.16
CI_h=TiN(min([find(cumsum(pTiN)>0.84,1,'first') 101 ]));
else
    CI_h=TiN(101);
end




for i=1:slength(b)
    for j=1:slength(b)
        overlap_matrix(i,j)=min([sum(b.pTiN_conv(i,b.pTiN_conv(i,:)>1e-2&b.pTiN_conv(j,:)>1e-2));sum(b.pTiN_conv(j,b.pTiN_conv(i,:)>1e-2&b.pTiN_conv(j,:)>1e-2))]);
    end
end

cluster_ix=1;
index=1;
cluster_matrix=ones(size(overlap_matrix));
b.clust_assignment=zeros(slength(b),1);
[v l]=min(overlap_matrix);
[i x]=min(v);
if overlap_matrix(x,l(x))<.15
    b.clust_assignment(x,1)=1;
    b.clust_assignment(l(x),1)=2;
    for i=1:slength(b)
        if b.clust_assignment(i)==0
            if overlap_matrix(i,x)>overlap_matrix(i,l(x)) && overlap_matrix(i,x)>=.15
                b.clust_assignment(i)=1;
            elseif overlap_matrix(i,x)<overlap_matrix(i,l(x)) && overlap_matrix(i,l(x))>=.15
                b.clust_assignment(i)=2;
            end
        end
    end
else
    b.clust_assignment=ones(slength(b),1);
end

while sum(b.clust_assignment==0)>0
    x=find(b.clust_assignment==0,1);
    b.clust_assignment(x,1)=max(b.clust_assignment)+1;
    for i=1:slength(b)
        if b.clust_assignment(i)==0 && overlap_matrix(i,x)>=.15
            b.clust_assignment(i,1)=b.clust_assignment(x);
        end
    end
end




if max(b.clust_assignment)>1
   
    disp('Outliers detected')
    number_of_clusters=max(b.clust_assignment);
    figure()
    plot(pTiN,'r--')
    hold on
    colors=jet(number_of_clusters);
    legend_cell{1}='Unclustered';
    for i=1:max(b.clust_assignment)
    outliers=reorder_struct(b,b.clust_assignment==i);
    pTiN_outliers=sum((outliers.wSeg.^2)'*outliers.pTiN_conv,1)/sum(outliers.wSeg.^2);
    
    pTiN_outliers=pTiN_outliers-min(pTiN_outliers);
    pTiN_outliers=(pTiN_outliers./sum(pTiN_outliers));
    [v l_out]=max(pTiN_outliers);
    plot(pTiN_outliers,'--','Color',colors(i,:))
    legend_cell{i+1}=sprintf('Cluster %d',i);
    end
    ylabel('Likelihood','FontSize',25)
    xlabel('Tumor in Normal Percent','FontSize',25)
    xlim([0 100])
    
    legend(legend_cell)
    b=rmfield(b,{'pTiN_down','pTiN_up','pTiN_sum','pTiN_conv'});
    save_struct(b,[pair_id 'deTiN_Segments.txt']);
    f=sprintf('%s_Likelihoodcurve',pair_id);
    print(gcf,'-dpng','-r400',[f '.png'])
    figure()
    imagesc(overlap_matrix)
    f=sprintf('%s_ClusteringMatrix',pair_id);
    print(gcf,'-dpng','-r400',[f '.png'])
else
figure()
plot(pTiN)
ylabel('Likelihood','FontSize',25)
xlabel('Tumor in Normal Percent','FontSize',25)
xlim([0 100])
f=sprintf('%s_Likelihoodcurve.png',pair_id);
print(gcf,'-dpng','-r400',[f '.png'])
b=rmfield(b,{'pTiN_down','pTiN_up','pTiN_sum','pTiN_conv'});
    save_struct(b,[pair_id 'deTiN_Segments.txt']);
end

[v l]=max(pTiN);
TiN=TiN(l);
CIs=[CI_l,CI_h];
TinN=TiN;
save(sprintf('%s.CItxt',pair_id),'CIs','-ascii', '-double');
save(sprintf('%s.TinN_txt',pair_id),'TinN','-ascii', '-double')
TiNCurve.allele_shift_pTiN=pTiN';
save_struct(TiNCurve,(sprintf('%s.pTiN_Curve_txt',pair_id)));




% Plotting used regions for estimate
figure()
hold on
if slength(call)<50000
plot(call.x,call.t_alt_count./(call.t_alt_count+call.t_ref_count),'o','Color',[.9 .9 .9])
%plot(call.x,call.n_alt_count./(call.n_alt_count+call.n_ref_count),'x','Color',[.9 .9 .9])
plot(call_filtered.x,(call_filtered.t_alt_count./(call_filtered.t_alt_count+call_filtered.t_ref_count)),'b.')
plot(call_filtered.x,(call_filtered.n_alt_count./(call_filtered.n_alt_count+call_filtered.n_ref_count)),'r.')
aL=(1:23);
xL=xhg19(aL,zeros(size(aL)));
xM=round(conv(xL,[1 1]/2, 'same')); xM(end)=[];
set(gca,'xtick',xM,'xticklabel',num2chrom(aL),'xlim',[1 max(xL)])
line([xL xL]',[0+0*xL 1+0*xL]','color',0.25*[1 1 1],'linestyle','--','LineWidth',.5)
else
    plot(call.x(1:20:end),call.t_alt_count(1:20:end)./(call.t_alt_count(1:20:end)+call.t_ref_count(1:20:end)),'o','Color',[.9 .9 .9])
%plot(call.x,call.n_alt_count./(call.n_alt_count+call.n_ref_count),'x','Color',[.9 .9 .9])
plot(call_filtered.x,(call_filtered.t_alt_count./(call_filtered.t_alt_count+call_filtered.t_ref_count)),'b.')
plot(call_filtered.x,(call_filtered.n_alt_count./(call_filtered.n_alt_count+call_filtered.n_ref_count)),'r.')
aL=num2chrom(1:23);
xL=xhg19(aL,zeros(size(aL)));
xM=round(conv(xL,[1 1]/2, 'same')); xM(end)=[];
set(gca,'xtick',xM,'xticklabel',aL,'xlim',[1 max(xL)])
line([xL xL]',[0+0*xL 1+0*xL]','color',0.25*[1 1 1],'linestyle','--','LineWidth',.5)
end
ylabel('allele fraction','FontSize',18)
h=gcf;
set(h,'PaperPositionMode','auto');         


f=sprintf('%s_RegionsUsed',pair_id);
legend('Background Tumor','Tumor LOH Hets','Normal LOH Hets')
print(gcf,'-depsc ',sprintf('%s_RegionsUsed.eps',pair_id));
print(gcf,'-dpng','-r400',[f '.png'])




else
    TiN=-1;
    pTiN=ones(101,1);
    CI_h=-1;
    CI_l=-1;
    TinN=-1;
    CIs=-1;
save(sprintf('%s.CItxt',pair_id),'CIs','-ascii', '-double');
save(sprintf('%s.TinN_txt',pair_id),'TinN','-ascii', '-double');
end
if isequal(firehose,'1')
    
    quit()
end
end
