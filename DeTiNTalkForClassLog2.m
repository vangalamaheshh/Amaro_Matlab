
call=load_struct('/xchip/cga/gdac-prod/cga/jobResults/CallSomaticMutations/BDD-CWGGV-TP-NT-SM-17M52-SM-17M53/18485730/BDD-CWGGV-TP-NT-SM-17M52-SM-17M53.call_stats.txt');
hets=load_struct('/local/cga-fh/cga/An_Breast_Exomeplus/Pair/BDD-CWGGV-TP-NT-SM-17M52-SM-17M53/jobs/hetpulldownpost/job.55118629/Tumor.cov');
b=load_struct('/local/cga-fh/cga/deTiN_Paper_analysis/Pair/BDD-137-TP-NT-SM-PPSX-SM-PPSY/jobs/capture/AllelicCapseg/job.64675426/results/BDD-137-Tumor-SM-PPSX.tsv');
hets.Start_position=str2double(hets.Start_position);
call.tmp=call.t_alt_count;
call.t_alt_count=call.n_alt_count;
call.n_alt_count=call.tmp;
call.tmp=call.t_ref_count;
call.t_ref_count=call.n_ref_count;
call.n_ref_count=call.tmp;
call.tmp=call.tumor_f;
call.tumor_f=call.normal_f;
call.normal_f=call.tmp;
call.contig=chromosome2num_legacy(call.contig);
call.position=str2double(call.position);
call.t_alt_count=str2double(call.t_alt_count);
call.t_ref_count=str2double(call.t_ref_count);
call.n_alt_count=str2double(call.n_alt_count);
call.n_ref_count=str2double(call.n_ref_count);
call.normal_f=str2double(call.normal_f);
call.tumor_f=str2double(call.tumor_f);
call.x=xhg19(call.contig,call.position);
call.total_reads=str2double(call.total_reads);

af_range=[0:.01:1];
aL=(1:23);
xL=xhg19(aL,zeros(size(aL)));
for i=1:slength(call)
    call.dist_to_telomere(i,1)=min(abs(xL-call.x(i,1)));
end
xTelo=call.dist_to_telomere>5000000;
xCoV=call.total_reads>30;
xAF=call.normal_f<.93&call.normal_f>.07&call.tumor_f<.99;
xShift=abs(.5-call.normal_f)<=abs(.5-call.tumor_f);
call=reorder_struct(call,xCoV&xAF&xTelo);
call=reorder_struct(call,ismember(call.position,hets.Start_position));

figure()
hold on
aL=num2chrom(1:23);
xL=xhg19(aL,zeros(size(aL)));
xM=round(conv(xL,[1 1]/2, 'same')); xM(end)=[];
set(gca,'xtick',xM,'xticklabel',aL,'xlim',[1 max(xL)])
for i=1:23
    line([xL(i),xL(i)],[0,1],'color',0.25*[1 1 1],'linestyle','--','LineWidth',.5)
end
t=1;
b.nMutRemovedByBinomail=zeros(slength(b),1);
bb=1;
b.overlap_len=zeros(slength(b),1);
ff=fieldnames(b)
if ~ismember('Num_Probes',ff)
    
    if ismember('width',ff)& ismember('cn',ff) &  ismember('start',ff)    &  ismember('end',ff)
        b.Num_Probes=cellstr(num2str(1+round(str2double(b.width)/1e3))); % hack 1 probe per 1kbp
        b.Chromosome=b.chrom;
        b.Start=b.start;
        b.End=b.end;
    elseif ismember('Startbp',ff)
        b.Start=b.Startbp;
        b.End=b.Endbp;
        b.Num_Probes=b.n_probes;
    else
        fprintf('require  Chromosome, Start,End, Num_Probes fields in seg file')
    end
end
b.Num_Probes=str2double(b.Num_Probes);
b=reorder_struct(b,b.Num_Probes>20);

b.Chromosome=chromosome2num_legacy(b.Chromosome);
b.Start=str2double(b.Start);
b.End=str2double(b.End);
b.xs=xhg19(b.Chromosome,b.Start);
b.xe=xhg19(b.Chromosome,b.End);
b.f=str2double(b.f);
b.n_hets=str2double(b.n_hets);
for i=1:slength(b)
    
    win=reorder_struct(call,call.x>=b.xs(i)&call.x<=b.xe(i));
    if (sum(win.tumor_f>.5)>0 && sum(win.tumor_f<.5)>0) && (sum(win.normal_f>.5)>0 && sum(win.normal_f<.5)>0)
        
        win.direction=win.normal_f>.5;
        win_median_low=median(win.normal_f(win.direction==0));
        win_median_high=median(win.normal_f(win.direction==1));
        counter=1;
        if sum(win.direction==0)>1 && sum(win.direction==1)>1
            p=1-binocdf(sum(win.direction==0),slength(win),.5);
            while p < .1 && sum(win.direction==1)>1 && sum(win.direction==0)>1
                [v l]=max(abs(win.normal_f(win.direction==0)-win_median_low));
                arr=win.position(win.direction==0);
                win=reorder_struct(win,~(win.position==arr(l)));
                counter=counter+1;
                p=1-binocdf(sum(win.direction==0),slength(win),.5);
            end
            while p > .9 && sum(win.direction==1)>1 && sum(win.direction==0)>1
                [v l]=max(abs(win.normal_f(win.direction==1)-win_median_high));
                arr=win.position(win.direction==1);
                win=reorder_struct(win,~(win.position==arr(l)));
                counter=counter+1;
                p=1-binocdf(sum(win.direction==0),slength(win),.5);
            end
            
            b.nMutRemovedByBinomail(i,1)=counter;
            call=reorder_struct(call,~((call.x>=b.xs(i)&call.x<=b.xe(i))&~ismember(call.x,win.x)));
        end
    end
    
    
    b.p_tumor_balanced(i,1)=binocdf(sum(win.t_alt_count),sum(win.t_alt_count+win.t_ref_count),.5);
    b.p_normal_balanced(i,1)=binocdf(sum(win.n_alt_count),sum(win.n_alt_count+win.n_ref_count),.5);
    
    b.nMutbelow(i,1)=sum(win.tumor_f<.5);
    b.nMutabove(i,1)=sum(win.tumor_f>.5);
    b.nMut(i,1)=slength(win);
end

b=reorder_struct(b,(b.f<.38&b.nMut>5) & b.nMutbelow>1 & b.nMutabove>1 );

if slength(b)==0
               disp('No LOH regions detected')
              pTiN=ones(101,1);

    TiN=-1;
    CI_l=-1;
    CI_h=-1;
    save(sprintf('%s.TinN_txt',pair_id),'TiN','-ascii','-double')
    
    
    TiNCurve.allele_shift_pTiN=pTiN;
    save_struct(TiNCurve,(sprintf('%s.pTiN_Curve_txt',pair_id)));
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
        if slength(cfilt)>9 && (cfilt.position(end)-cfilt.position(1))>10000
        call_filtered=mergeStruct(cfilt,call_filtered);
        call_filtered=rmfield(call_filtered,'N');
        end
    end
end

end  


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
    b.wSeg(seg,1)=(min(AFt)/.5)^2+(b.nMut(seg)/max(b.nMut));

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

    
    
    TumorSigma(seg,:)=normpdf([-5:1:5]/100,0,b.stTM(seg)+.01);
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
figure()
hold on
aL=num2chrom(1:23);
xL=xhg19(aL,zeros(size(aL)));
xM=round(conv(xL,[1 1]/2, 'same')); xM(end)=[];
set(gca,'xtick',xM,'xticklabel',aL,'xlim',[1 max(xL)])
for i=1:23
    line([xL(i),xL(i)],[0,1],'color',0.25*[1 1 1],'linestyle','--','LineWidth',.5)
end

plot(call.x,call.t_alt_count./(call.t_alt_count+call.t_ref_count),'o','Color',[.9 .9 .9])
%plot(call.x,call.n_alt_count./(call.n_alt_count+call.n_ref_count),'x','Color',[.9 .9 .9])
plot(call_filtered.x,(call_filtered.t_alt_count./(call_filtered.t_alt_count+call_filtered.t_ref_count)),'b.')
plot(call_filtered.x,(call_filtered.n_alt_count./(call_filtered.n_alt_count+call_filtered.n_ref_count)),'r.')

