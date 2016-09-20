function [TiN,log_odds_curve,CI_h,CI_l]=TiNLOH(germline_snps,somatic_snps,firehose,cyt_file,bed,cnv_blacklist,mode,pair_id)


if isequal(mode,'legacy')

disp('Merging Germline and Tumor HetSites')
m=load_struct(germline_snps);
if isfield(m,'genotype')
m=reorder_struct(m,ismember(m.genotype,'0/1'));
m=reorder_struct(m,~ismember(m.allelic_depth,''));
end
[sums names]=count(m.Start_position,-1);
m=reorder_struct(m,~ismember(m.Start_position,names(sums>1)));
mT=load_struct(somatic_snps);
[sums names]=count(mT.Start_position,-1);
mT=reorder_struct(mT,~ismember(mT.Start_position,names(sums>1)));
mT=reorder_struct(mT,ismember(mT.Start_position,m.Start_position));
m=reorder_struct(m,ismember(m.Start_position,mT.Start_position));


call.contig=chromosome2num_legacy(m.Chromosome);
call.position=str2double(m.Start_position);
call.t_alt_count=str2double(mT.i_t_alt_count);
call.t_ref_count=str2double(mT.i_t_ref_count);
if isfield(m,'genotype')
for i=1:slength(m)
acs=split(m.allelic_depth{i},',');
call.n_ref_count(i,1)=str2double(acs{1});
call.n_alt_count(i,1)=str2double(acs{2});
end
else
call.n_ref_count=str2double(m.n_ref_count);
call.n_alt_count=str2double(m.n_alt_count);
end

elseif isequal(mode, 'call_stats')
disp('using annotated callstats for LOH TiN')

call=load_struct(germline_snps);
disp('loaded call stats file preprocessing file')
call.contig=chromosome2num_legacy(call.contig);
call.position=str2double(call.position);
call.t_alt_count=str2double(call.t_alt_count);
call.t_ref_count=str2double(call.t_ref_count);
call.n_alt_count=str2double(call.n_alt_count);
call.n_ref_count=str2double(call.n_ref_count);
call.PoN_Germline=str2double(call.PoN_Germline);
call.alt_count_greater10_af_greater_20percent=str2double(call.alt_count_greater10_af_greater_20percent);
call.normal_f=str2double(call.normal_f);
if isfield(call,'total_reads')
call.total_reads=str2double(call.total_reads);
else
    call.total_reads=call.t_alt_count+call.t_ref_count+call.n_alt_count+call.n_ref_count;
end
call.observed_in_normals_count=str2double(call.observed_in_normals_count);


af_range=[0:.01:1];
call.PoN_Artifact=str2double(call.PoN_Artifact);
cut=quantile(call.total_reads,.8);
xCoV=call.total_reads<cut;
xPoN=(call.alt_count_greater10_af_greater_20percent>.05&call.alt_count_greater10_af_greater_20percent<.6&call.observed_in_normals_count>40);
xART=call.PoN_Artifact<.1;
xAF=call.normal_f<.95&call.normal_f>.05;
call=reorder_struct(call,xPoN&xAF&xART&xCoV);

end
call=reorder_struct(call,~isnan(call.contig));
call=reorder_struct(call,((call.t_alt_count+call.t_ref_count)>20)&((call.n_alt_count+call.n_ref_count)>20));
call=reorder_struct(call,~(call.contig==23|call.contig==24));
%cnv_blacklist=load_struct('/xchip/cga/reference/gistic2/CNV.hg19.bypos.111213.txt')
call.x=xhg19(call.contig,call.position);
disp('Done preprocessing file')


if exist(bed, 'file')
disp('Using Bed File for LOH detection')

    b=load_struct(bed);
    
      ff=fieldnames(b)
    if ~ismember('Num_Probes',ff)
        % Jabba seg.txt ? 
        if ismember('width',ff)& ismember('cn',ff) &  ismember('start',ff)    &  ismember('end',ff)      
            b.Num_Probes=cellstr(num2str(1+round(str2double(b.width)/1e3))); % hack 1 probe per 1kbp
            b.Chromosome=b.chrom;   
            b.Start=b.start;   
            b.End=b.end;   
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
    call.tumor_f=str2double(call.tumor_f);
    for i=1:slength(b)
    win=reorder_struct(call,call.x>=b.xs(i)&call.x<=b.xe(i));
    llikelihoods=(binocdf(win.t_alt_count,(win.t_alt_count+win.t_ref_count),.5));
    llikelihoods_artifacts=(binocdf(win.t_alt_count,(win.t_alt_count+win.t_ref_count),.1));
    b.p_artifacts(i,1)=sum(llikelihoods_artifacts<.01)/length(llikelihoods_artifacts);
    b.p_loh(i,1)=sum(llikelihoods>.98|llikelihoods<.02)/length(llikelihoods);
    b.nMutbelow(i,1)=sum(win.tumor_f<.5);
    b.nMutabove(i,1)=sum(win.tumor_f>.5);
    
    end
    
    
    disp('Done with LOH detection')

b=reorder_struct(b,(b.p_loh>.8&~isnan(b.p_loh)&(b.nMutbelow./(b.nMutabove+b.nMutbelow))<.70));
if slength(b)==0
          disp('No LOH regions detected')
    TiN=-1;
    CI_l=-1;
    CI_h=-1;
    save(sprintf('%s.TinN.txt',pair_id),'TiN','-ascii','-double')
    
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
        if slength(cfilt)>10 && (cfilt.position(end)-cfilt.position(1))>10000
        call_filtered=mergeStruct(cfilt,call_filtered);
        call_filtered=rmfield(call_filtered,'N');
        end
    end
end
call_filtered.llikelihoods=(binocdf(call_filtered.t_alt_count,(call_filtered.t_alt_count+call_filtered.t_ref_count),.5));
call_filtered=reorder_struct(call_filtered,call_filtered.llikelihoods>.98|call_filtered.llikelihoods<.02);
end  
else
    
disp('Detecting Regions of LOH')

cytoband=load_struct(cyt_file);
cytoband.arm=cellfun(@(x) x(1),strcat(cytoband.name,'.'),'UniformOutput',0);
% unique(cellfun(@(x) x{1},split(cytoband.name,'.'),'UniformOutput',0))
armlvl.chrom{1}=cytoband.chrom{1};
armlvl.start{1}=cytoband.chromStart{1};
armlvl.arm{1}=cytoband.arm{1};
cytoband.armname=strcat(cytoband.chrom,cytoband.arm);
armlvl.armname{1}=cytoband.armname{1};
counter=1;

for i=2:slength(cytoband)
   if isequal(cytoband.armname{i},armlvl.armname{counter})
   else
       armlvl.end{counter,1}=cytoband.chromEnd{i-1};
       counter=counter+1;
       armlvl.arm{counter,1}=cytoband.arm{i};
       armlvl.chrom{counter,1}=cytoband.chrom{i};
       armlvl.start{counter,1}=cytoband.chromStart{i};
       armlvl.armname{counter,1}=cytoband.armname{i};
   end
end
armlvl.end{counter,1}=cytoband.chromEnd{i};

armlvl=reorder_struct(armlvl,~ismember(armlvl.chrom,'chrX'));
armlvl=reorder_struct(armlvl,~ismember(armlvl.chrom,'chrY'));

armlvl.chrom=chromosome2num_legacy(armlvl.chrom);
armlvl.start=str2double(armlvl.start);
armlvl.end=str2double(armlvl.end);
armlvl.xstart=xhg19(armlvl.chrom,armlvl.start);
armlvl.xend=xhg19(armlvl.chrom,armlvl.end);


call.tumor_f=call.t_alt_count./(call.t_alt_count+call.t_ref_count);
call.normal_f=call.n_alt_count./(call.n_alt_count+call.n_ref_count);
clear wins
% options = statset('Display','final','TolFun',1e-10,'MaxIter',10000);
counter=1;
for i=1:slength(armlvl)
    window=armlvl.xstart(i):1.5e7:armlvl.xend(i);
    if length(window)==1
        window(1)=armlvl.xstart(i);
        window(2)=armlvl.xend(i);
    else
    window(end)=armlvl.xend(i);
    end
    for wi=1:length(window)-1
    wins.start(counter,1)=window(wi);
    wins.end(counter,1)=window(wi+1);
    wins.armname{counter,1}=armlvl.armname{i};
    win=reorder_struct(call,call.x>wins.start(counter)&call.x<wins.end(counter));
    wins.nMut(counter,1)=slength(win);
%     if wins.nMut(counter)>10000
% %         if i==206
% %             keyboard()
% %         end
%         %g=gmdistribution.fit([win.tumor_f]*100,2,'Regularize',1,'options',options);
%         
%         gm(counter,1).mu=g.mu;
%         gm(counter,1).PComponents=g.PComponents;
%         wins.mu_lower(counter,1)=min(gm(counter,1).mu)/100;
%         wins.mu_upper(counter,1)=max(gm(counter,1).mu)/100;
%         wins.pp_lower(counter,1)=min(gm(counter,1).PComponents);
%         wins.pp_upper(counter,1)=max(gm(counter,1).PComponents);
%     else
%         wins.mu_lower(counter,1)=0;
%         wins.mu_upper(counter,1)=0;
%         wins.pp_lower(counter,1)=0;
%         wins.pp_upper(counter,1)=0;
%     end
%     
    llikelihoods=(binocdf(win.t_alt_count,(win.t_alt_count+win.t_ref_count),.5));
    llikelihoods_artifacts=(binocdf(win.t_alt_count,(win.t_alt_count+win.t_ref_count),.1));
    wins.p_artifacts(counter,1)=sum(llikelihoods_artifacts<.01)/length(llikelihoods_artifacts);
    wins.p_loh(counter,1)=sum(llikelihoods>.98|llikelihoods<.02)/length(llikelihoods);
    wins.bpd_normal(counter,1)=nanmean(binopdf(win.n_alt_count,(win.n_alt_count+win.n_ref_count),.5));
    wins.means_tumor(counter,1)=mean(win.tumor_f);
    wins.means_normal(counter,1)=mean(win.normal_f);
    wins.std_tumor(counter,1)=std(win.tumor_f);
    wins.std_normal(counter,1)=std(win.normal_f);
    wins.nMutbelow(counter,1)=sum(win.tumor_f<.5);
    wins.nMutabove(counter,1)=sum(win.tumor_f>.5);
    counter=counter+1;
    end
end
% wins=reorder_struct(wins,log(wins.bpd_tumor+1e-50)<-4);
% k_tumor=wins.means_tumor<(mean(wins.means_tumor)-(4*std(wins.means_tumor)))|wins.means_tumor>(mean(wins.means_tumor)+(4*std(wins.means_tumor)));
% k_normal=wins.means_normal<(mean(wins.means_normal)-(4*std(wins.means_normal)))|wins.means_normal>(mean(wins.means_normal)+(4*std(wins.means_normal)));
w_mean_n=nanmean(wins.means_normal);
w_mean_t=nanmean(wins.means_tumor);
w_std_n=nanstd(wins.means_normal);
w_std_t=nanstd(wins.means_tumor);

blacklist=reorder_struct(wins,wins.p_loh<.45|wins.nMut<10|wins.p_artifacts>.15|wins.means_tumor<(w_mean_t-(2*w_std_t))|wins.means_normal<(w_mean_n-(2*w_std_n))|wins.nMutbelow./(wins.nMut)>.65);
call_filtered=call;
call_filtered=reorder_struct(call_filtered,~(abs(call_filtered.tumor_f-.5)<abs(call_filtered.normal_f-.5)));

%call_filtered=reorder_struct(call_filtered,~(call_filtered.normal_f>.5&call_filtered.tumor_f<.5));
%call_filtered=reorder_struct(call_filtered,~(call_filtered.normal_f<.5&call_filtered.tumor_f>.5));
    
for i=1:slength(blacklist)
        win=reorder_struct(call,call.x>=blacklist.start(i)&call.x<=blacklist.end(i));
        call_filtered=reorder_struct(call_filtered,~ismember(call_filtered.x,win.x));
end
end
 figure()
    hold on
plot(call.x,call.t_alt_count./(call.t_alt_count+call.t_ref_count),'o','Color',[.8 .8 .8])

if slength(call_filtered)>0
   
plot(call_filtered.x,(call_filtered.t_alt_count./(call_filtered.t_alt_count+call_filtered.t_ref_count)),'b.')
plot(call_filtered.x,(call_filtered.n_alt_count./(call_filtered.n_alt_count+call_filtered.n_ref_count)),'r.')
legend('Background','Tumor LOH Hets','Normal LOH Hets')
ylabel('allele fraction','FontSize',18)

print(gcf,'-dpdf',sprintf('%s_RegionsUsed.pdf',pair_id))
% figure()
% subplot(3,1,1)
% title('Regions consitent with 10% af')
% hist(wins.p_artifacts,50)
% subplot(3,1,2)
% title('Windowed af means of hets in the tumor')
% hist(wins.means_tumor,50)
% subplot(3,1,3)
% title('Windowed af means of hets in the normal')
% hist(wins.means_normal,50)
% if isequal(firehose,'0')
% print(gcf,'-depsc','QC.eps')
% end


% af=[0:.01:1];
% call_filtered.arm=cell(slength(call_filtered),1);
% call_filtered=reorder_struct(call_filtered,~(call_filtered.contig==23|call_filtered.contig==24));
% for i=1:slength(armlvl)
%     k=call_filtered.contig==armlvl.chrom(i)&call_filtered.position<armlvl.end(i)&call_filtered.position>armlvl.start(i);
%     c=reorder_struct(call_filtered,call_filtered.contig==armlvl.chrom(i)&call_filtered.position<armlvl.end(i)&call_filtered.position>armlvl.start(i));
%     call_filtered.arm(k)=armlvl.armname(i);
%     afbeta=zeros(1,length(af));
%     for j=1:slength(c)
%         %get alt and ref for each mutation 
%         x2=betapdf(af,c.t_alt_count(j)+1,c.t_ref_count(j)+1); 
%         afbeta=afbeta+(x2/sum(x2));
%         %add them up and look for peaks (in cn consistent reigon (3p) )
%     end
%     
%     j=11:(length(afbeta)-11);
%     kp=af(find((afbeta(j-10)<=afbeta(j))&(afbeta(j+10)<afbeta(j))));
%     kp=kp(kp>.20&kp<.80);
%     armlvl.afbeta{i,1}=afbeta;
%     if length(kp)==1
%         armlvl.loh(i,1)=0;
%     elseif max(kp)>.51
%         armlvl.loh(i,1)=1;
%     end
% end
 
 
disp('Files Merged Fitting TiN on LOH')
iter=1;
[TiN(iter) tin_curve(iter,:) call_filtered]=loglikelihood_fit_LOH(call_filtered);
sg=1;
% while sg 
     for j=1:slength(call_filtered)
         call_filtered.pTiNs(j,:)=call_filtered.pTiNs(j,:)/sum(call_filtered.pTiNs(j,:));
 [v call_filtered.mode_pTin(j,:)]=max(call_filtered.pTiNs(j,:));
     end
figure()
plot(call_filtered.tumor_f,call_filtered.mode_pTin,'b.')
hold on
plot(call_filtered.normal_f,call_filtered.mode_pTin,'r.')
%     legend('Tumor','Normal')
%     xlabel('allele fraction')
%     ylabel('TiN_mode')
%    % left_tail=find(cumsum(sum(call_filtered.pTiNs,1)/sum(sum(call_filtered.pTiNs,1)))<.05,1,'last');
%     right_tail=find(cumsum(sum(call_filtered.pTiNs,1)/sum(sum(call_filtered.pTiNs,1)))>.95,1,'first');
%    % call_filtered=reorder_struct(call_filtered,call_filtered.mode_pTin>=left_tail);
%    call_filtered=reorder_struct(call_filtered,call_filtered.mode_pTin<=right_tail);
%     
%     disp(sprintf('TinN fit: %d',TiN(iter)));
%     iter=iter+1;
%     
%     [TiN(iter) log_odds_curve(iter,:) call_filtered]=loglikelihood_fit_LOH(call_filtered);
%     delta=abs(TiN(iter)-TiN(iter-1));
%     sg=(iter < 10)&(.015<delta);
% end
    
disp(sprintf('TinN fit: %d ',TiN(iter)/100));
    

figure()
hold on
plot([0:1:100],-log((tin_curve-min(tin_curve))./sum(tin_curve)))
% plot([0:0.01:1],log_odds_curve(end,:),'b--')
% plot([TiN(end)/100 TiN(end)/100],[min(log_odds_curve(end,:))-1 max(log_odds_curve(end,:))+1],'k--')
xlabel('Tumor In Normal Level','FontSize',18)
ylabel('LogLikelihood','FontSize',18)


print(gcf,'-png',sprintf('%s_LogOddsCurve.png',pair_id))

TinN=TiN(end);
log_odds_curve=-log((tin_curve-min(tin_curve))./sum(tin_curve));
CI_l=max([0 find(log_odds_curve(1:TiN+1)>log_odds_curve(TiN+1)+1,1,'last')]);
CI_h=min([100 find(log_odds_curve(TiN+1:end)>log_odds_curve(TiN+1)+1,1,'first')+TiN]);
n_log_c=exp(-log_odds_curve(end,:))/sum(exp(-log_odds_curve(end,:)));
%[CI_l,CI_h]=findCIp(n_log_c,.05,.001);
CIs=[CI_l,CI_h];
TinN=TinN/100;
CIs=CIs./100;
save(sprintf('%s.CI.txt',pair_id),'CIs','-ascii', '-double');
 save(sprintf('%s.TinN.txt',pair_id),'TiN','-ascii', '-double')
save(sprintf('%s_log_odds_curve.mat',pair_id),'log_odds_curve');

else
    disp('No LOH regions detected')
    TiN=-1;
    CI_l=-1;
    CI_h=-1;
    save(sprintf('%s.TinN.txt',pair_id),'TiN','-ascii','-double')
    
    log_odds_curve=zeros(101,1);
    save('log_odds_curve.mat','log_odds_curve');
end

if isequal(firehose,'1')
    
    quit()
end

 


end

function test %#ok<DEFNU>
TiNLOH('/Users/amaro/Downloads/CLL-GCLL-0056-Normal-SM-41QD2.maf','/Users/amaro/Downloads/CLL-GCLL-0056-Tumor-SM-41QD2.maf',0,'/Users/amaro/Documents/cytobandWitharm.txt');
end