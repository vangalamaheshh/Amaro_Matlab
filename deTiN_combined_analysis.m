function [call, TinN, Iter_s]=deTiN_combined_analysis(call_stats_file,pair_id,firehose,iter_threshold,f_c,max_iter,chip,force_call_bed,mode,min_points,inital_tin,bed)
% Changing cell fields to numbers
if isstr(chip)
    chip=str2double(chip); % This is the window for TiN calculation
end

if isequal(firehose,'1') % lots of specifics to running in firehose
    iter_threshold=str2double(iter_threshold);
    f_c=str2double(f_c);
    max_iter=str2double(max_iter);
    set_verbose_level(30);
    figure('visible','off')
    min_points=str2double(min_points);
    inital_tin=str2double(inital_tin);
end

loglikeli=1;

%Checking call stats input
if ~isstruct(call_stats_file)
    
    warning('off','all');
    disp('Loading call stats file');
    
    
    call=load_struct(call_stats_file);
    call.t_ref_count=str2double(call.t_ref_count);
    call.t_alt_count=str2double(call.t_alt_count);
    call.n_ref_count=str2double(call.n_ref_count);
    call.n_alt_count=str2double(call.n_alt_count);
    call.tumor_f=call.t_alt_count./(call.t_alt_count+call.t_ref_count);
    call.normal_f=call.n_alt_count./(call.n_alt_count+call.n_ref_count);
    call.contig=chromosome2num_legacy(call.contig);
    call.position=str2double(call.position);
    call.PoN_Germline=str2double(call.PoN_Germline);
    call.alt_count_greater10_af_greater_20percent=str2double(call.alt_count_greater10_af_greater_20percent);
    call.PoN_Artifact=str2double(call.PoN_Artifact);

    if isfield(call,'total_reads')
        call.total_reads=str2double(call.total_reads);
    else
        call.total_reads=call.t_alt_count+call.t_ref_count+call.n_alt_count+call.n_ref_count;
    end
    call.observed_in_normals_count=str2double(call.observed_in_normals_count);


    af_range=[0:.01:1];


    
    
    
    
else
    call=call_stats_file;
    af_range=[0:.01:1];

end
% store full call_stats file prior to filtering
all_calls=call;
call_pre=sum(ismember(call.judgement,'KEEP'));

%limiting to 16x plus sites
call=reorder_struct(call,((call.t_alt_count+call.t_ref_count)>15)&((call.n_alt_count+call.n_ref_count)>15));

%LOH specific cut offs
xCoV=call.total_reads>30;
xPoN=(call.alt_count_greater10_af_greater_20percent>.05&call.alt_count_greater10_af_greater_20percent<.5&call.observed_in_normals_count>30);
xART=call.PoN_Artifact<.1;
xAF=call.normal_f<.93&call.normal_f>.07;
allele_shift_call=reorder_struct(call,xPoN&xAF&xART&xCoV);
cutd=min([quantile(allele_shift_call.normal_f,.15) .15]);
cutu=max([quantile(allele_shift_call.normal_f,.85) .85]);
if (sum(allele_shift_call.normal_f<cutd)/slength(allele_shift_call))>.05
  for i=cutd:-.01:0
      if (sum(allele_shift_call.normal_f<i)/slength(allele_shift_call))<=.08
          cutd=i;
          break
      end
  end
end
if (sum(allele_shift_call.normal_f>cutu)/slength(allele_shift_call))>.05
    for i=cutu:.01:1
        if (sum(allele_shift_call.normal_f>i)/slength(allele_shift_call))<=.08
          cutu=i;
          break
        end
    end
end
allele_shift_call=reorder_struct(allele_shift_call,allele_shift_call.normal_f<cutu&allele_shift_call.normal_f>cutd);

allele_shift_call=reorder_struct(allele_shift_call,~isnan(allele_shift_call.contig));
allele_shift_call=reorder_struct(allele_shift_call,((allele_shift_call.t_alt_count+allele_shift_call.t_ref_count)>20)&((allele_shift_call.n_alt_count+allele_shift_call.n_ref_count)>20));
allele_shift_call=reorder_struct(allele_shift_call,~(allele_shift_call.contig==23|allele_shift_call.contig==24));
%cnv_blacklist=load_struct('/xchip/cga/reference/gistic2/CNV.hg19.bypos.111213.txt')
allele_shift_call.x=xhg19(allele_shift_call.contig,allele_shift_call.position);
disp('Done preprocessing file')

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
    for i=1:slength(b)
    win=reorder_struct(allele_shift_call,allele_shift_call.x>=b.xs(i)&allele_shift_call.x<=b.xe(i));
    llikelihoods=(binocdf(win.t_alt_count,(win.t_alt_count+win.t_ref_count),.5));
    llikelihoods_artifacts=(binocdf(win.t_alt_count,(win.t_alt_count+win.t_ref_count),.1));
    b.p_artifacts(i,1)=sum(llikelihoods_artifacts<.01)/length(llikelihoods_artifacts);
    b.p_loh(i,1)=sum(llikelihoods>.98|llikelihoods<.02)/length(llikelihoods);
    b.nMutbelow(i,1)=sum(win.tumor_f<.5);
    b.nMutabove(i,1)=sum(win.tumor_f>.5);
    b.nMoreExtreme(i,1)=sum(abs(win.tumor_f-.5)<abs(win.normal_f-.5));
    end
    
    
    disp('Done with LOH detection')

b=reorder_struct(b,(b.p_loh>=.75&~isnan(b.p_loh)&(b.nMutbelow./(b.nMutabove+b.nMutbelow))<.70)&b.nMutbelow>1&(b.nMoreExtreme./(b.nMutbelow+b.nMutabove)<.5));
b.length=b.End-b.Start;
b.nMut=b.nMutbelow+b.nMutabove;
% b.density=b.nMut./b.length;
% cut=quantile(b.density,.95);
% b=reorder_struct(b,(b.nMut./b.length)<cut);

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
    allele_shift_call_filtered.empty=[];
else



    
    allele_shift_call_filtered=reorder_struct(allele_shift_call,allele_shift_call.x>=b.xs(1)&allele_shift_call.x<=b.xe(1));
    if slength(allele_shift_call_filtered)<10
        allele_shift_call_filtered=reorder_struct(allele_shift_call_filtered,zeros(slength(allele_shift_call_filtered),1)==1);
    end
if slength(b)>1
    for i=2:slength(b)
        cfilt=reorder_struct(allele_shift_call,allele_shift_call.x>=b.xs(i)&allele_shift_call.x<=b.xe(i));
        if slength(cfilt)>10 && (cfilt.position(end)-cfilt.position(1))>10000
        allele_shift_call_filtered=mergeStruct(cfilt,allele_shift_call_filtered);
        allele_shift_call_filtered=rmfield(allele_shift_call_filtered,'N');
        end
    end
end
allele_shift_call_filtered.llikelihoods=(binocdf(allele_shift_call_filtered.t_alt_count,(allele_shift_call_filtered.t_alt_count+allele_shift_call_filtered.t_ref_count),.5));
allele_shift_call_filtered=reorder_struct(allele_shift_call_filtered,allele_shift_call_filtered.llikelihoods>.98|allele_shift_call_filtered.llikelihoods<.02);
allele_shift_call_filtered.n_shift=abs(allele_shift_call_filtered.normal_f-.5);
allele_shift_call_filtered.t_shift=abs(allele_shift_call_filtered.tumor_f-.5);
mnshift=median(allele_shift_call_filtered.n_shift);
stdnshift=std(allele_shift_call_filtered.n_shift);
allele_shift_call_filtered=reorder_struct(allele_shift_call_filtered,allele_shift_call_filtered.n_shift<mnshift+2*stdnshift);
allele_shift_call_filtered=reorder_struct(allele_shift_call_filtered,allele_shift_call_filtered.n_shift>mnshift-2*stdnshift);

end  
    
if  slength(allele_shift_call_filtered)>0
    
% Plotting used regions for estimate
figure()
hold on
plot(allele_shift_call.x,allele_shift_call.t_alt_count./(allele_shift_call.t_alt_count+allele_shift_call.t_ref_count),'o','Color',[.9 .9 .9])
plot(allele_shift_call.x,allele_shift_call.n_alt_count./(allele_shift_call.n_alt_count+allele_shift_call.n_ref_count),'x','Color',[.9 .9 .9])
plot(allele_shift_call_filtered.x,(allele_shift_call_filtered.t_alt_count./(allele_shift_call_filtered.t_alt_count+allele_shift_call_filtered.t_ref_count)),'b.')
plot(allele_shift_call_filtered.x,(allele_shift_call_filtered.n_alt_count./(allele_shift_call_filtered.n_alt_count+allele_shift_call_filtered.n_ref_count)),'r.')

ylabel('allele fraction','FontSize',18)
h=gcf;
set(h,'PaperPositionMode','auto');         


lohd=unique(allele_shift_call_filtered.contig);
disp('Using the following regions to estimate TiN:')
strcat(num2str(lohd))
TiN=0:.01:1;
 
allele_shift_call_filtered.direction=allele_shift_call_filtered.tumor_f>.5;
%Segment by segment
pTiN=ones(slength(b), 101);
b.pTiN_sum=ones(slength(b),101);
b.wSeg=zeros(slength(b),1);
b.pTiN_conv=ones(slength(b),101);
b.modal_TiN=NaN(slength(b),1);
for seg=1:slength(b)
    cfilt=reorder_struct(allele_shift_call_filtered,allele_shift_call_filtered.x>=b.xs(seg)&allele_shift_call_filtered.x<=b.xe(seg));
    if slength(cfilt)>1
    b.isused(seg,1)=1;
    nshiftm=median(abs(.5-cfilt.normal_f));
    nstdm=std(abs(.5-cfilt.normal_f));


    cfilt=reorder_struct(cfilt,~(abs(.5-cfilt.normal_f)>=nshiftm+2*nstdm));
    cfilt=reorder_struct(cfilt,~(abs(.5-cfilt.normal_f)<=nshiftm-2*nstdm));
    nshiftm=median(abs(.5-cfilt.tumor_f));
    nstdm=std(abs(.5-cfilt.tumor_f));
    
    cfilt=reorder_struct(cfilt,~(abs(.5-cfilt.tumor_f)>=nshiftm+2*nstdm));
    cfilt=reorder_struct(cfilt,~(abs(.5-cfilt.tumor_f)<=nshiftm-2*nstdm));
    AFt(1)=median(cfilt.tumor_f(cfilt.direction==1))-.5;
    AFt(2)=.5-median(cfilt.tumor_f(cfilt.direction==0));
    
    % CHANGE THIS LINE TO BE STDEV ON MEDIAN*2 =
    % 2*1.253*sigma/np.sqrt(n)
    b.stT(seg,1)=std(abs(.5-cfilt.tumor_f));
    b.stTM(seg,1)=1.253*b.stT(seg,1)./sqrt(slength(cfilt));
    b.wSeg(seg,1)=mean(AFt)/.5;
 %   errorbar(median(cfilt.x),median(cfilt.tumor_f(cfilt.direction==1)),std(cfilt.tumor_f(cfilt.direction==1)),'.','Color',[255/255,127/255,0],'MarkerSize',20)
 %   errorbar(median(cfilt.x),median(cfilt.tumor_f(cfilt.direction==0)),std(cfilt.tumor_f(cfilt.direction==0)),'.','Color',[51/255,160/255,44/255],'MarkerSize',20)
 %   errorbar(median(cfilt.x),median(cfilt.normal_f(cfilt.direction==1)),std(cfilt.normal_f(cfilt.direction==1)),'.','Color',[255/255,127/255,0],'MarkerSize',20)
  %  errorbar(median(cfilt.x),median(cfilt.normal_f(cfilt.direction==0)),std(cfilt.normal_f(cfilt.direction==0)),'.','Color',[51/255,160/255,44/255],'MarkerSize',20)



    for i=1:slength(cfilt)
       if cfilt.direction(i)==1
          pTiN(i,:)=betapdf(.5+(AFt(1)*TiN),cfilt.n_alt_count(i)+1,cfilt.n_ref_count(i)+1);
       else
          pTiN(i,:)=betapdf(.5-(AFt(2)*TiN),cfilt.n_alt_count(i)+1,cfilt.n_ref_count(i)+1);
       end
     
        pTiN(i,:)=pTiN(i,:)./sum(pTiN(i,:));
              
        
    end
    b.pTiN_sum(seg,:)=sum(log(pTiN),1);
    pTiN=b.pTiN_sum(seg,:)+(1-max(b.pTiN_sum(seg,:)));
    pTiN=exp(pTiN);
    TumorSigma(seg,:)=normpdf([-5:1:5]/100,0,b.stTM(seg));
    TumorSigma(seg,:)=TumorSigma(seg,:)/sum(TumorSigma(seg,:));
    b.pTiN_conv(seg,:)=conv(pTiN,TumorSigma(seg,:),'same');
    [v l]=max(b.pTiN_conv(seg,:));
    b.modal_TiN(seg)=TiN(l);
   % b.pTiN_sum_conv(seg,:)=conv(b.pTiN_sum(seg,:),TumorSigma,'same');
    else
         b.isused(seg,1)=0;
    end
    

end
b=reorder_struct(b,~isnan(b.modal_TiN));
b=reorder_struct(b,b.modal_TiN<=(median(b.modal_TiN))+2*std(b.modal_TiN));
b=reorder_struct(b,b.modal_TiN>=(median(b.modal_TiN))-2*std(b.modal_TiN));
b=reorder_struct(b,~isnan(b.wSeg));

b.wSeg=b.wSeg./max(b.wSeg);
pTiN=sum((b.wSeg.^2)'*b.pTiN_conv,1)/sum(b.wSeg.^2);

pTiN=pTiN-min(pTiN);
pTiN=(pTiN./sum(pTiN));
[v l]=max(pTiN);
f=sprintf('%s_RegionsUsed',pair_id);
legend('Background Tumor','Background Normal','Tumor LOH Hets','Normal LOH Hets','Outliers')
print(gcf,'-depsc ',sprintf('%s_RegionsUsed.eps',pair_id));
saveas(gcf,[f '.png'],'png')
CI_l=TiN(max([find(cumsum(pTiN)<0.16,1,'last') 1 ]));
CI_h=TiN(min([find(cumsum(pTiN)>0.84,1,'first') 101 ]));

TiN=TiN(l);
CIs=[CI_l,CI_h];
TinN=TiN;
save(sprintf('%s.CItxt',pair_id),'CIs','-ascii', '-double');
save(sprintf('%s.TinN_txt',pair_id),'TinN','-ascii', '-double')
TiNCurve.allele_shift_pTiN=pTiN';
save_struct(TiNCurve,(sprintf('%s.pTiN_Curve_txt',pair_id)));
% full file at once

%for TiN_i=1:length(TiN) 
%p(TiN_i)=betapdf(allele_shift_call_filtered.tumor_f*TiN(TiN_i)

%end
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

%plotting

%fitting germline sites to a line
g = fittype('a*x');

%checking that samples are related
k=(call.normal_f>.4&call.tumor_f>.4);
if sum(k)/length(k)<.2
    disp('samples likely not related!')
    linfit.a=-1;
else
    
    linfit=fit(call.tumor_f(k),call.normal_f(k),g);
end

if linfit.a >1.2 || linfit.a <.8
    
    disp('samples flipped? Not related?')
    linfit.a=-1;
end


%Drawing somatic triangle and only including points which are rejected for
%normal alleles (i.e. not poor mapping etc.)

sites_to_consider=(ismember(call.judgement,'KEEP')|...
    ismember(call.failure_reasons,'germline_risk,normal_lod,alt_allele_in_normal')| ismember(call.failure_reasons,'normal_lod,alt_allele_in_normal') |...
    ismember(call.failure_reasons,'normal_lod,alt_allele_in_normal') | ismember(call.failure_reasons,'normal_lod') | ismember(call.failure_reasons,'alt_allele_in_normal') |...
    ismember(call.failure_reasons,'germline_risk')&~ismember(call.failure_reasons,',PoN') );


% making figure and filtering callstats file to candidate somatic sites
f=figure();
hold on
plot(call.tumor_f(1:20:slength(call)),call.normal_f(1:20:slength(call)),'.','color',[190/255 190/255 190/255],'marker','o');

%filters
call=reorder_struct(call,sites_to_consider);
%This region is designed to exclude sites that might be 100% af in the tumor
%and 50% in the normal (since these might include sites that are under LOH)
% Further we want to exclude very low allele fraction sites <.15% purity
% due to sequencer errors...



if isequal(mode,'deep_seq')
    k=(call.normal_f<(.4*call.tumor_f)-.05)&call.tumor_f>.15&call.normal_f<.4;
else
    k=(call.normal_f<(.4*call.tumor_f)-.05)&call.tumor_f>.15&call.tumor_f<.8|(call.tumor_f>.8&call.normal_f<.1);
end


% new point criteria... 

 for f=1:slength(call)
            %Only consider sites that have a greater tumorf than normalf
            %saves time.
            if binocdf(call.n_alt_count,60,.5)<.02
                call.logOdds(f,1)=-Inf;
                call.prob(f,1)=0;
            else
                [call.logOdds(f,1) call.prob(f,1)]=normal_contamination_estimate_fitting(call.t_alt_count(f),call.t_ref_count(f)...
                    ,call.n_alt_count(f),call.n_ref_count(f),.1,linfit.a,.001)  ;
            end
 end

k=call.logOdds>iter_threshold&call.tumor_f>.20;

refline(linfit.a,0);
ylim([0,1]);
xlim([0,1])
plot(call.tumor_f,call.normal_f,'bo')
plot(call.tumor_f(ismember(call.judgement,'KEEP')),call.normal_f(ismember(call.judgement,'KEEP')),'.','color',[50/255 190/255 50/255],'marker','o','LineWidth',2);




%only run on tumors that have more than 10 sites in the somatic
%triangle and that pass the germline presence filter
if sum(k)>=min_points && linfit.a~=-1
    disp('fitting possible tumor in normal')
    plot(call.tumor_f(k),call.normal_f(k),'r.','MarkerSize',10)
    
    if loglikeli% This if statement is just here from old code...
        
        
        call.dont_include=zeros(slength(call),1);
        iter=1;
        [Iter_s.C_fit(iter,1) Iter_s.Curve{iter} call_with_ll]=loglikelihood_fit(call,k);
        
        sg=1; % var for while loop: 's'till 'g'oing
        slopes(1,1)=Iter_s.C_fit(iter,1);
        N_k(1,1)=sum(k);
        
        while sg
            disp(sprintf('iteration %d TinN fit: %d',iter,Iter_s.C_fit(iter)));
            
            if Iter_s.C_fit(iter)>=.35 %The mutation estimate is not reliable above this point and we should not return any estimate
                Iter_s.C_fit(iter)=-1;
                disp(sprintf('TiN might be too high to use mutations returning -1'));
                sg=0
                
            elseif Iter_s.C_fit(iter)<.35
                kk=(call.normal_f< ((Iter_s.C_fit(iter)+.2)*call.tumor_f) & call.tumor_f>.15 & ((call.tumor_f./call.normal_f) > 3));
                %((call.tumor_f>.8&call.normal_f<.2)|call.tumor_f<.8)&call.dont_include==0); %resampling points to include based on the previous TiN estimate
                Iter_s.nMut(iter)=length(kk);
                iter=iter+1;
                call_k=call;
                k=kk;
                if sum(k)==0
                    disp( 'no mutations remaining using previous iteration')
                    sg=0;
                    iter=iter-1;
                else
                    
                    [Iter_s.C_fit(iter,1) Iter_s.Curve{iter}]=loglikelihood_fit(call_k,k);
                    delta=Iter_s.C_fit(iter)-Iter_s.C_fit(iter-1);
                    sg=(iter < max_iter)&(f_c<abs(roundsig(delta,2)));
                    tin=round(Iter_s.C_fit(iter,1)*100)+2;
                    nk=find(k==1);
                    % remove outlier poins to avoid including SNPs in
                    % estimate
                    for i=1:length(nk)
                        c_d=betacdf(af_range,call.n_alt_count(nk(i))+1,call.n_ref_count(nk(i))+1);
                        if c_d(tin)<.05
                            call.dont_include(nk(i))=1;
                        end
                    end
                    disp(strcat(num2str(sum(call.dont_include)),' Mutations removed as outliers'))
                end
                TinNfit=Iter_s.C_fit(end);
                Iter_s.nMut(iter)=length(kk);
                
            end
            
            
        end
        
        %plotting TiN line
        hline=refline(Iter_s.C_fit(end),0);
        set(hline,'color','k')
        text(.1,.9,num2str(Iter_s.C_fit(end)));
        
        disp(strcat('estimated  ', num2str(TinNfit), ' tumor in normal'));
        
        
        %reintegrating low coverage sites for recovery
        call=all_calls;
        
        
        %calculate the log odds of each call belonging to the germline fit (regression line of non somatic events) or
        %the somatic fit (TiN estimate)
        disp('calculating the log odds of being somatic of each possible mutation')
        sites_to_consider=(ismember(call.judgement,'KEEP')|...
            ismember(call.failure_reasons,'germline_risk,normal_lod,alt_allele_in_normal')| ismember(call.failure_reasons,'normal_lod,alt_allele_in_normal') |...
            ismember(call.failure_reasons,'normal_lod,alt_allele_in_normal') | ismember(call.failure_reasons,'normal_lod') | ismember(call.failure_reasons,'alt_allele_in_normal') |...
            ismember(call.failure_reasons,'germline_risk') );
        for f=1:slength(call)
            %Only consider sites that have a greater tumorf than normalf
            %saves time.
            if call.normal_f(f)>call.tumor_f(f)||sites_to_consider(f)==0
                call.logOdds(f,1)=-Inf;
                call.prob(f,1)=0;
            else
                [call.logOdds(f,1) call.prob(f,1)]=normal_contamination_estimate_fitting(call.t_alt_count(f),call.t_ref_count(f)...
                    ,call.n_alt_count(f),call.n_ref_count(f),TinNfit,linfit.a,.001)  ;
            end
        end
        
        
        %Sites we want to recover must have only "germline allele" failure
        %reasons and have a logOdds greater than the threshold of being
        %somatic.
        normlod_filter=( call.logOdds>iter_threshold & (ismember(call.failure_reasons,'germline_risk,normal_lod,alt_allele_in_normal')| ...
            ismember(call.failure_reasons,'normal_lod,alt_allele_in_normal') | ismember(call.failure_reasons,'normal_lod') | ismember(call.failure_reasons,'alt_allele_in_normal') |...
            ismember(call.failure_reasons,'germline_risk')| ismember(call.failure_reasons,'alt_allele_in_normal')&~ismember(call.failure_reasons,',PoN')));
        nk=find(normlod_filter);
        af_range=[0:.01:1];
        TiN_ind=round(TinNfit*100);
        %for the recovered sites we remove outliers again. %NOTE: I SHOULD
        %FIX THIS TO STORE cdf dists in the call struct
        for i=1:length(nk)
            c_d=betacdf(af_range,call.n_alt_count(nk(i))+1,call.n_ref_count(nk(i))+1);
            if c_d(TiN_ind+1)<.05
                normlod_filter(nk(i))=0;
            end
            
        end
        
        % Generate 1 Sigma confidence intervals on the observed allele
        % fraction for plotting
        call.ci_l=zeros(slength(call),1);
        call.ci_h=zeros(slength(call),1);
        for nnk=1:length(nk)
            i=nk(nnk);
            alleles=binoinv([0.16 0.84],call.n_alt_count(i)+call.n_ref_count(i),call.normal_f(i));
            call.ci_l(i,1)=call.normal_f(i)-(alleles(1)/(call.n_alt_count(i)+call.n_ref_count(i)));
            call.ci_h(i,1)=call.normal_f(i)-(alleles(2)/(call.n_alt_count(i)+call.n_ref_count(i)));
        end
        %plot recovered points.
        errorbar(call.tumor_f(normlod_filter),call.normal_f(normlod_filter),call.ci_l(normlod_filter),call.ci_h(normlod_filter),'.','color',[0 0 0],'marker','o')
        % saving and writing files to upload in firehose
        if sum(normlod_filter)>0
            call.judgement(normlod_filter)={'KEEP'};
        end
        
        TinN(1,1)=TinNfit;
        call_post=sum(ismember(call.judgement,'KEEP'));
        save_struct(call,strcat(pair_id,'.call_stats.TiN.txt'));
        
        
        
        
        
        
        
        % plotting
        xlabel('Tumor Allele Fraction');
        ylabel('Normal Allele Fraction');
        title(pair_id);
        legend('Background','Germline','Sites Considered','Default Keep','Used in Inital Fit','TiN Fit','Recovered Calls')
        saveas(gcf,sprintf('NormalEst_%s.fig',pair_id));
        h=gcf;
        set(h,'PaperPositionMode','auto');
        % set(h,'PaperOrientation','landscape');
        f=sprintf('TiN_Model_%s',pair_id);
        %print(h,'-dpdf',sprintf('TiN_Model_%s.pdf',pair_id));
        saveas(gcf,[f '.png'],'png')
        print(gcf,'-depsc',[f '.eps'])
        
        save(sprintf('TinN_%s.txt',pair_id),'TinN','-ascii', '-double')
        save(sprintf('calls_pre_%s.txt',pair_id),'call_pre','-ascii', '-double');
        save(sprintf('calls_post_%s.txt',pair_id),'call_post','-ascii', '-double');
        if exist('Iter_s')
            figure()
            plot(Iter_s.Curve{end})
            title('Loglikelihood Curve')
            
            saveas(gcf,sprintf('LogLikeCurve_%s.png',pair_id),'png')
            
            % print(gcf,'-dpng',sprintf('LogLikeCurve_%s.png',pair_id));
            
            llcurve=Iter_s.Curve{end};
            
            TiNdist=[0:.005:1];
            llcurve2=-llcurve+(1-max(-llcurve));
            n_llcurve=exp(llcurve2)/sum(exp(llcurve2));
            Iter_s.CI_l=TiNdist(max([find(cumsum(n_llcurve)<0.16,1,'last') 1 ]));
            Iter_s.CI_h=TiNdist(min([find(cumsum(n_llcurve)>0.84,1,'first') 201 ]));
            %[Iter_s.CI_l,Iter_s.CI_h]=findCIp(n_llcurve,.05,.005);
            
            TiNCurve.mutations_pTiN=n_llcurve;
            save_struct(TiNCurve,'ll_curve_muts.txt');
            
            CIs=[Iter_s.CI_l,Iter_s.CI_h];
            save(sprintf('%s.CItxt',pair_id),'CIs','-ascii', '-double');
            save('deTiNMut.mat','Iter_s')
            
        else
            Iter_s=0;
        end
        %save('slopes.txt','slopes','-ascii','-double');
        %save('MutsUsed.txt','N_k','-ascii','-double');
        if isequal(firehose,'1')
            
            quit()
        end
    elseif isequal(mode,'force') %Estimate based on cannonical sites (developed for JMML)
        display('Using Force call mode')
        bed=load_struct_noheader(force_call_bed);
        call_pre=sum(ismember(call.judgement,'KEEP'));
        call_post=sum(ismember(call.judgement,'KEEP'));
        k=call.tumor_f>.05&ismember(call.position,bed.col2)&ismember(call.contig,bed.col1);
        if sum(k)>0
            iter=1;
            [Iter_s.C_fit(iter,1) Iter_s.Curve{iter}]=loglikelihood_fit(call,k);
            
            errorbar(call.tumor_f(k),call.normal_f(k),call.ci(k),'.','color',[0 0 0],'marker','o')
            llcurve=Iter_s.Curve{end};
            
            
            n_llcurve=exp(-llcurve)/sum(exp(-llcurve));
            TiNCurve.mutations_pTiN=n_llcurve';
            [Iter_s.CI_l,Iter_s.CI_h]=findCIp(n_llcurve,.05,.0005);
            CIs=[Iter_s.CI_l,Iter_s.CI_h];
            save(sprintf('%s.CItxt',pair_id),'CIs','-ascii', '-double');
            TinNfit=Iter_s.C_fit(end);
            TinN(1,1)=TinNfit;
            save_struct(TiNCurve,'ll_curve_muts.txt');
            save(sprintf('TinN_%s.txt',pair_id),'TinN','-ascii', '-double')
            save(sprintf('calls_pre_%s.txt',pair_id),'call_pre','-ascii', '-double');
            
            save(sprintf('calls_post_%s.txt',pair_id),'call_post','-ascii', '-double');
            save_struct(call,strcat(pair_id,'.call_stats.TiN.txt'));
        else
            
            call=all_calls;
            Iter_s=-1;
            TinN(1,1)=-1  ;
            call_pre=sum(ismember(call.judgement,'KEEP'));
            call_post=sum(ismember(call.judgement,'KEEP'));
            save_struct(call,strcat(pair_id,'.call_stats.TiN.txt'));
            save(sprintf('TinN_%s.txt',pair_id),'TinN','-ascii', '-double')
            save(sprintf('calls_pre_%s.txt',pair_id),'call_pre','-ascii', '-double');
            save(sprintf('calls_post_%s.txt',pair_id),'call_post','-ascii', '-double');
            save_struct(call,strcat(pair_id,'.call_stats.TiN.txt'));
            
        end
    end
else
    h=gcf;
    set(h,'PaperPositionMode','auto');
    %set(h,'PaperOrientation','landscape');
    f=sprintf('TiN_Model_%s',pair_id);
    %print(h,'-dpdf',sprintf('TiN_Model_%s.pdf',pair_id));
    saveas(gcf,[f '.png'],'png')
    print(gcf,'-depsc',[f '.eps'])
    call=all_calls;
    Iter_s=-1;
    TinN(1,1)=-1  ;
    call_pre=sum(ismember(call.judgement,'KEEP'));
    call_post=sum(ismember(call.judgement,'KEEP'));
    save_struct(call,strcat(pair_id,'.call_stats.TiN.txt'));
    save(sprintf('TinN_%s.txt',pair_id),'TinN','-ascii', '-double')
    save(sprintf('calls_pre_%s.txt',pair_id),'call_pre','-ascii', '-double');
    save(sprintf('calls_post_%s.txt',pair_id),'call_post','-ascii', '-double');
    
    if isequal(firehose,'1')
        
        quit()
    end
    
    
end

if isequal(firehose,'1')
    
    quit()
end
end


