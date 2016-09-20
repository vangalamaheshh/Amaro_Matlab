%prototype of LOH + Mutation De:Tin
% sample with nearly 100% TiN
seg=load_struct('/Users/amaro/Downloads/CLL-GCLL-0056-Tumor-SM-41JNK.tsv');
% %Call Stats File Code:
% call=load_struct('/Users/amaro/Downloads/CLL-GCLL-0056-TP-NT-SM-41JNK-SM-41QD2.call_stats.txt');
% call.t_ref_count=str2double(call.t_ref_count);
% call.t_alt_count=str2double(call.t_alt_count);
% call.n_ref_count=str2double(call.n_ref_count);
% call.n_alt_count=str2double(call.n_alt_count);
% call.tumor_f=str2double(call.tumor_f);
% call.normal_f=str2double(call.normal_f);
% call.position=str2double(call.position);
% call.contig=chromosome2num(call.contig);
% call.init_n_lod=str2double(call.init_n_lod);
% call.init_t_lod=str2double(call.init_t_lod);
% call=reorder_struct(call,~isnan(call.contig));
% call=reorder_struct(call,((call.t_alt_count+call.t_ref_count)>15)&((call.n_alt_count+call.n_ref_count)>15));


%symmetry filter mean (af(hets)) stderr(af(hets) 50% within stderr

%Germline Maf Code
m=load_struct('/Users/amaro/Downloads/CLL-GCLL-0056-Normal-SM-41QD2.maf');
m=reorder_struct(m,ismember(m.genotype,'0/1'));
m=reorder_struct(m,~ismember(m.allelic_depth,''));

call.contig=chromosome2num(m.Chromosome);
call.position=str2double(m.Start_position);
for i=1:slength(m)
acs=split(m.allelic_depth{i},',');
call.n_ref_count(i,1)=str2double(acs{1});
call.n_alt_count(i,1)=str2double(acs{2});
end

seg.f=str2double(seg.f);
seg.Startbp=str2double(seg.Startbp);
seg.Endbp=str2double(seg.Endbp);
seg.Chromosome=chrom2num(seg.Chromosome);
seg.length=str2double(seg.length);
seg.muminor=str2double(seg.muminor);
seg.mumajor=str2double(seg.mumajor);
seg=reorder_struct(seg,~isnan(seg.f));
seg=reorder_struct(seg,~(seg.f==0));
seg=reorder_struct(seg,(seg.length>1*10^6));
seg.sigmamajor=str2double(seg.sigmamajor);
seg.W=seg.length./(seg.sigmamajor.^2); % using sigmas from allelic
%capseg as theyre: stdev(allelic shift probability * copy ratio data)??
seg.delta=seg.mumajor-seg.muminor;
af=[0:.01:1];



for i=1:slength(seg)
    s=reorder_struct(call,(call.position>seg.Startbp(i)&call.position<seg.Endbp(i)&call.contig==seg.Chromosome(i)));
    seg.nHet(i,1)=slength(s);
    seg.meanNcov(i,1)=mean(s.n_alt_count+s.n_ref_count);
    afbeta=zeros(length(af),1)';
    for j=1:slength(s)
        j_beta=betapdf(af,s.n_alt_count(j)+1,s.n_ref_count(j)+1); %normalize
        afbeta=afbeta+(j_beta/sum(j_beta));
    end
   
    j=2:(length(afbeta)-1);
    kp=af(find((afbeta(j-1)<=afbeta(j))&(afbeta(j+1)<afbeta(j))));
    kp=kp(kp>.20&kp<.95);
    seg.normal_mumajor(i,1)=1+(max(kp)-.5);
    seg.normal_delta(i,1)=(seg.normal_mumajor(i)-1)*2;
    
end


 g = fittype('a*x');
 linfit=fit(seg.delta,seg.normal_delta,g,'Weights',seg.W);
 hold on
 BubblePlot(seg.delta,seg.normal_delta,seg.length,[0 0 1],300)
 hline=refline(linfit.a,0);
 set(hline,'Color',[0 0 0]);