function [ notes ms3 ] = validation_powered_script(pileup_file, maffile, sample_name, p_min_alt, min_normal_coverage)

%addpath('/xchip/cga2/mara/CancerGenomeAnalysis/CancerGenomeAnalysis/trunk/matlab/');
%addpath('/xchip/cga2/mara/CancerGenomeAnalysis/CancerGenomeAnalysis/trunk/matlab/mike/');

ms = load_struct(pileup_file);
%b = char2cell(nan(slength(ms),1));
%bb = [ms.Chromosome b ms.Start_position];

z=repmat('|',length(ms.Chromosome),1);
ms.site = regexprep(cellstr([char(ms.Chromosome) z char(ms.Start_position)]),' ','');
%ms.site = arrayfun(@(r) [bb{r,:}], 1:size(bb,1),'Unif',false)';

maf_o = load_struct(maffile);
%a = char2cell(nan(slength(maf_o),1));
%aa = [maf_o.Chromosome a maf_o.Start_position];
%maf_o.site = arrayfun(@(r) [aa{r,:}], 1:size(aa,1),'Unif',false)';
z = repmat('|', length(maf_o.Chromosome),1);
maf_o.site = regexprep(cellstr([char(maf_o.Chromosome) z char(maf_o.Start_position)]),' ','');

[o ind1 ind2] = intersect(ms.site, maf_o.site);
ms = reorder_struct(ms, ind1);
maf = reorder_struct(maf_o, ind2);

[overlap ind1 ind2] = intersect(ms.site, maf.site, 'stable');

ms = reorder_struct(ms, ind1);
maf = reorder_struct(maf, ind2);

ms=make_numeric(ms,{'T_CountA','T_CountC','T_CountG','T_CountT',...
        'N_CountA','N_CountC','N_CountG','N_CountT',...
        'Tv_CountA','Tv_CountC','Tv_CountG','Tv_CountT',...
        'Nv_CountA','Nv_CountC','Nv_CountG','Nv_CountT'});
   
notes=struct();

%calc new fields    
ms.val=nan(slength(ms),8);
ms.ref=nan(slength(ms),1);
ms.var=nan(slength(ms),1);
ms.ord=nan(slength(ms),4);
ms.TTotCovVal=nan(slength(ms),1);
ms.NTotCovVal=nan(slength(ms),1);
ms.TVarCov=nan(slength(ms),1);
ms.TTotCov=nan(slength(ms),1);
ms.NTotCov=nan(slength(ms),1);

ms.v=ms.val;
ms.e=nan(slength(ms),8); % only tot and var counts in MAF, assume alt1=alt2=0
basetab=zeros(255,1);
basetab('A')=1;
basetab('C')=2;
basetab('G')=3;
basetab('T')=4;

notes.variant_more_thean_one_tumor_seq_allele1 = 0
for i=1:slength(ms)
    ms.ref(i)=basetab(ms.Reference_Allele{i}); %applying numbers for reference
    if size(basetab(maf.Tumor_Seq_Allele1{i}),1)>1
        display('more than one tumor_seq_allele1. skipping');
        notes.variant_more_thean_one_tumor_seq_allele1=notes.variant_more_thean_one_tumor_seq_allele1 + 1;
    else
        ms.var(i)=basetab(maf.Tumor_Seq_Allele1{i}); %applying numbers for tumor allele2 (alternate) Should this be Tumor_Seq_Allele1?
        if ms.ref(i)==ms.var(i), keyboard; error('cannot be 1');  end

        ms.val(i,:)=[ms.Tv_CountA(i) ms.Tv_CountC(i) ms.Tv_CountG(i) ms.Tv_CountT(i) ...
           ms.Nv_CountA(i) ms.Nv_CountC(i) ms.Nv_CountG(i) ms.Nv_CountT(i)]; %valuation values for alleles in T and normal
        ms.TTotCovVal(i)= sum(ms.val(i,1:4));
        ms.NTotCovVal(i)= sum(ms.val(i,5:8));

        s=setdiff(1:4,[ms.ref(i) ms.var(i)]); %alleles that are NOT present in either refrence or allele2
        if length(s)~=2, keyboard; error('cannot be 2'); end;

        summed_counts=ms.val(i,1:4)+ms.val(i,5:8); %Count all alleles present in combined tumor and normal
        [sorted_counts,sort_ord]=sort(summed_counts(s),2,'descend'); %Sorts the counts for the alelles NOT in reference or allele2

        ms.ord(i,:)=[ ms.ref(i) ms.var(i) s(sort_ord)];

        ms.v(i,:)=ms.val(i,[ms.ord(i,:) 4+ms.ord(i,:)]); %reference count, variant count, alt1 count, alt2 count

        ms.exp(i,:)=[ms.T_CountA(i) ms.T_CountC(i) ms.T_CountG(i) ms.T_CountT(i) ...
           ms.N_CountA(i) ms.N_CountC(i) ms.N_CountG(i) ms.N_CountT(i)]; %valuation values for alleles in T and normal

        ms.e(i,:)= ms.exp(i,[ms.ord(i,:) 4+ms.ord(i,:)]);
        ms.TVarCov(i,1) = ms.e(i,2);
        ms.TTotCov(i,1) = sum(ms.exp(i,1:4));
        ms.NTotCov(i,1) = sum(ms.exp(i,5:8));
        ms.TVarCovVal(i,1) = sum(ms.v(i,1:4));
    end
    
end

% test allele fraction the same

e_af = ms.e(:,2)./(ms.e(:,1)+ms.e(:,2)); %allele fraction of original
v_af = ms.v(:,2)./(ms.v(:,1)+ms.v(:,2)); %allele fraction of validation
if sum(~isnan(v_af))==0
    display('Zero reads covering the sites therefore infinite allele fraction for everything and no validation power at all!')
else
not_nan=~isnan(e_af) & ~isnan(v_af);
h=plot(e_af(not_nan),v_af(not_nan),'.');
add_x_equ_y_line
title('Allele Fraction in Discovery and Validation')
xlabel('Discovery AF')
ylabel('Validation AF')
saveas(h,'AlleleFractionDistribution','png')
%TO DO: Fisher Exact test to see if really are similar
end

% focus on sites with >=30 reads in normal to avoid f.p due to germline 
min_normal_coverage = str2num(min_normal_coverage);
e_n_ge_30=ms.NTotCov>=min_normal_coverage;

notes.normal_sites_more_30_disc=nnz(e_n_ge_30)

% find good validation data  val: t>=30, n>=30

good_val = sum(ms.v(:,1:4),2)>=30 & sum(ms.v(:,5:8),2)>=30;

notes.good_val_Tmore30_Nmore30=nnz(good_val)

notes.min_good_val_tumor_cov=min(ms.TTotCovVal)  
notes.min_good_val_normal_cov_val=min(ms.NTotCovVal) 

% should not be any site that is a germline het in this set ....

ms2=reorder_struct(ms,find(e_n_ge_30)); %reodering to find all normals above coverage of 30

very_unlikely=find(ms2.v(:,6)>(ms2.v(:,5)-3*sqrt(ms2.v(:,5))));  % require higher than 3 stds below the ref count %variant in normal greater than (reference normal - 3*sqrt(refrence normal)
% let's look at lower normals (below 30)
check_all=find(ms.v(:,6)>(ms.v(:,5)-3*sqrt(ms.v(:,5))))

notes.very_unlikely=size(very_unlikely,1)

% remove these ones. should have only atrifacts at this point

ms3=reorder_struct(ms2,setdiff(1:slength(ms2),very_unlikely));
ms3=reorder_struct(ms3,find(ms3.TVarCov>0)); 

% calculate noise level at each site
% calculate the number of observations that are very likely to be signal
% for each site, look at validation normal and exp normal and calculate
% the probability to see variants at that site (using laplace correction)
% take maximal probability 

ms3.p_e=[(ms3.e(:,1:4)+1)./repmat((sum(ms3.e(:,1:4),2)+4),1,4) ...
    (ms3.e(:,5:8)+1)./repmat((sum(ms3.e(:,5:8),2)+4),1,4)];  
ms3.p_v=[(ms3.v(:,1:4)+1)./repmat((sum(ms3.v(:,1:4),2)+4),1,4) ...
    (ms3.v(:,5:8)+1)./repmat((sum(ms3.v(:,5:8),2)+4),1,4)];
ms3.p_both=[(ms3.v(:,1:4)+ms3.e(:,1:4)+1)./repmat((sum(ms3.v(:,1:4)+ms3.e(:,1:4),2)+4),1,4) ...
    (ms3.v(:,5:8)+ms3.e(:,5:8)+1)./repmat((sum(ms3.v(:,5:8)+ms3.e(:,5:8),2)+4),1,4)];

%ms3.pvar_max=max([ms3.p_e(:,6:8) ms3.p_v(:,6:8)],[],2); % look in normal but change so only coming from validation exp

ms3.pvar_max=max(ms3.p_v(:,6:8),[],2) %pvalue max of variant - only look at vaalidation data because use different sequencing technologies

% ms3.pvar_max=max([ms3.p_both(:,6:8)],[],2); % look in normal

p_min_alt = str2num(p_min_alt)

ms3.min_val_count=nan(slength(ms3),1);
for i=1:slength(ms3)
%    ms3.min_val_count(i)=binoinv(1-0.05/slength(ms3),ms3.TTotCovVal(i),ms3.pvar_max(i));
    ms3.min_val_count(i)=binoinv(p_min_alt,ms3.TTotCovVal(i),ms3.pvar_max(i)); %minval count is the minimum number of variants need to see at that site 
    if ms3.min_val_count(i)<2
        ms3.min_val_count(i)=2;
    end
end
%normally use 2

vc=plot(ms3.min_val_count,ms3.v(:,2),'.');
%add_x_equ_y_line

hold on
x=1:max(ms3.v(:,2))
plot(x,x)
axis([0 100 0 2000])
title('Minimum Number Alleles Required to be Seen in Validation at Different Variant Counts')
ylabel('Number of Reads Supporting the Variant in Tumor')
xlabel('Minimum Number of Supporting Reads Required for Validation')
saveas(vc,'MinValidationCount','png')
legend('Calls','Calls above or on this line (x=y) would be validated')

% calculate sites with at least 99% power to have min_val_count or more

ms3.power_of_site=nan(slength(ms3),1);
for i=1:slength(ms3)
    ms3.power_of_site(i)=1-hyge2cdf(ms3.min_val_count(i)-1,ms3.TTotCovVal(i),ms3.TVarCov(i),ms3.TTotCov(i)); 
end

ms3.power_of_site(ms3.TVarCovVal==0)=0; %If have no covarage in the validation tumor, then can't have any power.
ms3.power_of_site(isnan(ms3.TVarCovVal))=0; 

[mn,mni]=min(ms3.power_of_site)

powered_sites=ms3.power_of_site>=0.99;

notes.number_power_sites=nnz(powered_sites) % 1483 %999
notes.percent_powered_sites=nnz(powered_sites)/slength(ms3) % 91.88 %92.33


% Graphing Power Data:
ms3.e_af = ms3.e(:,2)./(ms3.e(:,1)+ms3.e(:,2)); %allele fraction of original
ms3.v_af = ms3.v(:,2)./(ms3.v(:,1)+ms3.v(:,2)); %allele fraction of validation

%plot3(ms3.TTotCov,ms3.e_af, ms3.power_of_site,'.')
figure(5);
a = plot3(ms3.TTotCovVal(ms3.power_of_site>=0.99),ms3.v_af(ms3.power_of_site>=0.99), ms3.power_of_site(ms3.power_of_site>=0.99),'.b');
hold on
plot3(ms3.TTotCovVal(ms3.power_of_site<0.99),ms3.v_af(ms3.power_of_site<0.99), ms3.power_of_site(ms3.power_of_site<0.99),'.r');
title('Power to Validate by Allelic Fraction and Tumor Depth in Validation Data')
xlabel('Tumor Depth');
ylabel('Allele Fraction');
zlabel('Power');
legend('Powered (>=0.99)','Unpowered (<0.99)','Location','NorthEast')
grid
%saveas(a,'PowerAFDepth', 'png')
print_D('PowerAFDepth',{{'fig'},{'pdf'},{'png','-r180'}});

ms3.validated=ms3.v(:,2)>=ms3.min_val_count;

ms3.validated(ms3.TVarCovVal==0)=0;  %If have no covarage in the validation tumor, then force it to not validate!
ms3.validated(isnan(ms3.TVarCovVal))=0; 


notes.percent_validated = nnz(ms3.validated)/slength(ms3);

% Write Out Data
[o ind_maf_o ind_maf] = intersect(maf_o.site, ms3.site);
maf_o.validation_judgement = nan(slength(maf_o),1);
maf_o.validation_judgement(ind_maf_o) = ms3.validated(ind_maf);
maf_o.validation_power = nan(slength(maf_o),1);
maf_o.validation_power(ind_maf_o) = ms3.power_of_site(ind_maf);
maf_o.validation_tumor_alt_count = nan(slength(maf_o),1);
maf_o.validation_tumor_alt_count(ind_maf_o,1) = ms3.v(ind_maf,2);
maf_o.validation_tumor_ref_count = nan(slength(maf_o),1);
maf_o.validation_tumor_ref_count(ind_maf_o,1) = ms3.v(ind_maf,1);

save_struct(notes,'notes.txt')
save_struct(ms3,'coverage_detail.txt')

maf_o = rmfield(maf_o,'site');
save_struct(maf_o, [sample_name, '.oxoG3.maf.annotated.validated']);

end
