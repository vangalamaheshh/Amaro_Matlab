function [ notes ms3 ] = power_validation_script_indels(pileup_file, maffile, sample_name, min_normal)

%addpath('/xchip/cga_home/mara/CancerGenomeAnalysis/trunk/matlab/');
%addpath('/xchip/cga_home/mara/CancerGenomeAnalysis/trunk/matlab/mike/');

%SETTING UP MATRIX STRUCTURE
msp = load_struct(pileup_file);
msp.ins(strcmp(msp.ins,'-1')) = {'0'};
msp.del(strcmp(msp.del,'-1')) = {'0'};
msp.ref(strcmp(msp.ref,'-1')) = {'0'};

z=repmat('|',length(msp.Chromosome),1);
msp.site = regexprep(cellstr([char(msp.Chromosome) z char(msp.Start_position)]),' ','');

maf_o = load_struct(maffile);
z = repmat('|', length(maf_o.Chromosome),1);
maf_o.site = regexprep(cellstr([char(maf_o.Chromosome) z char(maf_o.End_position)]),' ',''); %NOTE: the preprocessing step uses the end position rather than the start position


disc_maf_tumor = reorder_struct(msp, strcmp(msp.sample_type,'discovery_tumor'));
disc_maf_normal = reorder_struct(msp, strcmp(msp.sample_type,'discovery_normal'));
val_maf_tumor = reorder_struct(msp, strcmp(msp.sample_type, 'validation_tumor'));

if sum(strcmp(msp.sample_type, 'validation_normal'))==0
    display('no normal validation present. using discovery normal for evaluation of noise') 
    val_maf_normal = disc_maf_normal;
else
    val_maf_normal = reorder_struct(msp, strcmp(msp.sample_type,'validation_normal'));
end

ms.Chromosome = maf_o.Chromosome;
ms.Start_position = maf_o.Start_position;
ms.site = maf_o.site;
ms.Variant_Type = maf_o.Variant_Type;

ms.T_CountI = disc_maf_tumor.ins; %insertion
ms.T_CountD = disc_maf_tumor.del; %deletion
ms.T_CountR = disc_maf_tumor.ref; %reference
ms.N_CountI = disc_maf_normal.ins;
ms.N_CountD = disc_maf_normal.del;
ms.N_CountR = disc_maf_normal.ref;
ms.Tv_CountI = val_maf_tumor.ins;
ms.Tv_CountD = val_maf_tumor.del;
ms.Tv_CountR = val_maf_tumor.ref;
ms.Nv_CountI = val_maf_normal.ins;
ms.Nv_CountD = val_maf_normal.del;
ms.Nv_CountR = val_maf_normal.ref;


ms = make_numeric(ms,{'T_CountI','T_CountD','T_CountR',...
    'N_CountI','N_CountD','N_CountR',...
    'Tv_CountI','Tv_CountD','Tv_CountR',...
    'Nv_CountI','Nv_CountD','Nv_CountR'});
    
notes=struct();

%calc new fields    
ms.TTotCovVal=nan(slength(ms),1);
ms.NTotCovVal=nan(slength(ms),1);
ms.TVarCov=nan(slength(ms),1);
ms.TTotCov=nan(slength(ms),1);
ms.NTotCov=nan(slength(ms),1);

ms.v=nan(slength(ms),6);
ms.e=nan(slength(ms),6); % only tot and var counts in MAF, assume alt1=alt2=0

for i=1:slength(ms)    
    if strcmp(ms.Variant_Type(i),'INS')
        ms.v(i,:) = [ms.Tv_CountR(i) ms.Tv_CountI(i) ms.Tv_CountD(i)...
            ms.Nv_CountR(i) ms.Nv_CountI(i) ms.Nv_CountD(i)];
        ms.e(i,:) = [ms.T_CountR(i) ms.T_CountI(i) ms.T_CountD(i)...
            ms.N_CountR(i) ms.N_CountI(i) ms.N_CountD(i)];
    end
    
    if strcmp(ms.Variant_Type(i),'DEL')
        ms.v(i,:) = [ms.Tv_CountR(i) ms.Tv_CountD(i) ms.Tv_CountI(i)...
            ms.Nv_CountR(i) ms.Nv_CountD(i) ms.Nv_CountI(i)];
        ms.e(i,:) = [ms.T_CountR(i) ms.T_CountD(i) ms.T_CountI(i)...
            ms.N_CountR(i) ms.N_CountD(i) ms.N_CountI(i)];
    end
      
    ms.TTotCovVal(i)= sum(ms.v(i,1:3));
    ms.NTotCovVal(i)= sum(ms.v(i,4:6));

    ms.TVarCov(i,1) = ms.e(i,2);
    ms.NVarCov(i,1) = ms.e(i,5);
    ms.TTotCov(i,1) = sum(ms.e(i,1:3));
    ms.NTotCov(i,1) = sum(ms.e(i,4:6));
end

% test allele fraction the same
e_af = ms.e(:,2)./(ms.e(:,1)+ms.e(:,2)); %allele fraction of original
v_af = ms.v(:,2)./(ms.v(:,1)+ms.v(:,2)); %allele fraction of validation

if sum(~isnan(v_af))==0
    data=0;
    display('Zero reads covering the sites therefore no validation power at all! Choose to fill these columns of maf with NA')
    maf_o.validation_judgement = nan(slength(maf_o),1);
    maf_o.validation_power = nan(slength(maf_o),1);
else
    data=1;
    not_nan=~isnan(e_af) & ~isnan(v_af);
    h=plot(e_af(not_nan),v_af(not_nan),'.');
    add_x_equ_y_line
    title('Allele Fraction in Discovery and Validation')
    xlabel('Discovery AF')
    ylabel('Validation AF')
    saveas(h,'AlleleFractionDistribution','png')
    %TO DO: Fisher Exact test to see if really are similar
end

if data==1
    % focus on sites with >=min_normal reads in normal to avoid f.p due to germline
    min_normal = str2num(min_normal);
    e_n_ge_30=ms.NTotCov>=min_normal;
    notes.normal_sites_more_30_disc=nnz(e_n_ge_30);

    % should not be any site that is a germline het in this set ....
    ms2=reorder_struct(ms,find(e_n_ge_30)); %reodering to find all normals above coverage provided as input (30 is default)
    ms3=reorder_struct(ms2,find(ms2.TVarCov>0));
    ms3=reorder_struct(ms2,find(ms2.TTotCovVal>0));

    % calculate noise level at each site
    % calculate the number of observations that are very likely to be signal
    % for each site, look at validation normal and exp normal and calculate
    % the probability to see variants at that site (using laplace correction)
    % take maximal probability 

    ms3.p_e=[(ms3.e(:,1:3)+1)./repmat((sum(ms3.e(:,1:3),2)+3),1,3) ...
        (ms3.e(:,4:6)+1)./repmat((sum(ms3.e(:,4:6),2)+3),1,3)];  
    ms3.p_v=[(ms3.v(:,1:3)+1)./repmat((sum(ms3.v(:,1:3),2)+3),1,3) ...
        (ms3.v(:,4:6)+1)./repmat((sum(ms3.v(:,4:6),2)+3),1,3)];
    ms3.p_both=[(ms3.v(:,1:3)+ms3.e(:,1:3)+1)./repmat((sum(ms3.v(:,1:3)+ms3.e(:,1:3),2)+3),1,3) ...
        (ms3.v(:,4:6)+ms3.e(:,4:6)+1)./repmat((sum(ms3.v(:,4:6)+ms3.e(:,4:6),2)+3),1,3)];

    %ms3.pvar_max=max([ms3.p_e(:,6:8) ms3.p_v(:,6:8)],[],2); % look in normal but change so only coming from validation exp
    ms3.pvar_max=max(ms3.p_v(:,5:6),[],2); %pvalue max of variant - only look at vaalidation data because use different sequencing technologies
    ms3.min_val_count=nan(slength(ms3),1);

    %pval_norm_cont = str2num(pval_norm_cont);
    for i=1:slength(ms3)
    %    ms3.min_val_count(i)=binoinv(1-0.05/slength(ms3),ms3.TTotCovVal(i),ms3.pvar_max(i));
        ms3.min_val_count(i)=binoinv(0.99,ms3.TTotCovVal(i),ms3.pvar_max(i)); %minval count is the minimum number of variants need to see at that site?
        if ms3.min_val_count(i) < 2
            ms3.min_val_count(i)=2; %won't validate with anything less than 2 reads.
        end
    end

    vc=plot(ms3.min_val_count,ms3.v(:,2),'.');
    %add_x_equ_y_line

    hold on
    x=1:max(ms3.v(:,2));
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
        %ms3.power_of_site(i) = 1 - binocdf(ms3.min_val_count(i)-1, ms3.TTotCovVal(i), ms3.TVarCov(i)/ms3.TTotCov(i));
    end

    ms3.power_of_site(ms3.TTotCovVal==0)=0; %If have no covarage in the validation tumor, then force it to say have no power!
    ms3.power_of_site(isnan(ms3.TTotCovVal))=0; 

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

    ms3.validated(sum(ms3.TTotCovVal)==0)=0;  %If have no covarage in the validation tumor, then force it to not validate!
    ms3.validated(isnan(ms3.TTotCovVal))=0; 


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
end

maf_o = rmfield(maf_o,'site');
save_struct(maf_o, [sample_name, '.indel.maf.annotated.validated']);

exit
end

