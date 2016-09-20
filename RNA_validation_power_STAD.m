function maf=RNA_validation_power_STAD(MAF,DiscoveryMaf)

maf=load_struct(MAF);
Dmaf=load_struct(DiscoveryMaf);


min_count_of_reads=2;

for i=1:slength(maf)
    maf.barcode{i,1}=maf.sample_id{i}(1:end-4);
    maf.id{i,1}=strcat(maf.barcode{i},'-',maf.pos{i});
end

for i=1:slength(Dmaf)
    Dmaf.barcode{i,1}=Dmaf.patient{i}(6:17);
    Dmaf.id{i,1}=strcat(Dmaf.barcode{i},'-',Dmaf.Start_position{i});
end

for i=1:slength(maf)
    k=find(ismember(Dmaf.id,maf.id{i}));
    maf.discovery_alt_count(i,1)=str2num(Dmaf.alt_count{k});
    maf.discovery_coverage(i,1)=str2num(Dmaf.alt_count{k})+str2num(Dmaf.ref_count{k});
    maf.RNA_call{i,1}=Dmaf.RNA_verified{k};
end



for i=1:slength(maf)
    strs=regexp(maf.reads_by_letter{i},'"','split');
    vals=strs{2};
    vals=vals(2:end-1);
    vals=regexp(vals,',','split');
    maf.A(i,1)=str2num(vals{1});
    maf.C(i,1)=str2num(vals{2});
    maf.G(i,1)=str2num(vals{3});
    maf.T(i,1)=str2num(vals{4});
    
    maf.power_of_site(i,1)=1-hyge2cdf(min_count_of_reads-1,(maf.A(i)+maf.C(i)+maf.G(i)+maf.T(i)),...
        maf.discovery_alt_count(i),maf.discovery_coverage(i));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Code from Mara
%ms3.power_of_site(i)=1-hyge2cdf(ms3.min_val_count(i)-1,ms3.TTotCovVal(i),ms3.TVarCov(i),ms3.TTotCov(i));
%ms3.min_val_count would be 2 in your case. 
%This is the lower bound that I use for validating a site. MuTect has a theoretical power calculation where it uses 3 as the lower bound. 
%ms3.TTotCovVal(i) = tumor total coverage validation
%ms3.TVarCov(i) = allele count in discovery data
%ms3.TTotCov(i) = total coverage of discovery data
%%%%%%%%%%%%%%%%%%%%%%%%%%



end


function test

MAF='/Users/amaro/Downloads/FullSNPSvalRNA.txt';
DiscoveryMaf='/Users/amaro/Downloads/SigGeneMaf_RNAverified_RPKM_columns_20140204.tsv';


maf_with_validation=RNA_validation_power_STAD(MAF,DiscoveryMaf);


end
