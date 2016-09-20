TGCT_CN = readtable('~/Documents/CN_by_histology.txt');
hist_data = readtable('STable1_ForComp.txt');


for i=1:height(TGCT_CN)
    TGCT_CN.sample{i,1}=regexp(TGCT_CN.samples{i},'DFCI_[0-9]+','match');
    TGCT_CN.sample{i,1}=char(TGCT_CN.sample{i,1});
end


[i m]=ismember(hist_data.patient_id,TGCT_CN.sample);



UK = readtable('~/Documents/CN_by_histology_uk.txt');
UK_hist = readtable('~/Documents/UK_histology_map.txt');