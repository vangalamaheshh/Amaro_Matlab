% Test suite for seg_data class.
function test_suite = test_seg_data
  initTestSuite;
end


function test_create_object

    sd = seg_data();
    assertFalse (isempty(sd));
end

function test_set_base_dir
    sd = seg_data();
    basedir = sd.basedir;
    
    if ispc
        expected_basedir = '//thumper12/xchip_tcga_scratch/DCC_downloads';
    else
        expected_basedir = '/xchip/tcga_scratch/DCC_downloads';
    end
    
   assertEqual(basedir,expected_basedir);
  
end

function test_get_latest_open_dir_name
    sd = seg_data();
    
    %TODO perhaps recode to grab latest directory rather than hardcoded, or
    %point to test directory structure.
    
    opendir = sd.get_latest_open_dir();
    
    expected_opendir_addon = '20090415-open';
    expected_opendir = [sd.basedir '/' expected_opendir_addon];
      
    assertFalse (strcmp(opendir,expected_opendir));
    
    
    sd.cutoff_date = 20090415;
    opendir = sd.get_latest_open_dir();
    assertEqual (opendir,expected_opendir);
end

function test_get_center_dir
    sd = seg_data();
    center_dir = sd.get_center_dir('harvard',9);
    
    expected_center_dir_addon = ['tcga/tumor/ov/cgcc/' ...
        'hms.harvard.edu/hg-cgh-244a/cna'];
    expected_center_dir = [sd.get_latest_open_dir '/' ...
        expected_center_dir_addon];
    
    assertEqual(center_dir, expected_center_dir);

end

function test_get_batch_dir_name
    sd = seg_data();
    
    batch_dir_name1 = sd.get_batch_dir_name('harvard',9);
    
    expected_batch_dir_name_addon1 = 'hms.harvard.edu_OV.HG-CGH-244A.1.0.0';
    expected_batch_dir_name1 = [sd.get_center_dir('harvard',9) '/' ...
        expected_batch_dir_name_addon1];
    assertEqual (batch_dir_name1,expected_batch_dir_name1);    
        
        
    assertExceptionThrown(@() sd.get_batch_dir_name('harvard',10), ...
        'seg_data:badHarvardbatchnum');

    
    batch_dir_name2 = sd.get_batch_dir_name('harvard',11);

    expected_batch_dir_name_addon2 = 'hms.harvard.edu_OV.HG-CGH-244A.2.2.0';
    expected_batch_dir_name2 = [sd.get_center_dir('harvard',11) '/' ...
        expected_batch_dir_name_addon2];
    assertEqual (batch_dir_name2,expected_batch_dir_name2);    

end


function test_get_batch_files

sd = seg_data();
hms_files = sd.get_batch_files('harvard',9);
hms_num_files = size(hms_files,2);
hms_num_files_expected = 87;
assertEqual(hms_num_files,hms_num_files_expected);
hms_first_file = hms_files{1};
hms_first_file_expected = [sd.get_latest_open_dir,'/',...
    'tcga/tumor/ov/cgcc/hms.harvard.edu/hg-cgh-244a/cna/hms.harvard.edu_OV.HG-CGH-244A.1.0.0/TCGA-01-0628-11A-01D-0360-02_Segment.tsv'];
assertEqual(hms_first_file, hms_first_file_expected);
hms_second_file = hms_files{2};
hms_second_file_expected = [sd.get_latest_open_dir,'/',...
    'tcga/tumor/ov/cgcc/hms.harvard.edu/hg-cgh-244a/cna/hms.harvard.edu_OV.HG-CGH-244A.1.0.0/TCGA-01-0629-11A-01D-0360-02_Segment.tsv'];
assertEqual(hms_second_file, hms_second_file_expected);

mskcc_files = sd.get_batch_files('mskcc',9);
mskcc_num_files = size(mskcc_files,2);
mskcc_num_files_expected = 1;
mskcc_files_expected = {[sd.get_latest_open_dir,'/',...
    'tcga/tumor/ov/cgcc/mskcc.org/hg-cgh-244a/cna/mskcc.org_OV.HG-CGH-244A.9.3.0/mskcc.org_OV.HG-CGH-244A.9.CBS.txt']};
assertEqual(mskcc_files,mskcc_files_expected);


stanford_files = sd.get_batch_files('stanford',9);
stanford_num_files = size(stanford_files,2);
stanford_num_files_expected = 2;
assertEqual(stanford_num_files, stanford_num_files_expected);
stanford_files_expected1 = [sd.get_latest_restricted_dir,'/',...
    'users/tcga4yeo/tumor/ov/cgcc/hudsonalpha.org/human1mduo/snp/hudsonalpha.org_OV.Human1MDuo.1.1.0/hudsonalpha.org_OV.Human1MDuo.1.1.0.seg.txt'];
assertEqual(stanford_files{1},stanford_files_expected1);
stanford_files_expected2 = [sd.get_latest_restricted_dir,'/',...
    'users/tcga4yeo/tumor/ov/cgcc/hudsonalpha.org/human1mduo/snp/hudsonalpha.org_OV.Human1MDuo.1.1.0/hudsonalpha.org_OV.Human1MDuo.1.1.0.segnormal.txt'];
assertEqual(stanford_files{2},stanford_files_expected2);


broad_files = sd.get_batch_files('broad',9);
broad_num_files = size(broad_files,2);
broad_num_files_expected = 77;
assertEqual(broad_num_files,broad_num_files_expected);
broad_first_file = broad_files{1};
broad_first_file_expected = [sd.get_latest_restricted_dir,'/',...
    'users/tcga4yeo/tumor/ov/cgcc/broad.mit.edu/genome_wide_snp_6/snp/broad.mit.edu_OV.Genome_Wide_SNP_6.9.5.0/SampleInfo.txt'];
assertEqual(broad_first_file, broad_first_file_expected);
broad_second_file = broad_files{2};
broad_second_file_expected = [sd.get_latest_restricted_dir,'/',...
    'users/tcga4yeo/tumor/ov/cgcc/broad.mit.edu/genome_wide_snp_6/snp/broad.mit.edu_OV.Genome_Wide_SNP_6.9.5.0/COTES_p_TCGAaffxB8_9a_S_GenomeWideSNP_6_F01_293016.seg.data.txt'];
assertEqual(broad_second_file, broad_second_file_expected);
end

%long running
function test_get_batch_data_harvard
    sd = seg_data();
    %sd.cutoff_date = 20090415;
    hms_batch9 = sd.get_batch_data('harvard',9);
   
    if ispc
        test_base_dir = '//thumper12/xchip_tcga_scratch/gsaksena';
    else
        test_base_dir = '/xchip/tcga';
    end
    test_sub_dir = 'CancerGenomeAnalysisData/trunk/test/seg_data';
    test_filename = 'harvard_seg_batch9.tsv.CBS.txt';
    test_file = [test_base_dir '/' test_sub_dir '/' test_filename];
 
    hms_batch9_expected = dataset ('file',test_file,...
        'format','%s%f%f%f%f%f','delimiter','\t',...
        'treatasempty',{'NA','na'});
    
    %test file may have been processed in a different order, so sort.
    hms_batch9 = sortrows(hms_batch9);
    hms_batch9_expected = sortrows (hms_batch9_expected);
    
    assertEqual(hms_batch9,hms_batch9_expected);
    
end

function test_filter_for_tumors
    if ispc
        test_base_dir = '//thumper12/xchip_tcga_scratch/gsaksena';
    else
        test_base_dir = '/xchip/tcga';
    end
    test_sub_dir = 'CancerGenomeAnalysisData/trunk/test/seg_data';
    test_filename = 'harvard_seg_batch9.tsv.CBS.txt';
    test_file = [test_base_dir '/' test_sub_dir '/' test_filename];
 
    hms_batch9 = dataset ('file',test_file,...
        'format','%s%f%f%f%f%f','delimiter','\t',...
        'treatasempty',{'NA','na'});
    
    hms_batch9_filtered = seg_data.filter_keep_tumors_only_ds(hms_batch9);
    
    samplenames = hms_batch9_filtered.Sample;
    num_samples = size(samplenames,1);
    
    %loop thru cell, slower method than in real code, but different.
    for i=1:num_samples
        tumor_codes{i}=samplenames{i}(14:15);
    end
    num_tumors = size(strmatch('01',tumor_codes),1);
    
    assertEqual(num_samples,num_tumors);
    
end

function test_filter_outliers
    if ispc
        test_base_dir = '//thumper12/xchip_tcga_scratch/gsaksena';
    else
        test_base_dir = '/xchip/tcga';
    end
    test_sub_dir = 'CancerGenomeAnalysisData/trunk/test/seg_data';
    test_filename = 'harvard_seg_batch9.tsv.CBS.txt';
    test_file = [test_base_dir '/' test_sub_dir '/' test_filename];
    
    expected_filename = 'hms_batch9_postqc.txt';
    expected_file = [test_base_dir '/' test_sub_dir '/' expected_filename];

 
    hms_batch9 = dataset ('file',test_file,...
        'format','%s%f%f%f%f%f','delimiter','\t',...
        'treatasempty',{'NA','na'});
    
    
    hms_batch9_filtered = ...
        seg_data.filter_remove_excessively_segmented_samples_ds(hms_batch9);
    
    hms_batch9_filtered_expected = dataset ('file',expected_file,...
        'format','%s%f%f%f%f%f','delimiter','\t',...
        'treatasempty',{'NA','na'});
    
    assertEqual(hms_batch9_filtered, hms_batch9_filtered_expected);
end

%long running
function test_export

    if ispc
        test_base_dir = '//thumper12/xchip_tcga_scratch/gsaksena';
    else
        test_base_dir = '/xchip/tcga';
    end
    test_sub_dir = 'CancerGenomeAnalysisData/trunk/test/seg_data';
    test_filename = 'harvard_seg_batch9.tsv.CBS.txt';
    test_file = [test_base_dir '/' test_sub_dir '/' test_filename];
 
    hms_batch9 = dataset ('file',test_file,...
        'format','%s%f%f%f%f%f','delimiter','\t',...
        'treatasempty',{'NA','na'});

    out_filename = ['test_export_' date '.txt'];
    out_file = [test_base_dir '/' test_sub_dir '/' out_filename];
    
    seg_data.export_ds(hms_batch9,out_file);
    
    hms_batch9_reimport = dataset ('file',out_file,...
        'format','%s%f%f%f%f%f','delimiter','\t',...
        'treatasempty',{'NA','na'});
    
    assertEqual(hms_batch9_reimport, hms_batch9);
    
end

%long running
function test_broad

sd = seg_data();
sd.cutoff_date = 20090415;
%batch_dir_name = get_batch_dir_name(obj,center, batch)
broad_dir9 = get_batch_dir_name(sd,'broad', 9);
broad_dir9_expected = [sd.basedir '/' ...
    '20090415-restricted/users/tcga4yeo/tumor/ov/cgcc/broad.mit.edu/genome_wide_snp_6/snp/broad.mit.edu_OV.Genome_Wide_SNP_6.9.5.0'];

assertEqual(broad_dir9,broad_dir9_expected);





batch_data = get_batch_data(sd,'broad',9);
num_points = size(batch_data,1);

assertEqual(num_points,67463);

end

function test_mskcc

sd = seg_data();
sd.cutoff_date = 20090415;
%batch_dir_name = get_batch_dir_name(obj,center, batch)
mskcc_dir9 = get_batch_dir_name(sd,'mskcc', 9);
mskcc_dir9_expected = [sd.basedir '/' ...
    '20090415-open/tcga/tumor/ov/cgcc/mskcc.org/hg-cgh-244a/cna/mskcc.org_OV.HG-CGH-244A.9.3.0'];

assertEqual(mskcc_dir9,mskcc_dir9_expected);

mskcc_dir11 = get_batch_dir_name(sd,'mskcc', 11);
mskcc_dir11_expected = [sd.basedir '/' ...
    '20090415-open/tcga/tumor/ov/cgcc/mskcc.org/cgh-1x1m_g4447a/cna/mskcc.org_OV.CGH-1x1M_G4447A.11.1.0'];
assertEqual(mskcc_dir11,mskcc_dir11_expected);

%TODO test totalfile length

batch_data = get_batch_data(sd,'mskcc',9);
num_points_batch9 = size(batch_data,1);
num_points_batch9_expected = 17489-1;
assertEqual(num_points_batch9, num_points_batch9_expected);

batch_data = get_batch_data(sd,'mskcc',11);
num_points_batch11 = size(batch_data,1);
num_points_batch11_expected = 38653-1;
assertEqual(num_points_batch11, num_points_batch11_expected);


end

function test_stanford

sd = seg_data();
sd.cutoff_date = 20090415;
%batch_dir_name = get_batch_dir_name(obj,center, batch)
stan_dir9 = get_batch_dir_name(sd,'stanford', 9);
stan_dir9_expected = [sd.basedir '/' ...
    '20090415-restricted/users/tcga4yeo/tumor/ov/cgcc/hudsonalpha.org/human1mduo/snp/hudsonalpha.org_OV.Human1MDuo.1.0.0'];
assertEqual(stan_dir9,stan_dir9_expected);

batch_data = get_batch_data(sd,'stanford',9);
num_points = size(batch_data,1);
num_points_expected = 96203;
assertEqual(num_points,num_points_expected);

%TODO test total length

end
