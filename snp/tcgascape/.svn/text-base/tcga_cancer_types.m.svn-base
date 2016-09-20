function [cancer_name cancer_tree] = tcga_cancer_types()
% TCGA_CANCER_TYPES return mappings of TCGA cancer types
%
%   [CANCER_NAME CANCER_TREE] = TCGA_CANCER_TYPES()
%
%  CANCER_NAME is a containers.Map from TCGA-defined names and superclass
%  names to the strings used in the Tumorscape data. CANCER_TREE is a
%  containers.Map of superclasses to their immediate subclasses in the
%  cancer type hierarchy, using TCGA-defined tissue type names for the leaf
%  nodes. The root name in the hierarch is 'all_cancers'.

%% define classes of cancer with "generative grammar"
cancer_tree = containers.Map;
cancer_tree('all_cancers')    = {'all_epithelial','all_neural','all_blood'};
cancer_tree('all_epithelial') = {'all_colorectal','all_lung','all_kidney',...
                                 'stad','ov','ucec','brca','hnsc','thca',...
                                 'cesc','blca','paad','prad','lihc'};
cancer_tree('all_colorectal') = {'coad','read'};
cancer_tree('all_kidney')     = {'kirc','kirp'};
cancer_tree('all_lung')       = {'luad','lusc'};
cancer_tree('all_neural')     = {'gbm','lgg'};
cancer_tree('all_blood')      = {'laml','dlbc'};

%% map TCGA abreviations to text displayed in Tumorscape
cancer_name = containers.Map;
cancer_name('all_cancers') =    'all_cancers';
cancer_name('all_epithelial') = 'Epithelial';
cancer_name('all_kidney') =     'Kidney';
cancer_name('all_lung') =       'Lung';
cancer_name('all_colorectal') = 'Colorectal';
cancer_name('all_neural') =     'Neural';
cancer_name('all_blood') =      'Blood';
%-------------
cancer_name('blca') =           'Bladder urothelial carcinoma';
cancer_name('brca') =           'Breast invasive adenocarcinoma';
cancer_name('cesc') =           'Cervical squamous cell carcinoma';
cancer_name('coad') =           'Colon adenocarcinoma';
cancer_name('dlbc') =           'Diffuse large B-cell lymphoma';
cancer_name('gbm') =            'Glioblastoma multiforme';
cancer_name('hnsc') =           'Head and neck squamous cell carcinoma';
cancer_name('laml') =           'Acute myeloid leukemia';
cancer_name('kirc') =           'Kidney renal clear cell carcinoma';
cancer_name('kirp') =           'Kidney renal papillary cell carcinoma';
cancer_name('lgg') =            'Brain lower grade glioma';
cancer_name('lihc') =           'Liver hepatocellular carcinoma';
cancer_name('luad') =           'Lung adenocarcinoma';
cancer_name('lusc') =           'Lung squamous cell carcinoma';
cancer_name('ov') =             'Ovarian serous cystadenocarcinoma';
cancer_name('paad') =           'Pancreatic adenocarcinoma';
cancer_name('prad') =           'Prostate adenocarcinoma';
cancer_name('read') =           'Rectum adenocarcinoma';

cancer_name('stad') =           'Stomach adenocarcinoma';
cancer_name('thca') =           'Thyroid carcinoma';
cancer_name('ucec') =           'Uterine corpus endometrioid carcinoma';



