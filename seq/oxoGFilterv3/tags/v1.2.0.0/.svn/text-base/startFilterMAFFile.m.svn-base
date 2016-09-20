function [] = startFilterMAFFile(mafFilename, outputMAFFilename, isCheckForMatFile)
% startFilterMAFFile(mafFilename)
% 
% isCheckForMatFile -- looks for a 
%
%   Entry point for running the third incarnation of the C>A/G>T artifact.
%
%   This script will perform the filtering and produce the following:
%   	MAF Files:  Filtered and Unfiltered
%           These will be named as <outputMAFFilename> and
%           <outputMAFFilename>.reject.maf, respectively.
%       Count Files:  Only contain the values of the total number of pass
%           and rejected mutations for the input maf file.
%           <outputMAFFilename>.pass_count.txt
%           <outputMAFFilename>.reject_count.txt
%           
%       Figures
%       Tables
%
%       The latter two can be used for further reporting downstream.
%
%   Note:  This script will preserve the order of mutations with the
%       following exception:
%
%       If multiple cases are included in the input maf file and the
%       mutations of cases are interleaved.  In this case, this script will
%       produce maf files with the mutations appearing clustered by case.
%
%       Simple example: 
%           Case input (3 mutations):
%        1 ...   i1-Tumor  i1-Normal ...
%        2 ...   i2-Tumor  i2-Normal ...
%        3 ...   i1-Tumor  i1-Normal ...
%
%           Case output (3 mutations):
%        1 ...   i1-Tumor  i1-Normal ...
%        3 ...   i1-Tumor  i1-Normal ...
%        2 ...   i2-Tumor  i2-Normal ...
%
% isCheckForMatFile (default: 0) -- If value = 1, check for a mat file with name
%   mat/<mafFilename>.mat to load mafTable.
%

[~, filename, ext] = fileparts(mafFilename);
if ~exist('mat', 'dir')
   mkdir('.', 'mat');
end
matFilename = ['mat/' filename ext '.mat'];

if nargin < 3
    isCheckForMatFile = 0;
end

if ~exist(mafFilename, 'file') && (isCheckForMatFile && ~exist(matFilename, 'file'))
   error(['Given maf file does not exist: '  mafFilename])
end

% Load the maf file
if ~isCheckForMatFile
    disp(['Loading ' mafFilename ' ...'])
    [mafTable] = loadMAFTable(mafFilename);
else
    if ~exist(matFilename, 'file')
        disp(['Loading ' mafFilename ' ...'])
        [mafTable] = loadMAFTable(mafFilename);
        disp(['Saving ' matFilename ' ... '])
        save(matFilename, 'mafTable');
    end
    disp(['Loading ' matFilename ' ...'])
    load(matFilename);
end

% Retrieve all cases (one entry for each unique case, not one entry per
%   mutation)
disp(['Retrieving cases for ' mafFilename ' ...'])
[pairs] = retrieveUniqueCaseControlBarcodes(mafTable);
disp(['Retrieved ' num2str(length(pairs)) ' cases'])

% Target rate for artifacts that escape filtering
fdrThresh = .01;  

% Measured p used in the B_A(ac,p) distibution in the binomial mixture
% model.
PoxoG = .96;

% Alt Allele Counts to use
acs = [3:50];

% Null distribution (real mutation) binomial parameter, p
p0 = .5;

% Create the output directories
outputFigureDir = 'figures/';
if ~exist(outputFigureDir, 'dir')
   mkdir('.', outputFigureDir); 
end

outputTableDir = 'tables/';
if ~exist(outputTableDir, 'dir')
   mkdir('.', outputTableDir); 
end

tfid = fopen([outputTableDir 'caseTableData.tsv'] , 'w');
headers = [{'case'}; {'control'}; {'N'}; {'N_A'}; {'N_nA'};{'N_filtered'};{'N_oxo'}; {'eN_oxo'}; {'N_mut'}; {'RR'}; {'RR_CI_low'}; {'RR_CI_hi'}; {'PA'}; {'PA_CI_low'}; {'PA_CI_hi'};{'PoxoG'};{'PoxoG_Est_CI_low'};{'PoxoG_Est_CI_hi'};];
fprintf(tfid, [strjoin('\t', headers) '\n'])

totalRejectCount = 0;
totalPassCount = 0;

% Initialize a map for all of the mutations
% This map will be used to add new columns to the maf files.
mutationMap = containers.Map('KeyType', 'char', 'ValueType', 'any');

% For each unique case, run the filter
for i = 1:length(pairs)

    pair = pairs(i);    
    unameIndices = retrieveMutationsForGivenCase(pair.case, pair.control, mafTable);
    unameMafTable = trimStruct(mafTable, unameIndices);
    caseName = [pair.case ' v. ' pair.control];
    disp(['Case ' num2str(i) ': ' caseName])
    
    % Create a list of the mutations segregated by artifact mode
    %  M_A are the artifact mode mutations
    %  M_nA are the non-artifact mode mutations
    [M, M_A, M_nA] = generateMutations(mafTable, unameIndices);
        
    % Create a cut using False Discovery Rate 
    fdrStructs(i) = cutLineFDR(M, fdrThresh, PoxoG);
        
    %% Update pass and reject counts
    totalPassCount = totalPassCount + sum(~fdrStructs(i).cut);
    totalRejectCount = totalRejectCount + sum(fdrStructs(i).cut);
    
    %% Estimate the a cutline for FoxoG x alt count
    % acVal ./ acs is FoxoG.  The line is (FoxoG(i), acs(i))
    [acVal, ~, ~, ~] = estimateCutLineFromFdrOutput(fdrStructs(i), acs, PoxoG);
    cutLineFoxoGs = (acVal' ./ acs);
    
    %% Plot the artifact mode mutations
    plotArtifactMode(outputFigureDir, [pair.case ' v ' pair.control], M_A, fdrStructs(i), cutLineFoxoGs, acs)
    
    %% Plot the non-artifact mode mutations
    plotNonArtifactMode(outputFigureDir, [pair.case ' v ' pair.control], M_nA, acs)
    
    %% Estimate num oxo artifacts
    [Noxo,Nmut,alf,alfci,eNoxo] = oxogBinomialMixture(M_A.alt_read_count, round(M_A.foxog .* M_A.alt_read_count), acs, PoxoG, p0);

    %% Generate Lego plots (before and after)
    Pbefore.label = caseName;
    plotMutationSpectrumCategLego1(unameMafTable,'unit',[outputFigureDir caseName '_lego_before.png'], Pbefore);
    Pafter.label = [caseName ' Filtered'];
    plotMutationSpectrumCategLego1(trimStruct(unameMafTable, ~fdrStructs(i).cut),'unit',[outputFigureDir caseName '_lego_after.png'], Pafter)

    %TODO: Fix the plots.
    
    %% Estimate false positive (i.e. rejected real)
    % TODO: What about FP or FN numbers that clearly cannot happen?  Such
    %   as CESC-HSCX1025
    [FP,FPci,FN,FNci] =  estimateFalseRatesForOxoGfdr(M.alt_read_count, M.isArtifactMode, fdrStructs(i), Noxo, [(Noxo - eNoxo) (Noxo + eNoxo)], PoxoG);
    
    % TODO:  Do not let FN be >0 when there are no artifact mutations
    % TODO: What about RR that is higher than the number of filtered
    %   mutations?
    
    %% Estimate PoxoG error
    [PoxoGEst,PoxoGEstci]=PoxoGestimate(M.alt_read_count, round(M.foxog .* M.alt_read_count),fdrStructs(i).iART);
    
    %% Write table entry
    fprintf(tfid, [strjoin('\t', [{pair.case}; {pair.control};{num2str(length(M.foxog))}; {num2str(length(M_A.foxog))}; ...
        {num2str(length(M_nA.foxog))}; {num2str( sum(fdrStructs(i).cut))}; {num2str(Noxo)}; {num2str(eNoxo)}; {num2str(Nmut)}; ...
        {num2str(FP)}; {num2str(FPci(1))}; {num2str(FPci(2))};{num2str(FN)}; {num2str(FNci(1))}; {num2str(FNci(2))}; ...
        {num2str(PoxoGEst)}; {num2str(PoxoGEst-PoxoGEstci)}; {num2str(PoxoGEst+PoxoGEstci)} ...
        ]) '\n'] )
    
    %% Generate temporary table for later use in appending (column-wise) to output maf files
    % TODO: Finish this: pox_cutoff still needed.
    caseMutMap = addMutationEntries(pair, unameMafTable, M, fdrStructs(i), acVal, acs);
    mutationMap = [mutationMap; caseMutMap];
end
fclose(tfid);

%% Write pass and reject counts to a file.
writeFileWithSingleNumber([outputMAFFilename '.pass_count.txt'], totalPassCount);
writeFileWithSingleNumber([outputMAFFilename '.reject_count.txt'], totalRejectCount);

%% Write out new maf files (one for pass and one for reject
writeFilterMAFFiles(mafTable, outputMAFFilename, mutationMap);

disp(' ')