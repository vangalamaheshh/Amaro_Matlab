function [] = startFilterMAFFile(mafFilename, outputMAFFilename, outputDir, isCheckForMatFile, isGeneratingPlots, globalPoxoG, artifactThresholdRate)
% [] = startFilterMAFFile(mafFilename, outputMAFFilename, outputDir, isCheckForMatFile, globalPoxoG, artifactThresholdRate, isGeneratingPlots)
% 
% mafFilename -- input maf filename to be filtered.  This file will not be
%   modified.
% 
% outputMAFFilename -- filtered maf filename to be generated.
%
% outputDir (default: './') -- outputDir for filter results that are 
%   used for the report.
%   Note:  outputMAFFilename, figures, and tables will be located 
%       in outputDir.  The mat files will not be located in the outputDir.
% 
% globalPoxoG (default: .96) -- poxoG value to use in all estimations of
%   the given mafFile.  
%
% artifactThresholdRate (default .01) -- Expected proportion of mutations
%   that are artifacts.  e.g. .01 --> 1%  
%
% isCheckForMatFile (default: 0) -- looks for a mat file that contains the
%   mafTable already.  This is just to save on load time.  If value = 1, 
%   check for a mat file with name  mat/<mafFilename>.mat to load mafTable.
% 
% isGeneratingPlots (default: 0) -- true/false for whether this script
%   should generate the orchestra and lego png files.  These will be put in
%   the [outputDir]/figures/ .
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
% See also load_table

[~, filename, ext] = fileparts(mafFilename);
if ~exist('mat', 'dir')
   mkdir('.', 'mat');
end
matFilename = ['mat/' filename ext '.mat'];

%% Parse arguments and determine which should use default values.
if nargin < 7 || ~exist('artifactThresholdRate') || isempty(artifactThresholdRate)
    artifactThresholdRate = .01;
end

if nargin < 6 || ~exist('globalPoxoG') || isempty(globalPoxoG)
    globalPoxoG = .96;
end

if nargin < 5 || ~exist('isGeneratingPlots') || isempty(isGeneratingPlots)
    isGeneratingPlots = 0;
end

if nargin < 4 || ~exist('isCheckForMatFile') || isempty(isCheckForMatFile)
    isCheckForMatFile = 0;
end

if nargin < 3 || ~exist('outputDir') || isempty(outputDir)
    outputDir = './';
end

if ~exist(mafFilename, 'file') && (isCheckForMatFile && ~exist(matFilename, 'file'))
   error(['Given maf file does not exist: '  mafFilename])
end

%% Load the maf file
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

%% Retrieve all cases (one entry for each unique case, not one entry per
%   mutation)
disp(['Retrieving cases for ' mafFilename ' ...'])
[pairs] = retrieveUniqueCaseControlBarcodes(mafTable);
disp(['Retrieved ' num2str(length(pairs)) ' cases'])

% Target rate for artifacts that escape filtering
fdrThresh = artifactThresholdRate;  

% Measured p used in the B_A(ac,p) distibution in the binomial mixture
% model.
PoxoG = globalPoxoG;

% Alt Allele Counts to use
acs = [3:50];

% Null distribution (real mutation) binomial parameter, p
p0 = .5;

%% Create the output directories
if ~exist(outputDir, 'dir')
   mkdir(outputDir); 
end

outputFigureDir =[ outputDir '/figures/'];
if ~exist(outputFigureDir, 'dir')
   mkdir('.', outputFigureDir); 
end

outputTableDir = [outputDir '/'];
if ~exist(outputTableDir, 'dir')
   mkdir('.', outputTableDir); 
end

%% Initialize a file to list the cases in the input maf
caseFid = fopen([outputDir '/allCases.txt'], 'w');

%% Initialize a table for summary data where each case is a row in a table.
tfid = fopen([outputTableDir 'caseTableData.tsv'] , 'w');
headers = [{'case'}; {'control'}; {'N'}; {'N_A'}; {'N_nA'};{'N_filtered'};{'N_oxo'}; {'N_oxo_CI_low'}; {'N_oxo_CI_high'};{'N_mut'}; {'RR'}; {'RR_CI_low'}; {'RR_CI_hi'}; {'PA'}; {'PA_CI_low'}; {'PA_CI_hi'};{'PoxoG'};{'PoxoG_Est_CI_low'};{'PoxoG_Est_CI_hi'};];
fprintf(tfid, [strjoin('\t', headers) '\n'])

totalRejectCount = 0;
totalPassCount = 0;

%% Initialize a map for all of the mutations
% This map will be used to add new columns to the maf files.
mutationMap = containers.Map('KeyType', 'char', 'ValueType', 'any');

%% Initialize an array for storage of Passing Artifact & Rejected Real counts (and hi/low CIs) in a summary view
allPAs = zeros(length(pairs),3);
allRRs = zeros(length(pairs),3);
Ns = zeros(length(pairs),1);
%% For each unique case, run the filter
for i = 1:length(pairs)

    pair = pairs(i);    
    unameIndices = retrieveMutationsForGivenCase(pair.case, pair.control, mafTable);
    unameMafTable = trimStruct(mafTable, unameIndices);
    
    caseName = [pair.case ' v. ' pair.control];
    fprintf(caseFid, [caseName '\n']);
    
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

    %% Estimate num oxo artifacts
    [Noxo,Nmut,alf,alfci,NoxoCI] = oxogBinomialMixture(M_A.alt_read_count, round(M_A.foxog .* M_A.alt_read_count), acs, PoxoG, p0);

    %% If we are 95% certain that this case has less than 1% artifacts, do not filter.  Do this by setting 
    if NoxoCI(2) < (fdrThresh*length(M.alt_read_count))
       disp('We are 95% certain that this case has less than 1% artifacts.  Doing no filtering.')
       fdrStructs(i).cut = zeros(size(fdrStructs(i).cut));
    end
    
    %% Estimate the a cutline for FoxoG x alt count
    % acVal ./ acs is FoxoG.  The line is (FoxoG(i), acs(i))
    [acVal] = estimateCutLineFromFdrOutput(fdrStructs(i), acs, PoxoG);
    cutLineFoxoGs = (acVal' ./ acs);

    %% Plotting figures
    if isGeneratingPlots
        
        %% Plot the artifact mode mutations
        plotArtifactMode(outputFigureDir, caseName, M_A, fdrStructs(i), cutLineFoxoGs, acs)

        %% Plot the non-artifact mode mutations
        plotNonArtifactMode(outputFigureDir, caseName, M_nA, acs)

        %% Generate Lego plots (before and after)
        Pbefore.label = caseName;
        plotMutationSpectrumCategLego1(unameMafTable,'unit',[outputFigureDir caseName '_lego_before'], Pbefore);
        Pafter.label = [caseName ' Filtered'];
        plotMutationSpectrumCategLego1(trimStruct(unameMafTable, ~fdrStructs(i).cut),'unit',[outputFigureDir caseName '_lego_after'], Pafter)
    
    else
        disp('Figures are not being generated.')
    end
    
    %% Estimate false positive (i.e. rejected real)
    [FP,FPci,FN,FNci] = estimateFalseRatesForOxoGfdr(M.alt_read_count, M.isArtifactMode, fdrStructs(i), Noxo, NoxoCI, PoxoG);
    allPAs(i,:) = [FN FNci(1) FNci(2)]; 
    allRRs(i,:) = [FP FPci(1) FPci(2)]; 
    Ns(i) = length(M.alt_read_count);
    
    %% Estimate PoxoG error.  
    % Note that this estimate is not actually used.  This estimate is
    %  recorded for later analysis.
    [PoxoGEst,PoxoGEstci] = PoxoGestimate(M.alt_read_count, round(M.foxog .* M.alt_read_count),fdrStructs(i).iART, [], Noxo, NoxoCI);
    
    %% Write stdout
    tbl_entry = [strjoin('\t', [{pair.case}; {pair.control};{num2str(length(M.foxog))}; {num2str(length(M_A.foxog))}; ...
        {num2str(length(M_nA.foxog))}; {num2str( sum(fdrStructs(i).cut))}; {num2str(Noxo)}; {num2str(NoxoCI(1))}; {num2str(NoxoCI(2))}; {num2str(Nmut)}; ...
        {num2str(FP)}; {num2str(FPci(1))}; {num2str(FPci(2))};{num2str(FN)}; {num2str(FNci(1))}; {num2str(FNci(2))}; ...
        {num2str(PoxoGEst)}; {num2str(PoxoGEstci(1))}; {num2str(PoxoGEstci(2))} ...
        ]) '\n'];
    fprintf(1,  tbl_entry)
    
    %% Write table entry to table file
    fprintf(tfid, tbl_entry )
    
    %% Generate temporary table for later use in appending (column-wise) to output maf files
    caseMutMap = addMutationEntries(pair, unameMafTable, M, fdrStructs(i), acVal, acs, PoxoG);
    mutationMap = [mutationMap; caseMutMap];
end
fclose(tfid);
fclose(caseFid);

%% Write pass and reject counts to a file.
writeFileWithSingleNumber([outputDir '/' outputMAFFilename '.pass_count.txt'], totalPassCount);
writeFileWithSingleNumber([outputDir '/' outputMAFFilename '.reject_count.txt'], totalRejectCount);

%% Write out new maf files (one for pass and one for reject
disp('Writing MAF files...')
writeFilterMAFFiles(mafTable, [outputDir '/' outputMAFFilename], mutationMap);

%% Plot summary data across all cases
if isGeneratingPlots
    disp(' Generating summary plots for all cases...')
    [dummy, basename, ext] = fileparts(outputMAFFilename);
    plotFalseRates([outputFigureDir '/' outputMAFFilename '.RR'], [basename ' Rejected Reals'], 'Rejected Real Count', Ns, allRRs);
    plotFalseRates([outputFigureDir '/' outputMAFFilename '.PA'], [basename ' Passing Artifacts'], 'Passing Artifact Count', Ns, allPAs);
    disp('False rate plots completed.')
    
    % Plot summary orchestra plot of all artifact mode mutations
    allFdrStructs = fdrStructs(1);
    for k = 2:length(fdrStructs)
        allFdrStructs = mergeStruct(allFdrStructs, fdrStructs(k));
    end
    [M_summary, M_A_summary, M_nA_summary] = generateMutations(mafTable, [1:length(mafTable.Start_position)]);
    plotArtifactMode(outputFigureDir, [outputMAFFilename ' All Artifact Mode'], M_A_summary, allFdrStructs , [], acs);
    plotNonArtifactMode(outputFigureDir, [outputMAFFilename ' All Non-Artifact Mode'], M_nA_summary, acs)
    disp('Summary Orchestra plots completed.')
end

%% Clean up
close all    % close all figures so xvnc can terminate
disp('Done.')