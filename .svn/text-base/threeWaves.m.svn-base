function t=threeWaves(q,iter)

tic;
basicDir='/xchip/data/gadgetz/new_partial/packSS/';

addpath([basicDir,'main'],...
[basicDir,'suppMatlabFiles/Gbp2DimPackage'],...
[basicDir,'suppMatlabFiles/IMAGES/'],...
[basicDir,'suppMatlabFiles/GENERAL_PURPOSE_M_FILES/'],...
[basicDir,'suppMatlabFiles/matlabSS/boykovMatlab/'],...
[basicDir,'suppMatlabFiles/matlabSS/BP/'],...
[basicDir,'suppMatlabFiles/matlabSS/joachimsMatlab'],...
[basicDir,'suppMatlabFiles/matlabSS/WBP'],...
[basicDir,'suppMatlabFiles/matlabSS/'],...
[basicDir,'suppMatlabFiles/stand_alone/'],...
[basicDir,'suppMatlabFiles/GENERAL_PURPOSE_M_FILES/spectral_clustering_toolbox']        )
[basicDir,'suppMatlabFiles/'];
% q=2
if q==2
  inputInfo=struct;
  inputInfo.saveName=['tt2-' num2str(iter)];
  inputInfo.saveDir=[basicDir,'resTmpSemi/'];
  inputInfo.qpotts=2;
  inputInfo.backboneFile=[basicDir,'backboneData/','backbone_st2'];
  inputInfo.dir=struct;
  inputInfo.dir.basicDir=basicDir;
  inputInfo.dir.imageDir=[basicDir,'imageDir/'];
  inputInfo.dir.nameData='handICML';
  inputInfo.numRep=5; % was 100
  inputInfo.deltaT=0.02; % was 0.05
  runAllMethods2(inputInfo,1);
end
% q=3
if q==3
  inputInfo=struct;
  inputInfo.saveName=['tt3-' num2str(iter)];
  inputInfo.saveDir=[basicDir,'resTmpSemi/'];
  inputInfo.qpotts=3;
  inputInfo.backboneFile=[basicDir,'backboneData/','backbone_st3'];
  inputInfo.dir=struct;
  inputInfo.dir.basicDir=basicDir;
  inputInfo.dir.imageDir=[basicDir,'imageDir/'];
  inputInfo.dir.nameData='handICML';
  inputInfo.numRep=5; % was 100
  inputInfo.deltaT=0.02; % was 0.05
  runAllMethods2(inputInfo,1);
end

t=toc;













