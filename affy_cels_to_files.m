function affy_cels_to_files(flist)

addpath ~/matlab/GenePattern
if ~exist('my_gp','var')
  global my_gp;
  my_gp=[];
  my_gp=GenePatternServer('http://genepattern:8080','gadgetz@broad.mit.edu');
end

