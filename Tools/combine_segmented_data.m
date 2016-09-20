function n = combine_segmented_data(D,directory,prefix,suffix,mergedfile,segheaderlines)
% COMBINE_SEGMENTED_DATA combine the seg files in a directory into one file
% and insert sample names as given by D
%
%   N = COMBINE_SEGMENTED_DATA(D,DIRECTORY,PREFIX,SUFFIX,MERGEDFILE)
%           N is number of seg files combined
%           D is datastructure
%           DIRECTORY is the directory with the seg files to be combined
%           PREFIX is the prefix of the seg files (defaults to 'Sample')
%
cwd = pwd;

if ~exist('prefix','var') || isempty(prefix)
    prefix = 'Sample';
end

if ~exist('directory','var') || isempty(directory)
    directory = pwd;
end

if ~exist('suffix','var') || isempty(suffix)
    suffix = '.seg.dat';
end

if ~exist('segheaderlines','var') || isempty(segheaderlines)
  segheaderlines = 1;
end


if ischar(D.sdesc)
    D.sdesc = cellstr(D.sdesc);
end

if size(D.sdesc,1)>size(D.sdesc,2)
    D.sdesc = D.sdesc';
end



directory = add_slash_if_needed(directory);

cd(directory)

files = ls(['*' suffix]);
files = sort(regexp(files,['\<' prefix '[\S]+'],'match'));
n = length(files);
if ~exist('mergedfile','var')
    mergedfile = [directory 'merged_' prefix '.seg.dat'];
end
fid1 = fopen(mergedfile,'w');
fprintf(fid1,'%s\t%s\t%s\t%s\t%s\t%s\n','Sample','Chromosome','Start.bp','End.bp','NUM.SNPs','Seg.CN');
for fl = files
    
    fid = fopen(char(fl));
    towrite = textscan(fid,'%s\t%s\t%s\t%s\t%s\t%s\n');
    fclose(fid);
    nn = regexp(char(fl),'[\d]+','match');
    
    sampname = repmat(D.sdesc(str2double(nn{:})),length(towrite{1}),1);
    towrite{1} = sampname;
    towrite = reshape([towrite{:}]',1,6*length(towrite{1}));
    towrite = towrite(6*segheaderlines+1:end);
    fprintf(fid1,'%s\t%s\t%s\t%s\t%s\t%s\n',towrite{:});
    
end

fclose(fid1);

cd(cwd)
