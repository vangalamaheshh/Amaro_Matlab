function categ2wig(wigfile,categdir,outfile)
% given an example wiggle file,
% reformats a set of per-chromosome category files in categdir
% into a wiggle file with the same intervals etc.
%
% procedure:
%   (1) read through wiggle file
%       -- save all "track" and "fixedStep" lines
%       -- give error if any "variableStep" lines are encountered: these are not supported yet
%       -- for each fixedStep line, save the number of positions listed after it
%   (2) allocate memory to store the category of every position in the wiggle file
%   (3) load each chromosome's category file, extract the positions
%   (4) write the category wiggle file
%
% Mike Lawrence 2009-11-17

[a b] = system(['grep -n ''^\D'' ' wigfile]);
b = split(b,char(10));
b = setdiff(b,{''});

track = grep('^\d*:track',b,1);
fs = grep('^\d*:fixedStep',b,1);
vs = grep('^\d*:variableStep',b,1);
other = setdiff(1:length(b),[track;fs;vs]);

if ~isempty(vs), error('variableStep not supported'); end
if ~isempty(other), fprintf('Invalid lines:\n');disp(b(other));error('Aborting.');end

c = parse(b(fs),'^(\d*):fixedStep chrom=chr(\S*) start=(\d*) step=(\d*)$',{'line','chr','start','step'});


keyboard


%%%%%%%%%%%%%%  not finished

%%%% obsolete: replaced by
%
%   /xchip/tcga_scratch/lawrence/db/allcateg/txt2wig.cmd

