function N = add_cnvs_to_dstruct(N,cnvfile)
% ADD_CNVS_TO_DSTRUCT adds a field "iscnv" to put datastructure
%
%   N = ADD_CNVS_TO_DSTRUCT(N,CNVFILE)



%get current known CNVs
fid = fopen(cnvfile);
CNVS = textscan(fid,'%s\t%s\t%s\t');
fclose(fid);

chrn = cellfun(@(x) regexp(x,'[\d]+','match'),CNVS{1}(2:end));
chrn = cellfun(@str2num,chrn);

startpos = cellfun(@str2num,CNVS{2}(2:end));
endpos = cellfun(@str2num,CNVS{3}(2:end));

cnvstartpos = chrn*1e11 + startpos;
cnvendpos = chrn*1e11 + endpos;


% Sort cnvs by position
[cnvstartpos,si] = sort(cnvstartpos);
cnvendpos = cnvendpos(si);

cnvstartval = ones(size(cnvstartpos));
cnvendval = -ones(size(cnvendpos));
cpos = double(N.chrn)*1e11 + double(N.pos);
cposval = zeros(size(cpos));


allpos = [cpos; cnvstartpos; cnvendpos];
allval = [cposval; cnvstartval ; cnvendval];

[allpos,posidx] = sort(allpos);
allval = allval(posidx);

allval = cumsum(allval);
allval(posidx>(length(cnvstartpos)+length(cpos)))=allval(posidx>(length(cnvstartpos)+length(cpos))) + 1;


cposval = allval(posidx<=length(cpos));

m2 = cposval>=1;

N.iscnv = zeros(length(N.pos),1);
N.iscnv(m2) = 1;

    