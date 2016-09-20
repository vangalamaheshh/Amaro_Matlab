function peaks = read_gistic_peaks_file(fname)
% READ_GISTIC_PEAKS_FILE read "ragged column" GISTIC peaks output into a struct
%
% PEAKS = read_gistic_peaks_file(FNAME)
%
% Read tab-delimited "amp genes" or "del genes" file and return a struct
% array PEAKS which has one element for each peak. The fields of PEAKS are:
%
%   'cytoband'  cytoband text (from row 1)
%   'qv'        q-value (row 2)
%   'resid_qv'  residual q-value (row 3)
%   'loctext'   genomic location as text (row 4)
%   'chr'       parsed chromosome part of location
%   'chrn'      numeric chromosome
%   'start'     base position of peak start
%   'end'       base position of peak end

dps = read_dlm_file(fname,char(9));
dps = vertcat(dps{:});

% get q-value row and eliminate columns where it is non-numeric
qv = str2double(dps(2,:));
goodcols = ~isnan(qv);
qv = qv(goodcols);
dps = dps(:,goodcols);

cytoband = dps(1,:);
resid_qv = str2double(dps(3,:));
loctext = dps(4,:);

parselocs = regexp(loctext,'chr(.+):([0-9]+)-([0-9]+)','tokens','once');
parselocs(cellfun(@(x) length(x)~=3,parselocs))=[];
parselocs = vertcat(parselocs{:});
pchr = parselocs(:,1)';
pchrn = chromosome2num(pchr)';
pstart = str2num(char(parselocs(:,2)))';
pend = str2num(char(parselocs(:,3)))';

% extract the ragged columns
genelistlist = cell(size(qv));
for k = 1:length(genelistlist)
    genelist = dps(5:end,k);
    genelist(cellfun(@isempty,genelist)) = [];
    genelistlist{k} = genelist;
end

peaks = struct('cytoband',  cytoband, ...
               'qv',        num2cell(qv), ...
               'resid_qv',  num2cell(resid_qv), ...
               'loctext',   loctext, ...
               'chr',       pchr, ...
               'chrn',      num2cell(pchrn), ...
               'start',     num2cell(pstart), ...
               'end',       num2cell(pend), ...
               'genes',     genelistlist );
