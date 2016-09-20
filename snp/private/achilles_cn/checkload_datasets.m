% TODO - fix the description below, which is totally wrong
% If achilles and gistic CLE datasets are not in memory, load them.
% The cell lines are reduced to just those that match.
% After invocation: 
%      G.sdesc  Names of samples in CLE with copy number data
%      rg    struct array containing the reference genome
%      cyto  struct array of cytoband information
%      acx   maps indexes matching AC.dat samples to correspond
%            with CL.dat
%      
%!        etc...

% manual cell line name corrections
clmap_file = 'CLE-Achilles_CLmatch_101012.tab';

%% map cell lines between Achilles and CLE datasets

verbose('=> Mapping cell lines.',10);

% read cell line to array name mapping information

% read as SIS
%! (All we need from here is cline, all_cline to determine mismatch cause,
%! everything else is done in process_cle)

%!SI = read_sample_info_file(sampinfo_file);

all_cline = {SI.celllineprimaryname};
[~,matches] = match_string_sets_hash({SI.array},G.sdesc);
cline = all_cline(matches);
%!!!get_cleCL.sis = SI(matches);
cline = regexprep(lower(cline),'[ :-]+','');


acline = regexprep(lower(AC.sampleID),'[ :-]+','');
if exist('clmap_file','var') && ~isempty(clmap_file)
    % use a custom translation file to fix any trivial cell line name mismatches
    % (column 1 is gistic/CLE name; column 2 is hairpin header name)
    fid = fopen([input_dir clmap_file],'r');
    textscan(fid,'%s%s',1,'Delimiter',char(9));
    clmap = textscan(fid,'%s%s','Delimiter',char(9));
    fclose(fid);
    xlat = regexprep(lower(clmap{2}),'[ :-]+','');
    [~,mA,mB] = match_string_sets_hash(acline,xlat);
    acline(mA) = regexprep(lower(clmap{1}(mB)),'[ :-]+','');
end

% use the translated cell line names to map cell lines
[M,gcx,acx] = match_string_sets(cline,acline);
verbose('...matched %d of %d cell lines',10,...
        length(gcx),length(AC.sampleID));

% create & write list of unmatched cell lines
file = [output_dir DATA_ID '.unmatched_cell_lines.tab'];
unmatched_clines = AC.sampleID(~any(M,2))';
reason = cellstr(repmat('no CLE match',length(unmatched_clines),1));
[~,nodata] = match_string_sets(unmatched_clines,all_cline);
reason(nodata) = {'no CN data'};
write_filtered_tabcols(file,[],{'Achilles name',unmatched_clines},{'reason',reason});

%% reduce copy number data samples to those that match experimental samples
