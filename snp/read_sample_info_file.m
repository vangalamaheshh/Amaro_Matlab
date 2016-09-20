function s=read_sample_info_file(fname,force_text,opts);
%READ_SAMPLE_INFO_FILE read sample info file and return sample information.  
%
%  S = READ_SAMPLE_INFO_FILE(FNAME,FORCE_TEXT) returns a structure array S
%  with fieldnames specified by the column headings in the sample info file 
%  FNAME.  FNAME is a tab-delimited file or an excel file (with .xls 
%  extension).  The first line of FNAME contains the tab-delimited field 
%  labels.  The subsequent lines of FNAME contain the tab-delimited field 
%  data.  FORCE_TEXT is an optional logical input that, when true, forces 
%  READ_SAMPLE_INFO_FILE to treat the info file as tab-delimited. Optional
%  input opts.  Set opts.allow_repeats to true to allow non-unique array
%  column in sample info file (default: false).
%
%
%  Updates:
%
%       18 Sept 2007:  
%               -Added conversion to char format for all .xls files read.
%               -Added row identification in error message.
%               -Changed .xls file read (removed for loops)
%
%       23 Oct 2007:
%               -Added field 'Uniq ID' ==> 'uniqid' field to allow
%               cross-platform matching by unique id.
%
%       24 Oct 2007:  
%               -Revised to allow for struct fields not defined in COLS;
%               also, file field matching can be to COLS_NAMES or COLS.
%               Jen Dobson (jdobson@broad.mit.edu)   
%       
%       21 Feb 2012L
%               -Added option, opts.allow_repeats, to allow non-unique
%               array column without generating an error.
%               Barbara Tabak (barbarat@broadinstitute.org)
%
%---
% $Id$
% $Date$
% $Last Changed By$
% $Rev$

if ~exist('opts','var')
    opts = struct;
end
opts = impose_default_value(opts,'allow_repeats', false);
%% Initialization


%Structure field names
cols={'array','name','type','cellxeno','primarymet','hormone','wga','gender',...
      'ploidy','lohctrl','dup','core','good','past_qc','rep','in100','loh','platform','batch','exp1','exp2',...
      'autopsy','age','survival','histology','alt_hist','egfr_mut','tumor_content','T','N','M','pack_years','signature','notes',...
      'gcm','p_100X','p_100H','p_250S',...
      'paired','cell_type','wga','lcm','provider','filename','nsp','sty',... %line is added for Barbara's sample info file
        'uniqid'};
    
%Info file field names
cols_names={'Array(text)','sample(text)','Type','Cell line or xenograft?','Primary/met?','Hormone status','WGA?',...
            'Gender','Ploidy(numeric)','LOH control?','duplicate of(text)','core','good quality?','Force past QC?',...
            'repeat instance of individual''s cancer?','In 100K array list?','LOH?','Platform','Batch',...
            'Expression array #1','Expression array #2','Autopsy?',...
            'Age(numeric)','survival(numeric)','Histology','Alternative Histology','EGFR mutant?','tumor content(numeric)',...
            'T','N','M',...
            'pack-years(numeric)','signature','Notes',...
            'For GCM?','Paired 100X sample','Paired 100H sample','Paired 250Sty sample',...
            'Paired','Cell type','WGA','LCM','Provider','Filename','Nsp','Sty','Unique ID'};

if ~exist('force_text','var')
  force_text=0;
end

%% Read .xls (excel) file

if strcmp(file_ext(fname),'xls') & ~force_text
  [nmr,txt,raw]=xlsread(fname);
  raw = cellfun(@num2str,raw,'UniformOutput',0);  %make sure everything's type char
  nan_pos=zeros(size(raw));
  
  nancell = cellfun(@isnan,raw,'UniformOutput',0);
  lencell = cellfun(@length,raw);
  nancell(find(lencell~=1))={0};
  
  nan_pos = cellfun(@double,nancell);

  raw(find(nan_pos))=cellstr(repmat('EMPTY',length(find(nan_pos)),1));
  f{1}=raw(1,:);
  is_xls=1;
else
  is_xls=0;
end
 
%% Read Tab-Delimited data

if ~is_xls
  [f,fid]=read_dlm_file(fname,char(9),1);
end

%% Map fields in file to standard field names

standardfields = regexprep(lower(cols_names),'\(.*\)','');
filefields = regexprep(lower(f{1}),'\(.*\)','');

%match by infofile column names
[Mi,mi1,mi2]=match_string_sets(standardfields,filefields);
%match by structure names
[Ms,ms1,ms2] = match_string_sets(cols,filefields);
%the remaining fields
matchedstructids = [mi1'  ms1'];
[mffidx,ui,uj] = unique([mi2' ms2']);
matchedfilefields = filefields(mffidx);
matchedstructfields = cols(matchedstructids(ui));

[unmatchedfields, umffidx] = setdiff(filefields,matchedfilefields);
goodidx = ~cellfun(@isempty,unmatchedfields);
unmatchedfields = unmatchedfields(goodidx);
umffidx = umffidx(goodidx);
formatedumf = regexprep(lower(unmatchedfields),'\W','');


allstructfields = [matchedstructfields formatedumf];
allffidx = [mffidx umffidx];

allfilefields = [matchedfilefields unmatchedfields];

verbose('...Matching fields to Sample Info File:',30);
verbose([repmat('      ',length(allstructfields),1) strvcat(allfilefields) repmat('  --->  ',length(allstructfields),1) strvcat(allstructfields)],30);

%remainingfilefields = setdiff(filefields,filefields(m2));


%% Make Data Structure

if ~is_xls
  form=[repmat('%s',1,length(f{1}))];
  
  %read M columns of sample info file into cells of F_dat (cells are Nx1)
  F_dat=textscan(fid,form,'bufSize',1000000,'delimiter','\t','emptyValue',NaN);

  %horizontally concatinate M cells of F_dat (new F_dat is N x M cell)
  F_dat=cat(2,F_dat{:});
  
  empty_pos=find(cellfun('isempty',F_dat));
  
  F_dat(empty_pos)=cellstr(repmat('EMPTY',length(empty_pos),1));


else
  F_dat=raw(2:end,:);
end

if ~strcmp(allstructfields{1},'array')
  allstructfields=[ 'array' allstructfields ];
  allffidx=[ 1 allffidx ];
  verbose('Cannot find an "Array" column ... using first column as array',30);
end

%s is structure whos fields are the column names of F_dat
s=cell2struct(F_dat(:,allffidx),allstructfields,2);

  
%% Error catching: Make sure all sample names are unique
if ~opts.allow_repeats
    [UList,ULidx,SLidx] = unique({s.array});
    
    try ULidx == SLidx;
    catch
        N = hist(ULidx,max(ULidx));
        notunique = find(N~=1);
        error(['The array column of the sample info file is not unique.'...
            '  Please see rows: ' num2str(notunique+1)])
    end
end
