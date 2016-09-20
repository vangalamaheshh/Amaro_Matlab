function M=read_modelled_data_file(fname,nlines,affy_call_names,has_chr_pos,has_allele_data,just_calls,extra_cols,header_rows,sz)
%READ_MODELLED_DATA_FILE convert snp data file into matlab
%structure.
%   M = READ_MODELLED_DATA_FILE
%   (FNAME,NLINES,AFFY_CALL_NAMES,HAS_CHR_POS,HAS_ALLELE_DATA,JUST_CALLS,EXTR
%   A_COLS,HEADER_ROWS,SZ)
%
%   The input parameters are:
%
%       FNAME -- the name of the input file.  FNAME is a tab-delimited file
%       containing SNP array data.  After the header rows, it contains one
%       row for each SNP.  After prefix columns, the columns are the intensity values for each
%       sample.  The column format is: :
%
%
%          |<------------------------Prefix Columns-------------------->|<--Intensity data--->
%          SNP ID | CHROMOSOME(optional)| CHROMOSOME POSITION(optional) |      Columns
%
%          Intensity data columns are multiple columns as given below.  The INTENSITY DATA columns are
%          repeated n times, where n is the number of samples.
%
%          If (JUST_CALLS == 0) and (HAS_ALLELE_DATA == 0), INTENSITY DATA has the format:
%                     INTENSITY     
%
%          If (JUST_CALLS == 1) and (HAS_ALLELE_DATA ==0), INTENSITY DATA has the format:
%                     INTENSITY | CALL  
%
%          If HAS_ALLEEL_DATA = 1, INTENSITY DATA has the format:
%                     ALLELE 1 INTENSITY | CALL | ALLELE 2 INTENSITY 
%
%          If HAS_ALLELE_DATA = 2, INTENSITY DATA has the format:
%                     ALLELE 1 INTENSITY | ALLELE 2 INTENSITY | CALL
%
%       NLINES -- total number of lines in file.  When passed in, NLINES should be 
%       equivalent to the number of SNP locations PLUS the number of 
%       header lines. 
%
%       AFFY_CALL_NAMES -- Integer (-1,0, or 1) value indicates format of
%       CALLS in input file.  -1 signifies no calls made.  0 signifies that
%       the call made is 'A', 'B', 'AB', or 'NoCall'.  1 signifies that the call
%       format is 'AA', 'BB', or 'AB'.
%      
%       HAS_CHR_POS --  Integer (1 or 0) value indicates whether each SNPs
%       chromosome and chomosome position is given in input file.
%       Chromosome can be given if and only if chromosome position is
%       given.  (For gene pattern SNP files, HAS_CHR_POS = 1.)
%
%       HAS_ALLELE_DATA --  Integer (0,1, or 2) value indicates whether data
%       file gives intensities for each allele and indicates format.  See
%       FNAME description above.
%
%       JUST_CALLS --  Integer (1 or 0) value indicates whether data file
%       gives calls only (no intensity data) for each snp.  (For gene pattern SNP files,
%       JUST_CALLS = 0.)
%
%       EXTRA_COLS -- Integer value indicating the number of extra columns
%       between the chromosome position column and the first intensity value
%       in the data file.  (EXTRA_COLS can also be a character string giving
%       the format of all columns preceeding the first column of intensity
%       values.)
%
%       HEADER_ROWS -- (int) number of header rows in input file FNAME
%
%       SZ -- Number of lines of data file to read during each iteration (default 1000).
%
%   M is the matlab structure whose fields include:
%
%       M.DAT -- an m-by-n array of snp intensity data, where m is the
%       number of lines of data (number of SNPS) and n is the number of samples.
%
%       M.ADAT -- an m-by-n-by-p array of intensity data for each allele, where m is the number
%       of lines, n is the number of samples, and p is 2.  M.ADAT(:,:,1)
%       gives data for first allele, M.ADAT(:,:,2) gives data for second
%       allele. M.ADAT only exists if allele data are given in input file.
%
%       M.AFFY_CALLS -- an m-by-n array of integers indicating the call, 
%       where m and n are as given above.  1 indicates a call of 'AA', 2
%       indicates a call of 'BB', and 3 indicates a call of 'AB'.  NaN is
%       used to indicate NoCall.  (M.AFFY_CALLS exists only if call data are
%       given in input file.)
%
%       M.MARKER -- an m-by-1 cell array, m as above.  M.MARKER lists the
%       SNP marker names, taken from the first column of FNAME.
%
%       M.POS -- a column vector of length m, m as above.  Gives chromosome
%       positions for each SNP.  (Does not exist if HAS_CHR_POS = 0.)
%
%       M.CHR -- an m-by-1 cell array of string, m as above.  Gives chromosome name as a string
%       for each SNP.  (Does not exist if HAS_CHR_POS = 0.)
%
%       M.CHRN -- an m-by-1 cell array, m as above.  Gives chromosome number (int)
%       for each SNP.  (X = 23; Y = 24)  (Does not exist if HAS_CHR_POS =0.)
%
%       M.SDESC -- Cell array of strings giving column labels from data file.
%
%---
% $Id$
% $Date$
% $LastChangedBy$
% $Rev$



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% INITIALIZE PARAMETERS %%%%%%%%%%%%%%%%%%
%This section of code initializes the parameters needed for 
%the file read operation.  

if ~exist('header_rows','var')
  header_rows=0;
end
if ~exist('extra_cols','var')
  extra_cols=0;
end

%Get or define number of lines to read per iteration
if ~exist('sz','var') || isempty(sz)
  sz=1000;
end
verbose(['... Read ' num2str(sz) ' lines per iteration'],20);


% Get number of lines in file
if ~exist('nlines','var') || isempty(nlines) || nlines<0
  verbose('... Counting lines',20);
  nlines=line_count(fname);
  verbose(['... nlines =' num2str(nlines)],20);
  nlines=nlines-1;
  nlines=nlines-header_rows;
else
  nlines=nlines-1-header_rows;
end

if ~exist('affy_call_names','var')
  affy_call_names=0;
end

if ~exist('has_chr_pos','var')
  has_chr_pos=0;
end

if ~exist('has_allele_data','var')
  has_allele_data=0;
end

% number of columns to step over for next sample in input file array
step_size=1+double(affy_call_names>=0)+(has_allele_data>0);  


if ~exist('just_calls','var')
  just_calls=0;
end
if just_calls
  step_size=1;
  if (has_allele_data>0)
    error('cannot have allele data if just calls');
  end
end


% f{end} is a cell array of strings giving the labels of the data columns
% read_dlm_file leaves the file position indicator pointing to the first
% line of data
[f,fid]=read_dlm_file(fname,char(9),1+header_rows);

if isempty(f{end}{end})
  extra_tab=1;
  extra_tab_st='%*s';
else
  extra_tab=0;
  extra_tab_st='';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%  END INITIALIZE PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  SET FORMAT STRING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This section selects the correct formating string to use with TEXTSCAN, 
%according to the information given in the input file.
%Each line of FNAME is read using the format specified in FORM.
%(See matlab built-in function TEXTSCAN for more information.)
%EXTRA_FORM specifies prefix. 


if ~ischar(extra_cols)
  if has_chr_pos
    if extra_cols>0
      extra_form=[ '%s%s%s' repmat('%*s',1,extra_cols)]; %(SNP)(CHR)(CHR_POS)(EXTRA_COLS)
    else
      extra_form=[ '%s%s%s' ]; %(SNP)(CHR)(CHR_POS)
    end
  else
    if extra_cols>0
      extra_form=[ '%s' repmat('%*s',1,extra_cols)]; %(SNP)(EXTRA_COLS)
    else
      extra_form=[ '%s' ]; %(SNP)
    end
  end    
else
  extra_form=extra_cols;
  extra_cols=length(find(extra_cols=='%'))-1-2*has_chr_pos;
end

if has_allele_data
  if has_allele_data==2
    form=[ extra_form repmat('%f%f%s',1,floor((length(f{end})-extra_tab-1-2*has_chr_pos-extra_cols)/step_size)) ...
           extra_tab_st];
  else    
    form=[ extra_form repmat('%f%s%f',1,floor((length(f{end})-extra_tab-1-2*has_chr_pos-extra_cols)/step_size)) ...
           extra_tab_st];
  end    
else
  if affy_call_names>=0
    form=[ extra_form repmat('%f%s',1,floor((length(f{end})-extra_tab-1-2*has_chr_pos-extra_cols)/step_size)) extra_tab_st];
  else    
    form=[ extra_form repmat('%f',1,floor((length(f{end})-extra_tab-1-2*has_chr_pos-extra_cols)/step_size)) extra_tab_st];
  end
end
if just_calls
  form=[ extra_form repmat('%s',1,floor((length(f{end})-extra_tab-1-2*has_chr_pos-extra_cols)/step_size)) extra_tab_st];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END SET FORMAT STRING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% READ DATA FROM FILE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This section reads the data from FNAME into the cell array M_dat.  (NOT M.DAT!!)
%File is read using the matlab built-in TEXTSCAN.  File is read NEXT_SZ 
%lines at a time.  The nth cell of M_dat (M{n}) corresponds to the nth
%column in the data file.
%and deals the data to the 
%appropriate structures of M.

keep_reading=1;
M_dat_full={};
verbose([datestr(now) '-- Reading modelled data file (read_modelled_data_file.m)'],10);
tot=0;

%number of samples
ns=floor((length(f{end})-extra_tab-1-has_chr_pos*2-extra_cols)/step_size);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialize M structure %%%
M.dat=zeros(nlines,ns);
if has_allele_data
  M.adat=zeros(nlines,ns,2);
end
cur_calls=zeros(sz,ns);
if affy_call_names>=0
  M.affy_calls=zeros(nlines,ns);
end
M.marker=cell(nlines,1);
if has_chr_pos
  M.pos=zeros(nlines,1);
  M.chr=cell(nlines,1);
end
%%%% END INITIALIZE %%%%%%%%

%Conversion tables for CALL data
conv_tab=[1 3 0 0 0 NaN];
conv_tab_affy=[0 0 0 0 0 NaN];

prev=1;
next_sz=min(sz,nlines-tot+1);

%%%% DATA READ LOOP %%%%%%%%
%%%% read nlines at a time %
% NOTE: M_dat is different from M.dat!!

while (keep_reading)
  M_dat=textscan(fid,form,next_sz,'bufSize',50000000,'delimiter','\t','emptyValue',NaN,'treatAsEmpty','NA');
 
  if (length(M_dat{1})<next_sz)  %length(M_dat{1}) is number of lines read
    keep_reading=0;     %if there exist fewer lines than were read, stop reading
  end

  %if more lines were read than there are data, get rid of excess lines
  tot=tot+length(M_dat{1});
  if tot>=nlines
    if tot>nlines
      for kk=1:length(M_dat)     %length(M_dat) is number of columns read
        M_dat{kk}=M_dat{kk}(1:(end-(tot-nlines)));  
      end    
      tot=nlines;
    end
    keep_reading=0;
  end
  
 
  
  verbose(['...' num2str(tot) ' lines read'],20);

  %%%% DEAL DATA to M STRUCTURE %%%%
  
  M.marker(prev:tot)=M_dat{1};  
  
  if has_chr_pos
    M.chr(prev:tot)=M_dat{2};
    pos_w_nan=M_dat{3};
    idx_w_nan=find(cellfun('isempty',pos_w_nan));
    pos_w_nan(idx_w_nan)=cellstr(repmat('NaN',length(idx_w_nan),1));
    M.pos(prev:tot)=str2num(strvcat(pos_w_nan));
    start_col=4;
  else
    start_col=2;
  end
  
  if ~just_calls
    if has_allele_data
      if has_allele_data==2
        M.adat(prev:tot,:,1)=cat(2,M_dat{start_col:3:end});
        M.adat(prev:tot,:,2)=cat(2,M_dat{(start_col+1):3:end});
        M.dat(prev:tot,:)=sum(M.adat(prev:tot,:,:),3);
      else
        M.adat(prev:tot,:,1)=cat(2,M_dat{start_col:3:end});
        M.adat(prev:tot,:,2)=cat(2,M_dat{(start_col+2):3:end});
        M.dat(prev:tot,:)=sum(M.adat(prev:tot,:,:),3);
      end
    else
      M.dat(prev:tot,:)=cat(2,M_dat{start_col:step_size:end});
    end
  end
  

 %Read CALL data
  if affy_call_names>=0
    for j=(start_col+1+(has_allele_data==2)-just_calls):step_size:length(M_dat) %iterate over each sample
      xx=strvcat(M_dat{j}); % the calls for this sample (vertical array)
      xxl=cellfun('length',M_dat{j});
      if affy_call_names>0
        if size(xx,2)>1 
	  cur_calls(1:size(xx,1),(j-1-start_col+just_calls+step_size-(has_allele_data==2))/step_size) = conv_tab_affy(xxl)'+(xx(:,1)=='A')+(xx(:,2)=='B')*2;
	else 
	  cur_calls(1:size(xx,1),(j-1-start_col+just_calls+step_size-(has_allele_data==2))/step_size) = conv_tab_affy(xxl)'+(xx(:,1)=='A');
        end
      else
	cur_calls(1:size(xx,1),(j+1-start_col+just_calls-(has_allele_data==2))/step_size) = conv_tab(xxl)'+(xx(:,1)=='B');
      end
      %    cur_calls(:,(j-1)/2)=str2num(strvcat(regexprep(M_dat{j},...
      %        {'^NoCall$','^A$','^B$','^AB$'},{'NaN','1','2','3'})));
    end
    M.affy_calls(prev:tot,:)=cur_calls(1:size(xx,1),:);
  end
  prev=tot+1;
end 



if has_allele_data
  M.sdesc=f{end}((start_col+extra_cols):step_size:(end-extra_tab));
else
  M.sdesc=f{end}((start_col+extra_cols):step_size:(end-extra_tab));
end

if has_chr_pos
  M=add_chrn(M);  %adds field '.chrn' to M structure
end

fclose(fid);

if (0) % use impute_missing_values to handle NaNs
  nans=isnan(M.dat);
  if nnz(nans)>0
    bad_lines=find(any(nans,2));
    verbose([datestr(now) 'Replacing NaNs with 0 in lines:' num2str(bad_lines)],10);
    M.dat(find(nans))=0;
  end
end

verbose([datestr(now) ':Data Structure Created. (read_modelled_data_file.m)'],10)

