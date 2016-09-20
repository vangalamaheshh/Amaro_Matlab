function gp_snp_to_D(varargin)
% GP_SNP_TO_D gene pattern wrapper function to convert snp datafile to
% matlab structure.  
%
%      gp_snp_to_D(-i input_filename -o output_file_name -a has_allele_data
%      -nl number_of_lines -ec extra_column -h header_rows -sz read_chunks)
%
%      gp_snp_to_D reads intensity data from a SNP data file and writes the 
%      the data to a matlab structure saved in the output file.
%
%      The options passed to GP_SNP_TO_D are:
%                -i: (string) input file name (REQUIRED)
%                -o: (string) output file name (REQUIRED)
%                -a: (int) allele data (0 if no allele data; 1 if allele data
%                as [intensity 1; call; intensity 2], 2 if allele data given
%                as [intensity 1; intensity 2; call]
%                -nl: (int) number of lines in file
%                -ec: (int) number of extra columns in file between
%                chromosome position column and intensity data
%                 -h: number of header lines
%                -sz: number of lines of data to read at each iteration of
%                TEXTSCAN
%
%      Example: gp_snp_to_D('-i','testsnps','-o','mstructout.mat','-a','0','-h','0')
%
%      For file format specifications or additional information, see READ_MODELLED_DATA_FILE.
%
%      FUTURE: add skip_snp
%
% ---
% $Id: gp_snp_to_D.m 74 2007-09-15 12:02:38Z gadgetz $
% $Date: 2007-09-15 08:02:38 -0400 (Sat, 15 Sep 2007) $
% $LastChangedBy: gadgetz $
% $Rev: 74 $

verbose([datestr(now) ': Running GP_SNP_TO_D.'],10)


%addpath /xchip/gistic/Code/GISTIC1.0/matlab/
%addpath /xchip/gistic/Code/GISTIC1.0/matlab/gp_modules
%addpath /xchip/gistic/Code/GISTIC1.0/matlab/snp

addpath ~/GISTIC/Code/gp_modules
addpath ~/GISTIC/Code
addpath ~/GISTIC/Code/snp

%Populate the information stucture
a=handle_args({'i','o','a','nl','ec','h','sz'},varargin);



if isempty(a.i) || isempty(a.o)
  error('must have input and output files');    %Check that input and output files are specified
else
  infile=a.i;                                   %Assign filenames to variables INFILE and OUTFILE from information structure a
  outfile=a.o;
end
if ~isempty(a.a) 
  has_allele_data=str2num(a.a);                 %Assign whether allele data is provided
else
  has_allele_data=0;
end

if ~isempty(a.nl)                               %Assign number of lines to NLINES from A
  nlines=str2num(a.nl);
  if isempty(nlines)
    error('number of lines is not a number');
  end
else
  nlines=-1;
end

if ~isempty(a.ec)                               %Assign number of extra columns to EXTRA_COLUMNS
  extra_columns=str2num(a.ec);              
  if isempty(extra_columns) || extra_columns<0
    error('extra columns should be a non-negative number');
  end
else
  extra_columns=0;
end

if ~isempty(a.h) 
  header_rows=str2num(a.h);                     %Assign number of header rows
  if isempty(header_rows) || header_rows<0
    error('header rows should be a non-negative number');
  end
else
  header_rows=0;
end

if ~isempty(a.sz) 
  read_chunks=str2num(a.sz);                    %Assign read chunks
  if isempty(read_chunks) || read_chunks<0
    error('read chunks should be a non-negative number');
  end
else
  read_chunks=10000;
end


verbose([datestr(now) ': Passing ' infile ' to SNP_TO_D.M. (GP_SNP_TO_D.M)'],10);

M=snp_to_D(infile,nlines,has_allele_data,extra_columns,header_rows,read_chunks);

save(outfile,'M','-v7.3');
verbose([datestr(now) ': Data structure saved to' outfile '(GP_SNP_TO_D.M)'],10)
