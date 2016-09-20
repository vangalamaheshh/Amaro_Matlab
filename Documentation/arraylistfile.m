% Array List File
% 
% The array list file is a tab-delimited text file (*.txt), although
% .xls read is currently supported.  Each column has a heading as
% described below. The required information in the array list file is a
% column of array names to be included in gistic preprocessing.
% Subsequent columns of the array list file are optional and give
% information about which samples to include (include information) in
% different parts of the gistic_preprocessing algorithm as well as the location
% of the .mat file containing the D structure with the array's data.
% 
% 
% 
% The array include information is as follows:
% 
% 	  -1 in the include information indicates that the array
% 	  should be completely removed from the data structure at that
% 	  point in the analysis
% 
% 	  0 indicates that the array should not be included analyses
% 	  involving interdepence of the arrays
% 
% 	  1 indicates that the array should be used as normal, unless
% 	  a flag indicates otherwise
% 
% 	  2 indicates that use of the array in interdependent analyses
% 	  should be forced.
% 
% The column heading names and descriptions are:
% 
% 	 'array' -- required 
% 
% 		   This is a list of the names of the arrays to
% 	   	   include in the gistic analysis.  The names in the
% 	   	   array column must match the names in the **array**
% 	   	   column of the sample info file, as well as the
% 	   	   names in the '.sdesc' field of the data structure.
% 
% 	'dfile'  -- not required
% 
% 		   The location of the .mat file with the sample data,
% 		   array name, and and affy_calls.  (The location of
% 		   the output from snp_to_D.  
% 	   	   
% 	 'include in batch' -- not required
% 
% 		  -1,0,1,2 valued column giving data include
%                    information for batch effect correction.
% 
% 	 'include in merge' platforms -- not required
% 
% 		  -1,0,1,2 valued column giving data include
%                    information for merge platforms.  '0' means sample
%                    does not merge with other samples, but is carried
%                    through the analysis.  1 and 2 means the samples
%                    are merged if possible, but if there is no match
%                    the data for the other arrays receives NaNs.
% 
% 	 'include in core' -- not required
% 
% 		  -1,0,1,2 valued column giving data include
%                    information gistic core.
% 
% 	 'include in normalization' -- not required
% 
% 		  -1,0,1,2 valued column giving data include
%                    information for normalization.
% 	
