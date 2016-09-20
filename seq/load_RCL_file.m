function [P C] = load_RCL_file(filename)
% load_RCL_file(filename)
%
% Loads a file that contains the results of RegionCovPerLane.java
%            (no header line)
%            <targname/genename> <chr> <start> <end> <lane0_count> <lane1_count> ... <laneN_count>
%
% Returns:
%     P, a struct with name, chr, start, end
%     C, a matrix with the lanecounts
%
%     Both P and C are sorted by (chr,start)
%
% If collapse_flag is TRUE (default), collapses 
%
% Mike Lawrence 2009-10-16

if ~exist(filename,'file'), error('Not found: %s',filename); end

numcols = get_colcount(filename,0);
x = read_table(filename, ['%s' repmat('%f',1,numcols-1)],char(9),0);
P=[];
P.name = x.dat{1};
P.chr = x.dat{2};
P.start = x.dat{3};
P.end = x.dat{4};
C = cat(2,x.dat{5:end});

[P ord] = sort_struct(P,{'chr','start'});
C = C(ord,:);
