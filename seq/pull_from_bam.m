function [R B S ri bi si rname] = pull_from_bam(bamname,chr,st,en,params)
% [R B S ri bi si rname] = pull_from_bam(bamname,chr,st,en,params)
%
% "st" and "en" may be vectors of the same length,
%   in which case results will be concatenated in the order requested
%
% "chr" may not be a vector
%
% R = reads
% B = bases
% S = sequence of reference genome
%
% ri,bi,si = partitions by requested-region (useful in vectorized format)
%            each tells the *last* item belonging to that requested-region
%
% rname = java String array containing names of reads
%
% params.quiet -- quiet flag [default = 0]
% params.refdir -- reference directory [default = '/xchip/tcga/gbm/analysis/lawrence/genome/hg18']
%
% the columns of R are:
%  (1) readgroup (lane)
%  (2) readnumber (name hash)
%  (3) whichpairmate (1 or 2--or -1 for unpaired)
%  (4) readstart (1-based in chr) (-200=unmapped)
%  (5) readend (1-based in chr) (for unmapped, = readlen-199)
%  (6) strand (0=plus, 1=minus) (-1=unmapped)
%  (7) number of mismatches (-1=unmapped)
%  (8) mapping quality (-1=unmapped)
%  (9) index into B (1-based)
%      Note: index of first base = R(i,9)
%            index of last base  = R(i+1,9)-1 = R(i,9)+R(i,5)-R(i,4)
%  (10) pairmate chr (-1=unmapped/unpaired)
%  (11) pairmate start (1-based, -1=unmapped/unpaired)
%  (12) pairmate strand (0=plus, 1=minus, -1=unmapped/unpaired)
% [the following columns included only if P.include_aux_cols==true]
%  (13) insertsize  (total #bp from beginning of left read to end of right read, -1 if unpaired/unmapped)
%  (14) isdup  (1 if has been marked in the BAM as a duplicate, 0 otherwise)
% 
%
% the columns of B are:          (Note: insertions are not represented!)
%  (1) base    -100=deletion
%              -1=N   1=A  2=C  3=G  4=T   base=reference (or read is unmapped)
%              63=N  65=A 66=C 67=G 68=T   base=non-reference
%  (2) base quality (-100=deletion)
%  (3) read index back to R
%  (4) base position (1-based in chr; negative numbers for unmapped reads)
%
% for retrieving >~100Kb you need to increase the amount of memory Matlab gives to java
% by having the file "java.opts" in the startup directory, containing the line "-Xmx4g".
% type "java.lang.Runtime.getRuntime.maxMemory" to confirm
%
% Mike Lawrence and Gaddy Getz 2009-08-10

if nargin==4 && isstruct(en)
  % ALTERNATE CALLING SYNTAX: pull_from_bam(bamname,chr,st,params)
  params = en;
  clear en;
end

if ~exist('params','var'), params=[]; end
params=impose_default_value(params,'quiet',1);
params=impose_default_value(params,'build','hg18');
buildno = interpret_build(params.build);
if buildno==18
  params=impose_default_value(params,'refdir','/xchip/cga1/lawrence/xchip_tcga_gbm_analysis_lawrence/genome/hg18');
elseif buildno==19
  params=impose_default_value(params,'refdir','/xchip/cga1/annotation/db/ucsc/hg19');
else
  error('unknown build: please specify refdir');
end

params=impose_default_value(params,'blacklist','none');
params=impose_default_value(params,'maxreads',[]);
params=impose_default_value(params,'include_unmapped_reads',true);
params=impose_default_value(params,'include_duplicate_reads',false);
params=impose_default_value(params,'include_aux_cols',false);

javaclasspath(split(get_jcp,':'));
import org.broadinstitute.cga.tools.seq.*;
import net.sf.samtools.*;
import java.io.*;
import java.lang.*;

if ~exist('en','var'), en=st; end
if any(size(st)>1) || any(size(en)>1)   % vector mode
  error('pull_from_bam: vectorized mode not yet implemented');
end

demand_file(bamname);
x = org.broadinstitute.cga.tools.seq.BamGrasp();
if params.quiet, x.setQuietModeOn; else x.setQuietModeOff; end
x.openFile(String(bamname),String(params.blacklist),String(params.refdir));

[R B S ri bi si rname] = BamGrasp_load_region(x,chr,st,en,params);

x.closeFile();

