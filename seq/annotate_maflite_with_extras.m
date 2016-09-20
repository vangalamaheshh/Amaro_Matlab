function annotate_maflite_with_extras(infile,outfile,P)

if ~exist('P','var')
  P=[];
end

P=impose_default_value(P,'build','hg18'); % hg18_v2 hg19
P=impose_default_value(P,'dbsnp_rod','/xchip/cga1/annotation/db/dbsnp/dbsnp_130_hg18.rod');
P=impose_default_value(P,'refseq_shelf','/xchip/cga1/annotation/db/genbank/genbank_transcript_Release_39.shelf');
P=impose_default_value(P,'sequencing_source','Capture'); % Capture
P=impose_default_value(P,'phase','Phase_I'); % Phase_I, Phase_II ...
P=impose_default_value(P,'platform','Illumina_GAIIx'); % Illumina_GAIIx or Illumina_HiSeq


%#Launch Matlab and run the following commands:
%p=genpath('~/CancerGenomeAnalysis/trunk/matlab/')
%addpath(p)
%annotate_maflite(<input maf>, <output maf>, <build>)
%
%#Exit Matlab and run the following:
%set p = ~/CancerGenomeAnalysis/trunk/analysis_pipeline/genepattern/modules
%set s = $p/FormatMafliteToTCGASpecs/
%set dbsnp = /xchip/cga1/annotation/db/dbsnp/dbsnp_130_hg18.rod
%set refseq = /xchip/cga1/annotation/db/genbank/genbank_transcript_Release_39.shelf
%sh $p/FormatMafliteToTCGASpecs/Reformat_maf_to_TCGA.sh <input_maf> <output_maf> $s $dbsnp $refseq Capture Phase_I Illumina_GAIIx
%$p/AddUniProtToMaf/add_uniprot_to_maf.py <input_maf> <output_maf>
%cat <input_maf> | $p/AddRefseqCDDAlignmentsToMafAnnotated/add_refseq_cdd_alignments_to_mafannotated.py -x 0 | ...
%    $p/AddCosmicToMafAnnotated/add_cosmic_to_mafannotated.py -x > <output_maf>

%% annotate_maflite

final_outfile=outfile;

outfile_1=[ infile '.1'];
annotate_maflite(infile,outfile_1,P.build);

modules_path='~/CancerGenomeAnalysis/trunk/analysis_pipeline/genepattern/modules/';
spec_module_path=[ modules_path 'FormatMafliteToTCGASpecs/'];


outfile_2=[ infile '.2'];
unix_cmd=['sh ' modules_path 'FormatMafliteToTCGASpecs/Reformat_maf_to_TCGA.sh ' outfile_1 ' ' outfile_2 ' ' ...
          spec_module_path ' ' P.dbsnp_rod ' ' P.refseq_shelf ' ' P.sequencing_source ' ' P.phase ' ' P.platform];
disp(unix_cmd);
unix(unix_cmd);

disp('Hit a key');
pause 

outfile_3=[ infile '.3'];
unix_cmd=[ modules_path 'AddUniProtToMaf/add_uniprot_to_maf.py ' outfile_2 ' ' outfile_3 ];
disp(unix_cmd);
unix(unix_cmd);

unix_cmd=[ 'cat ' outfile_3 ' | ' modules_path 'AddRefseqCDDAlignmentsToMafAnnotated/' ...
           'add_refseq_cdd_alignments_to_mafannotated.py -x 0 | ' ...
           modules_path 'AddCosmicToMafAnnotated/add_cosmic_to_mafannotated.py -x > ' final_output ];
disp(unix_cmd);
unix(unix_cmd);

