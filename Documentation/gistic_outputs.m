Description of Outputs for run_gistic_from_seg_lite

All Lesions File:  Locations of gistic peaks with peak info and amplitudes 
       (see /xchip/cancergenome/data/Jen/Glioma/GISTICtest/all_lesions_file_.txt)
        
        Tab delimited file with (8 + Number of Samples) columns:
        
Lesion Name | Location Descriptor | Wide Limits | Peak Limits | Region Limits| Qvalues | Broad or Focal | Amplitude Threshold Key | Sample Data…

There are three groups of rows.  The first group contains the gistic peaks with a 0, 1, or 2 in the sample data to indicate whether that sample contains the peak.  (See amplitude threshold key for details.)  The second group contains the gistic peaks with the log2 ratios from gistic.  The third group contains only the peaks which are broad or both broad and focal with 1s and 0s similar to group 1 above.

Scores File:  Gistic output for all segments

	Tab delimited file with 8 columns:

	Type | Chromosome | Start | End | q-value | score | amplitude | frequency
