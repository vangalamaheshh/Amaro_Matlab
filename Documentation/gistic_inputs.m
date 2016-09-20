Description of Inputs for run_gistic_from_seg_lite


Segmentation File (required)  [segfile]

	6 column file with NO HEADERS

        Sample Name ||  Chromosome || Start Pos. || End Pos.  || #of Markers || Value
        
Markers File (required) [markersfile]

	3 column file with NO HEADERS

	Marker Name || Marker Chromosome || Marker Positions

Reference Gene File (required) [refgenefile]

	Matlab file containing struct objects cyto and rg

Array List File (not required) [array_list_file]

	1 Column file with HEADER = ‘array’

	(list of arrays to match to names in seg file and use for gistic)

CNV_FILE (not required) [cnv_file]

	2 Options:

1) 2 column file: (NO HEADER)

Marker Name || CNV Label

2) 6 column file: (WITH HEADER)

Marker || Chormosome || CNV Start || CNV End || CNV Flank Start || CNV Flank End


Other Inputs (not required):  t_amp (float), t_del (float), join_segment_size (int), ext (string), qv_thresh (float), save_seg_data (bool), res (?), remove_X (bool)
