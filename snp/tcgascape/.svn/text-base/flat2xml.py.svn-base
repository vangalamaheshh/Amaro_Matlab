#!/usr/local/bin/python

'''
flat2xml.py
Convert flat files for TumorscapeTCGA to XML for new portal

Nam Pho - nampho@broadinstitute.org
April 2012

schum May-June 2012
'''

import os
import sys
import string, re
import datetime
import glob
from xml.etree.ElementTree import ElementTree

#! temporary global gene hash set for testing
#! (filters the genes in peaks by genes that have flat files)
gene_set = set()

# 2 or 3 arguments allowed
if len(sys.argv) != 5 :
    script = os.path.basename(sys.argv[0])
    sys.stderr.write("Usage: python %s <input parameter xml file> <base dir> <input dir> <output dir>\n" % script)
    sys.exit(1)

##--- FUNCTION: output a single gene to the analysis and metadata XML files

# global regexp for parsing gene and region tokens from first line of gene files
gpattern = re.compile(r'Result:\s+(\S+)\s+\((.+)\)')

def gene2xml(gene_file):

    #-- helper subfunction: dump <gene-alteration> element
    def dump_gene_alteration(scnas,alt_type,alt_desc):
        #-- helper for resorting
        def ranker(roi):
            return float(roi[5]) - (1 if roi[2]=="Yes" else 0) - (2 if roi[0]=="all_cancers" else 0)
        # resort, remove redundant "all_cancers" (work-around for bungled deletions 2/2011-10/2012)
        scnas = sorted(scnas,key=ranker);
        if scnas[1][0]=="all_cancers":
            scnas.pop(1)

        ANAL.write("\t\t<gene-alteration alteration-type='%s' description='%s'>\n" % (alt_type,alt_desc))
        rank = 1
        for roi in scnas:
            tissue = roi[0].strip()
            if tissue_dict.has_key(tissue):
                inpeak = ""
                if roi[2] == "Yes":
                    inpeak = "true"
                else:
                    inpeak = "false"
                ANAL.write("\t\t\t<gene-alteration-analysis in-driver-peak='%s' q-value='%s' " % (inpeak, roi[5]) )
                ANAL.write("rank='%d' " % rank)
                ANAL.write("overall-freq='%s' focal-freq='%s' high-level-freq='%s'>\n" % (roi[6], roi[7], roi[8]) )
                ANAL.write("\t\t\t\t<tissue-ref tissue-type='%s'/>\n" % disease_map[tissue])
                if not roi[3] == 'No peak on chromosome': 
                    ANAL.write("\t\t\t\t<peak-ref region='%s'/>\n" % (roi[3]))
                ANAL.write("\t\t\t</gene-alteration-analysis>\n")
                rank = rank + 1
        ANAL.write("\t\t</gene-alteration>\n")

    # load gene flat file input    
    FILE = open(gene_file, "r")
#    flat_content = []
    
#!    try:
#!        flat_content = FILE.readlines()
#!    finally:
#!        FILE.close()

    # figure out number of cancer types from number of lines
#!    ntissues = (len(flat_content) - 15) / 2

    # parse the gene and region from the first line of the gene file
    gmatch = gpattern.match(FILE.readline())
#!    gmatch = gpattern.match(flat_content[0])
    gene = gmatch.group(1)
    region = gmatch.group(2)

    # write out the gene data
    ANAL.write("\t<gene-analysis>\n")
    ANAL.write("\t\t<gene-ref symbol='%s' />\n\n" % (gene))
    META.write("\t\t<gene symbol='%s' region='%s'/>\n" %(gene,region))
    
    gene_set.add(gene)    # add gene to filter

    # skip blank line and read summary
    FILE.readline()
    gene_amp_desc = FILE.readline().replace("Amplification Summary:\t","",1).strip()
    gene_amp_desc = gene_amp_desc.replace("subtypes","cancer types");
 
    # skip 2 blank, 3 header lines
    FILE.readline()
    FILE.readline()
    FILE.readline()
    FILE.readline()
    FILE.readline()

    # get gene amp data
    amps = []
    line = FILE.readline().strip()
    while len(line) > 0:
#!        print("A>> %s\n",line.split("\t")[0])
        amps.append(line.split("\t"))
        line = FILE.readline().strip()
#!    for cancer in flat_content[8:ntissues+8]:
#!        amps.append( cancer.strip().split("\t") )    
    # output gene amp alteration element
#!    gene_amp_desc = flat_content[2].replace("Amplification Summary:\t","",1).strip()
#!    gene_amp_desc = gene_amp_desc.replace("subtypes","cancer types");
    dump_gene_alteration(amps,"amplification",gene_amp_desc)

    # read deletion summary
    gene_del_desc = FILE.readline().replace("Deletion Summary:\t","",1).strip()
    gene_del_desc = gene_del_desc.replace("subtypes","cancer types");

    # skip 2 blank, 3 header lines
    FILE.readline()
    FILE.readline()
    FILE.readline()
    FILE.readline()
    FILE.readline()

    # get gene del data
    dels = []
    line = FILE.readline().strip()
    while len(line) > 0:
#!        print("D>> %s\n",line.split("\t")[0])
        dels.append(line.split("\t"))
        line = FILE.readline().strip()
#!    for cancer in flat_content[(ntissues+15):(2*ntissues+15)]:
#!        dels.append( cancer.strip().split("\t") )

    # output gene del data
#!    gene_del_desc = flat_content[ntissues+9].replace("Deletion Summary:\t","",1).strip()
#!    gene_del_desc = gene_del_desc.replace("subtypes","cancer types");
    dump_gene_alteration(dels,"deletion",gene_del_desc)

    ANAL.write("\t</gene-analysis>\n")
    FILE.close()

##--- FUNCTION: read tissue types
def suck_diseases(tissue_file):
    tissue_dict = {}
    FILE = open(tissue_file, "r")
    tissue_lines = []
    print("Processing %s..." % tissue_file)
    try: 
        tissue_lines = FILE.readlines()
    finally:
        FILE.close()
    tissue_lines.pop(0)  # remove header
    for tissue_record in tissue_lines:
        tissue = tissue_record.split("\t")[0].strip()
        tissue_desc = tissue_record.split("\t")[1].strip()
        tissue_dict[tissue] = tissue_desc
    return tissue_dict

##--- FUNCTION: read optional disease name map 
def suck_disease_map(disease_map_file,tissue_dict):
    disease_map = {}
    # default is identity map
    for key in tissue_dict.keys():
        disease_map[key] = key
    # if we have a file, read in the mapping
    if len(disease_map_file) == 1:
        FILE = open(disease_map_file.pop(0), "r")
        tissue_lines = []
        print("Processing disease map file %s..." % tissue_file)
        try: 
            text_lines = FILE.readlines()
        finally:
            FILE.close()
        text_lines.pop(0)  # remove header
        for mapping in text_lines:
            old_name = mapping.split("\t")[0].strip()
            new_name = mapping.split("\t")[1].strip()
            disease_map[old_name] = new_name
    else:
        if len(disease_map_file) > 1:
            print("More than one disease map file!")
    return disease_map

##--- FUNCTION: read all the peaks from amp or del file and group by tissue type
# returns a dictionary mapping disease name to a list of flat file lines
def suck_peaks(peak_file):
    FILE = open(peak_file, "r")
    peak_lines = []
    
    # load input peaks file
    print("Processing %s..." % peak_file)
    try: 
        peak_lines = FILE.readlines()
    finally:
        FILE.close()
    # remove header
    peak_lines.pop(0)
    # create lookup by type
    by_type = {}
    for peak in peak_lines:
        tissue = peak.split("\t")[0].strip()
        if by_type.has_key(tissue):
            by_type[tissue].append(peak)
        else:
            by_type[tissue] = [peak]
    return by_type

##--- FUNCTION: output the cancer type, its enclosed peaks, and (to metadata) their enclosed genes
peak_set = set()  # global peak set ensures unique peak regions in metadata
def dump_peaks(alt_type,tissue_peaks):
    ANAL.write("\t\t<tissue-alteration alteration-type='%s'>\n" % alt_type)
    for peak in tissue_peaks:
        p = peak.split("\t")
        region = p[1]
        ANAL.write("\t\t\t<tissue-alteration-analysis q-value='%s' overall-freq='%s' focal-freq='%s' high-level-freq='%s'>\n" %(p[2],p[3],p[4],p[5])) 
        ANAL.write("\t\t\t\t<peak-ref region='%s'/>\n" % region)
        ANAL.write("\t\t\t</tissue-alteration-analysis>\n")
        # emit unique peak regions to metadata
        if region not in peak_set:
            peak_set.add(region)
            META.write("\t<peak region='%s'>\n" % region) # (no description)
            peak_genes = p[6].split("; ")
            for peak_gene in peak_genes:
                if peak_gene in gene_set:  #! temporary filter
                    META.write("\t\t<gene-ref symbol='%s'/>\n" % peak_gene)
            META.write("\t</peak>\n")
    ANAL.write("\t\t</tissue-alteration>\n")

##---- FUNCTION write gistic-analysis description
def write_analysis_description(OUTFILE,tscape_vertext,gistic_vertext,genome_text) :
    OUTFILE.write("Tumorscape %s was run with the following GISTIC parameters:&lt;br/>" % tscape_vertext)
    OUTFILE.write("&lt;table>&lt;br/>&lt;tr>&lt;th>Parameter&lt;/th>&lt;th>Value&lt;/th>&lt;/tr>")

    # abbreviatory for html row
    def html_row2(col1,col2):
        OUTFILE.write("&lt;tr>&lt;td>%s&lt;/td>&lt;td>%s&lt;/td>&lt;/tr>" % (col1,col2) )
    # output the rows
    html_row2('core GISTIC version',gistic_vertext)
    html_row2('reference genome build',genome_text)
    html_row2('amplification threshold',gistic_params['t_amp'])
    html_row2('deletion threshold',gistic_params['t_del'])
    html_row2('high-level amplification threshold',tscape_params['ht_amp'])
    html_row2('high-level deletion threshold',tscape_params['ht_del'])
    html_row2('broad length cutoff',gistic_params['broad_len_cutoff'])
    html_row2('peak confidence level',gistic_params['conf_level'])
    html_row2('cap',gistic_params['cap'])
    html_row2('gene-GISTIC',gistic_params['do_gene_gistic'])
    html_row2('arm-level peel-off',gistic_params['arm_peeloff'])
    html_row2('significance threshold',gistic_params['qv_thresh'])
    html_row2('join segment size',gistic_params['join_segment_size'])
    html_row2('X chromosome removed',gistic_params['remove_X'])
    html_row2('maximum segments/sample',gistic_params['max_segs_per_sample'])
    html_row2('minimum samples/disease',tscape_params['min_samples'])
    OUTFILE.write("&lt;/table>&lt;br/>")

######################
###  MAIN PROGRAM  ###
######################

# command line arguments
parm_fname = sys.argv[1]
base_dir = sys.argv[2]
#!!!TODO next two args could default to standard subdirectories of base_dir
input_dir = sys.argv[3] # contains 20K+ flat files
output_dir = sys.argv[4] # contains metadata and analysis XML load files

# get input parameters from XML file
inparams = ElementTree()
inparams.parse(parm_fname)
inparam_root = inparams.getroot()
gistic_params = inparam_root.find('gistic-run-parameters').attrib
tscape_params = inparam_root.attrib

# read in tumorscape/gistic version text
VERF = open(base_dir + os.sep + "version.txt")
try:
    vertext = VERF.readlines()
finally:
    VERF.close()

run_id = tscape_params['run_id']
run_version = tscape_params['run_version']
run_description = tscape_params['run_description']

# get refgene text
igv_genome = 'hg19'
igv_genome2 = 'human genome reference build GRCh37/hg19'

# set genome tag in metadata analysis
if tscape_params.has_key('igv_genome'):
    igv_genome = tscape_params['igv_genome']
if tscape_params.has_key('igv_genome2'):
    igv_genome2 = tscape_params['igv_genome2']

# open analysis and metadata files
ANAL = open("%s/tumorscape-analysis-%s.xml" % (output_dir,run_id), "w")
META = open("%s/tumorscape-metadata-%s.xml" % (output_dir,run_id), "w")

# output analysis file header and open analysis tag
ANAL.write("<?xml version='1.0' encoding='UTF-8'?>\n\n")
ANAL.write("<!DOCTYPE gistic-analysis PUBLIC '-//Broad Institute//DTD Gistic Analysis 1.2//EN' 'http://www.broadinstitute.org/gistic/gistic-analysis-1.2.dtd'>\n\n")
ANAL.write("<gistic-analysis version='%s' description='" % run_version)
write_analysis_description(ANAL,vertext[0].strip(),vertext[1].strip(),igv_genome)
ANAL.write("' igv-dir='%s'>\n" % run_id)

# output metadata file header and open metadata tag
META.write("<?xml version='1.0' encoding='UTF-8'?>\n\n")
META.write("<!DOCTYPE gistic-metadata PUBLIC '-//Broad Institute//DTD Gistic Metadata 1.2//EN' 'http://www.broadinstitute.org/gistic/gistic-metadata-1.2.dtd'>\n\n")
META.write("<gistic-metadata version='%s' description='%s'>\n" %(run_version,run_description))
# TODO analysis/metadata run description

# load tissue types
tissue_file = glob.glob(input_dir + "/GCM_disease_summary*.txt").pop(0)
tissue_dict = suck_diseases(tissue_file)

# load name mapping file (deals with evolving nomenclature)
disease_map_file = glob.glob(input_dir + "/GCM_disease_name_map*.txt")
disease_map = suck_disease_map(disease_map_file,tissue_dict)

# read in the peaks (and tissue descriptions)
amp_peak_file = glob.glob(input_dir + "/GCM_amp_regions*.txt").pop(0)
amps_by_type = suck_peaks(amp_peak_file)

del_peak_file = glob.glob(input_dir + "/GCM_del_regions*.txt").pop(0)
dels_by_type = suck_peaks(del_peak_file)

## output gene data

META.write("\t<genome build='%s' description='%s'>\n" % (igv_genome,igv_genome2))
flat_gene_files = os.listdir(input_dir)
print("processing %d flat files" % len(flat_gene_files))
for gene in flat_gene_files:
    if re.search("GCM_", gene) is None and not re.search("\.txt$",gene) is None:
        # gene summary flat file
        gene2xml(input_dir + "/" + gene)
    else:
        # special amp peak, del peak, or tissue type flat file
        print("skipping %s..." % gene)

META.write("\t</genome>\n")

# output tissue type metadata
for tissue in tissue_dict.keys():
    disease_desc = tissue_dict[tissue].strip()
    disease_name = disease_map[tissue]
    # translate old disease names to standardized ones
    for oldname in disease_map.keys():
        newname = disease_map[oldname]
        # (NOTE: () >< qualifiers avoid multiple substitutons when disease names are substrings)
        disease_desc = disease_desc.replace("("+oldname+")","("+newname+")")
        disease_desc= disease_desc.replace(">"+oldname+"<",">"+newname+"<")
    # translate ">"s so they can be put in HTML
    disease_desc = disease_desc.replace("<","&lt;")
    # retroactive consistency of terminology
    disease_desc = disease_desc.replace("cancer subtype","cancer type")
    disease_desc = disease_desc.replace("Tissue Type","Cancer Type")
    # output disease entry to metadata
    META.write("\t<tissue tissue-type='%s' description='%s'/>\n" % (disease_name,disease_desc))

# output peak data for disease
for tissue in tissue_dict.keys():
    disease = disease_map[tissue]
    ANAL.write("\t<tissue-analysis>\n")
    ANAL.write("\t\t<tissue-ref tissue-type='%s'/>\n" %(disease))
    if amps_by_type.has_key(tissue):
        dump_peaks('amplification',amps_by_type[tissue])
    if dels_by_type.has_key(tissue):
        dump_peaks('deletion',dels_by_type[tissue])
    ANAL.write("\t</tissue-analysis>\n")

# close analysis/metadata tags/files
ANAL.write("</gistic-analysis>\n")
ANAL.close()
META.write("</gistic-metadata>\n")
META.close()
