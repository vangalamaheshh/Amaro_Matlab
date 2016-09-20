function a=parse_affy_annot(a)

probeset_col=strmatch('Probe Set ID',a.headers,'exact');
a.probeset=a.data(:,probeset_col);

refseq_transcript_col=strmatch('RefSeq Transcript ID',a.headers,'exact');
a.refgene=regexp(a.data(:,refseq_transcript_col),'(\w*)+','tokens');

description_col=strmatch('Target Description',a.headers,'exact');
a.desc=a.data(:,description_col);

symbol_col=strmatch('Gene Symbol',a.headers,'exact');
a.symb=a.data(:,symbol_col);
