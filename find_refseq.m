function i=find_refseq(rg,refseq)

i=strmatch(refseq,strvcat(rg(:).regseq));
