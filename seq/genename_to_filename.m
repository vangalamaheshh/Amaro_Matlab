function filename = genename_to_filename(gene)
    filename = upper(gene);
    filename(filename=='\')='_';
    filename(filename=='/')='_';
