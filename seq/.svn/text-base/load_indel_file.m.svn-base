function X = load_indel_file(fname)
% Mike Lawrence 2009-10-21

demand_file(fname);

tmp = load_lines(fname);
tmp = tmp(grep('.*SOMATIC.*CODING.*',tmp,1));
tmp = tmp(setdiff(1:length(tmp),grep('INVALID',tmp,1)));
X = parse(tmp,['^(chr\S*)\t(\d*)\t(\d*)\t(\S*)\t(SOMATIC)\t(\S*\t)?(CODING)\t(\S*)'],...
        {'chr','start','end','details1','status','details2','type','gene'});


