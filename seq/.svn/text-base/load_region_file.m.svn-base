function R = load_region_file(region_list);
% Mike Lawrence 2009-10-16

numcols = get_colcount(region_list,0);

x = read_table(region_list,['%s' repmat('%f',1,numcols-1)],char(9),0);
R = [];
R.name = x.dat{1};
R.chr = x.dat{2};
R.start = x.dat{3};
R.end = x.dat{4};
if numcols<5
  R.membership = ones(length(R.name),1);
else
  R.membership = x.dat{5};
end
