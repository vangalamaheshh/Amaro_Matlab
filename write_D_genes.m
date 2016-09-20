function write_D_genes(fname,D,add_fields)

f=fopen(fname,'w');
txt=as_column(D.gacc);
if isfield(D,'gsymb')
  txt=[ txt as_column(D.gsymb) ];
end
txt=[txt as_column(D.gdesc) ];

for i=1:length(add_fields)
  fld=getfield(D,add_fields{i});
  if iscell(fld)
    fldtxt=fld;
  else
    fldtxt=cellstr(num2str(as_column(fld)));
  end
  txt=[txt fldtxt]; 
end

txt=txt';
fprintf(f,['%s' repmat('\t%s',1,size(txt,1)-1) '\n'],txt{:});

fclose(f);

