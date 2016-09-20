function write_mit_odf_file(fname,odf)

if ~isfield(odf,'data')
  error('odf must have a data field');
end

nl=char([13 10]);


f=fopen(fname,'w');

fprintf(f,['ODF 1.0' nl]);

fields=fieldnames(odf);
fields=fields(setdiff(1:length(fields),strmatch('data',fields)));
fields=fields(setdiff(1:length(fields),strmatch('col_names',fields)));
fields=fields(setdiff(1:length(fields),strmatch('col_types',fields)));
fields=fields(setdiff(1:length(fields),strmatch('n_header',fields)));
nf=length(fields);

% add column_names, column_types, data_lines (+3)
fprintf(f,['HeaderLines= %d' nl],nf+3); 

fprintf(f,['COLUMN_NAMES:']);
for i=1:length(odf.col_names)
  fprintf(f,'\t%s',odf.col_names{i});
end
fprintf(f,nl);

fprintf(f,['COLUMN_TYPES:']);
for i=1:length(odf.col_types)
  fprintf(f,'\t%s',odf.col_types{i});
end
fprintf(f,nl);

for i=1:nf
  val=getfield(odf,fields{i});
  if ischar(val)
    fprintf(f,['%s= %s' nl],fields{i},val);
  elseif length(val)>1
    fprintf(f,['%s= %s' nl],fields{i},num2str(val));
  else
    fprintf(f,['%s= %f' nl],fields{i},val);
  end
end

fprintf(f,['DataLines= %d' nl],length(odf.data{1}));

for i=1:length(odf.data{1})
  for j=1:length(odf.data)
    switch lower(odf.col_types{j})
     case 'string' 
      fprintf(f,'%s',odf.data{j}{i});
     case 'float'
      fprintf(f,'%f',odf.data{j}(i));
     case 'int'
      fprintf(f,'%d',odf.data{j}(i));      
     case 'boolean'
      if odf.data{j}(i)
        fprintf(f,'true');
      else
        fprintf(f,'false');
      end
    end
    if j<length(odf.data)
      fprintf(f,'\t');
    end
  end
  fprintf(f,nl);
end

fclose(f);

% EXAMPLE
%
% SDF 1.0
% HeaderLines= 8
% COLUMN_NAMES:	Samples	True Class	Predicted Class	Confidence	Correct?
% COLUMN_TYPES:	String	String	String	float	boolean
% Model= Prediction Results
% PredictorModel= KNN
% NumFeatures= 10
% NumCorrect= 17
% NumErrors= 0
% DataLines= 17
% cup15	0	0	0.43092585	true
% cup13	1	1	0.70015657	true
% cup14	1	1	0.5717054	true
% cup17	1	1	0.63318145	true
% cup12	1	1	0.11933289	true
% cup9	1	1	0.61582136	true
% cup8	1	1	0.66873986	true
% cup6	1	1	0.6302405	true
% cup22	1	1	0.69220066	true
% cup7	1	1	0.14071608	true
% cup10	1	1	0.61542803	true
% cup4	1	1	0.526848	true
% cup5	1	1	0.6296716	true
% cup2	1	1	0.671063	true
% cup11	1	1	0.7921765	true
% cup3	1	1	0.7399955	true
% cup1	1	1	0.6953054	true
