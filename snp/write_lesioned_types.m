function write_lesioned_types(fname,C,regs,ts)
% ---
% $Id$
% $Date: 2008-05-15 10:22:33 -0400 (Thu, 15 May 2008) $
% $LastChangedBy: cmermel $
% $Rev$ 
  
  
  types = {C.sis.type};
  for i=1:length(types)
    cur_type = char(types(i));
    if cur_type(end) == ' '
      types(i) = cellstr(cur_type(1:length(cur_type)-1));
    end
  end
      
  unique_types = unique(types);
  num_each_type = zeros(length(unique_types),1);
  for j=1:length(unique_types)
    num_each_type(j) = length(strmatch(unique_types(j),types,'exact'));
  end
  
  
  for k=1:2
    sample_types{k} = zeros(length(regs{k}),length(unique_types));
    for i=1:length(regs{k})
      switch k
        case 1
         samples = find(C.dat(regs{k}(i).peak,:) > ts(k));
       otherwise
         samples = find(C.dat(regs{k}(i).peak,:) < -ts(k));
      end
      for j=1:length(unique_types)
        sample_types{k}(i,j) = length(strmatch(unique_types(j), types(samples),'exact'));
        sample_type_freq{k}(i,j) = sample_types{k}(i,j)/num_each_type(j);
      end
    end
  end
  
  for k=1:2
    count_types{k} = zeros(length(regs{k}),1);
    for i=1:length(regs{k})
      count = 0;
      for j=1:length(unique_types)
        if sample_type_freq{k}(i,j) > 0 
          count = count+1;
        end
      end
      count_types{k}(i) = count;
      
    end
  end
  
  f = fopen([fname],'w');
 
  for k=1:2
    switch k
     case 1
      fprintf(f,'Amplifications');
      fprintf(f,'\n');
      fprintf(f,'Chr: Start - End');     
      fprintf(f,'\t');
      fprintf(f,'resid_q-value');
      fprintf(f,'\n');
     otherwise
      fprintf(f,'\n');
      fprintf(f,'Deletions');
      fprintf(f,'\n');
      fprintf(f,'Chr: Start - End');
      fprintf(f,'\t');
      fprintf(f,'resid_q-value');
      fprintf(f,'\n');
    end
  
    for i = 1:length(regs{k})
      fprintf(f,'%s',strcat('chr', num2str(regs{k}(i).chrn), ':', ...
                  num2str(C.pos(regs{k}(i).peak_wide_st)), '-', ...
                      num2str(C.pos(regs{k}(i).peak_wide_en))));
      fprintf(f,'\t');
      fprintf(f,'%g',regs{k}(i).resid_qv);
      fprintf(f,'\t');
      [mx,mi] = sort(sample_type_freq{k}(i,:),'descend');
      ind = 1;
      while mx(ind) >0
        fprintf(f,'%s %s %2.3f %s %i',strvcat(cellstr(unique_types(mi(ind)))), ' | freq = ', ...
                mx(ind)*100, '% | num = ', sample_types{k}(i,mi(ind)));
        fprintf(f,'\t');
        ind=ind+1;
      end
      fprintf(f,'\n');
    end
    
  end
  
  fclose(f);
  
  
