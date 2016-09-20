function change = find_change_from_ref_tum1_tum2(st,en,ref,tum1,tum2)

if length(ref)~=length(tum1) || length(tum1)~=length(tum2) || length(tum2)~=length(st) || length(st)~=length(en), error('ref/tum1/tum2 different lengths'); end

fprintf('\n')
nx = length(ref);
change = cell(nx,1);
err_flag = 0;
for i=1:nx
  if ~mod(i,100000), fprintf('%d/%d ',i,nx); end
  isref1 = strcmp(tum1{i},ref{i});
  isref2 = strcmp(tum2{i},ref{i});
  hom = strcmp(tum1{i},tum2{i});
  isins = strcmp(ref{i},'-');
  isdel = ~isins & (strcmp(tum1{i},'-') | strcmp(tum2{i},'-'));
  coord_length_matches_ref_length = en(i)-st(i)+1==length(ref{i});

  if isref1 & isref2                  % non-mutation
    change{i} = ref{i};
  elseif isdel                        % deletion
  	if ~coord_length_matches_ref_length
  	  err_flag = 1;
  	  fprintf('Mutation coordinate error (line %d): Length of start/end coordinates must equal length of reference bases.\n%d\t%d\t%s\t%s\t%s\n',i,st(i),en(i),ref{i},tum1{i},tum2{i})
  	  change{i} = 'error';
  	else
  	  change{i} = ['-' ref{i}];
  	end
  elseif isins  % insertion
  	if en(i)-st(i) ~= 1
  	  err_flag = 1;
  	  fprintf('Mutation coordinate error (line %d): Insertion start/end coordinates must be adjacent bases!\n%d\t%d\t%s\t%s\t%s\n',i,st(i),en(i),ref{i},tum1{i},tum2{i})
  	  change{i} = 'error';
  	else
      if isref1
        change{i} = ['+' tum2{i}];
      else
        change{i} = ['+' tum1{i}];
      end
    end
  else                                % point-mutation
  	if ~coord_length_matches_ref_length
  	  err_flag = 1;
  	  fprintf('Mutation coordinate error (line %d): Length of start/end coordinates must equal length of reference bases.\n%d\t%d\t%s\t%s\t%s\n',i,st(i),en(i),ref{i},tum1{i},tum2{i})
  	  change{i} = 'error';
  	elseif length(tum1{i})~=length(ref{i}) || length(tum2{i})~=length(ref{i})
  	  err_flag = 1;
  	  fprintf('Mutation coordinate error (line %d): Length of alternate bases must equal length of refence bases.\n%d\t%d\t%s\t%s\t%s\n',i,st(i),en(i),ref{i},tum1{i},tum2{i})
  	  change{i} = 'error';
  	else
      if length(ref{i}) > 1
        if isref1
        	change{i} = ['~' tum2{i}];
        else
          change{i} = ['~' tum1{i}];
  	  end
      elseif isref1 & ~isref2
        change{i} = tum2{i};
      elseif isref2 & ~isref1
        change{i} = tum1{i};
      else % ~isref1 & ~isref2
        if hom   % tum1==tum2
          change{i} = tum1{i};
        else     % try both alleles
          change{i} = [tum1{i} '/' tum2{i}];
        end
      end
    end
  end
end
if err_flag
  fprintf('\nError in maflite format!  The faulty lines displayed above will not be included in output.\n');
end
if nx>=100000, fprintf('\n'); end
