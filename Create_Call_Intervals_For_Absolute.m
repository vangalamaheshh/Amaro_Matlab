function [ output_args ] = Create_Call_Intervals_For_Absolute( Individual_ids, OxoG_filtered_Mafs_paths, mouse , target, pair_mapping_file)

addpath /home/unix/amaro/
addpath /xchip/cga2/amaro/CancerGenomeAnalysis/trunk/matlab
addpath /xchip/cga2/amaro/CancerGenomeAnalysis/trunk/matlab/mike
mouse
if isequal(mouse,'mm9') || isequal(mouse,'mm10')
mouse=1;
else
mouse=0;
end


Individual_ids
if exist(Individual_ids, 'file')
Individual_ids=readtext(Individual_ids, char(9));
%Individual_ids={Individual_ids{:,2}}

Individual_ids={Individual_ids{1:end,2}};
else 
  Individual_ids=str2cell(Individual_ids);
end


%OxoG_filtered_Mafs_paths
OxoG_filtered_Mafs_paths=readtext(OxoG_filtered_Mafs_paths, char(9));
%OxoG_filtered_Mafs_paths={OxoG_filtered_Mafs_paths{:,2}};

%OxoG_filtered_Mafs_paths={OxoG_filtered_Mafs_paths{2:end}};
%[uni_individuals x y]=unique(Individual_ids);
pair_mapping_files=readtext(pair_mapping_file,char(9));
uni_individuals=Individual_ids;
counter=1;
Num_Samp_I=zeros(size(uni_individuals));

for i=1:size(pair_mapping_files,1)
    map=readtext(pair_mapping_files{i,2},char(9));
 for j=1:size(uni_individuals,2)
     for k=1:size(map,1)
         if isequal(map{k,2},uni_individuals{j})
             master_map{counter,1}=map{1,1};
             master_map{counter,2}=uni_individuals{j};
             counter=counter+1;
         end
     end
 end
 
end

for i=1:size(master_map,1)
    for j=1:size(OxoG_filtered_Mafs_paths,1)
        if ~isempty(strfind(OxoG_filtered_Mafs_paths{j,1},master_map{i,2}))
            master_map{i,3}=OxoG_filtered_Mafs_paths{j,2};
        end
    end
end
%Num_Samp_I=zeros(size(uni_individuals));


%for i =1:size(uni_individuals,2)
%for j=1:size(OxoG_filtered_Mafs_paths,2)
%if ~isempty(regexp(OxoG_filtered_Mafs_paths{j},uni_individuals{i}))
%Num_Samp_I(i)=Num_Samp_I(i)+1;
%end
%end
%end



%Num_Samp_I=count(y);
Num_Samp_I=count({master_map{:,1}});
uni_individuals=unique({master_map{:,1}});

vounter=1;
for i=1:size(uni_individuals,2)
  i  
  if Num_Samp_I(i)>0
  for j=1:Num_Samp_I(i)
        new_maf=master_map{vounter,3};
        vounter
  vounter=vounter+1;
        
if Num_Samp_I(i)==1 && ~isempty(new_maf)
	CreateForceCallIntervalList_Oxo(new_maf,new_maf,sprintf('%s.intervals',uni_individuals{i}),mouse);
	outfile{i}=sprintf('%s.intervals',uni_individuals{i}); 
    lines(i)=size(readtext(outfile{i}),1);
	  else
            if ~isempty(new_maf)
                if j==1;
                I_maf=load_struct(new_maf);
                end
                if j<Num_Samp_I(i)
                N_maf=load_struct(new_maf);
                I_maf=concat_structs_keep_all_fields({I_maf,N_maf});
            
                end
                if j==Num_Samp_I(i)
		  file_name=sprintf('%s.maf.oxoG',uni_individuals{i});
save_struct(I_maf,file_name);
               % other_maf_text=readtext(file_name,char(9));
               % maf2{1,1}='## Oncotator';
               
                      %  for h=1:size(other_maf_text,1)
                      %      for k=1:size(other_maf_text,2)
                     %           maf2{1+h,k}=other_maf_text{h,k};
                    %        end
                   %     end
                  %cell2csv(file_name,maf2,char(9));
file_name
new_maf
		     sprintf('%s.intervals',uni_individuals{i})
mouse
                CreateForceCallIntervalList_Oxo(file_name,new_maf,sprintf('%s.intervals',uni_individuals{i}),mouse);
                
outfile{i}=sprintf('%s.intervals',uni_individuals{i});
        lines(i)=size(readtext(outfile{i}),1);      %  unix(['awk {print "chr" $0}' outfile '>' outfile '.intervals' ]);
		delete(file_name);
                end
            end    
        end
         
  end
    
  clear I_maf N_maf new_maf 
  end
end
if mouse
%unix('chmod +x if_mouse.sh')
%unix('for file in *.intervals; do awk '{print "chr"$0'} $file > $file.mouse; done')
end

annotation_file{1,1}='individual_id';
annotation_file{1,2}=target;
line_file{1,1}='individual_id';
line_file{1,2}='number_of_intervals_force'
for i=1:size(uni_individuals,2)
    if Num_Samp_I(i)>0
annotation_file{i+1,1}=uni_individuals{i};
line_file{i+1,1}=uni_individuals{i};
%if ~mouse
annotation_file{i+1,2}=sprintf('%s/%s',pwd,outfile{i});
line_file{i+1,2}=lines(i);
%end
if mouse
%annotation_file{i+1,2}=sprintf('%s/%s.mouse',pwd,outfile{i});
end
    end
end
annotation_file
cell2csv('force_calls_annot.txt',annotation_file,char(9))
cell2csv('force_calls_lines.txt',line_file,char(9));

exit

end

