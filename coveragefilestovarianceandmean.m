 
total=zeros(8,1);

%ss_cov=load_struct('Coverage_Metrics1.txt.sample_summary');
%int_sum=load_struct('Coverage_Metrics1.txt.sample_interval_summary');


for f=1:8    
    
covsum=sprintf('Coverage_Metrics%s.txt.sample_summary',num2str(f))
intsum=sprintf('Coverage_Metrics%s.txt.sample_interval_summary',num2str(f))

ss_cov=load_struct(covsum);
int_sum=load_struct(intsum);


for i=1:slength(ss_cov)
    if ~isempty(regexp(ss_cov.sample_id{i},'-CN','match'))
    filt(i)=1;
    else   
    filt(i)=0;
    end
end

 ss_cov=reorder_struct(ss_cov,filt==1)
 ss_cov.mean=cellfun(@str2num,ss_cov.mean);
 
 
 
 for i=1:slength(int_sum)
    nn=regexp(int_sum.Target{i},':','split');
    n=regexp(nn{2},'-','split');    
    if size(n,2)>1
    s=str2num(n{1});
    e=str2num(n{2});
    total(f)=total(f)+(e-s);
    end
 end

 if f==1
 for i=1:slength(ss_cov)
     if ~isempty(regexp(ss_cov.sample_id{i},'CN-N','match'))
total_coverage.(sprintf('CN%sNormal',char(ss_cov.sample_id{i}(1:3))))=zeros(8,1);  
     else
         total_coverage.(sprintf('CN%s',char(ss_cov.sample_id{i}(1:3))))=zeros(8,1);  
     end
 end
 
 end
 
 for i=1:slength(ss_cov)
     if ~isempty(regexp(ss_cov.sample_id{i},'CN-N','match'))
total_coverage.(sprintf('CN%sNormal',char(ss_cov.sample_id{i}(1:3))))(f,1)=ss_cov.mean(i)
     else
         total_coverage.(sprintf('CN%s',char(ss_cov.sample_id{i}(1:3))))(f,1)=ss_cov.mean(i)
 end
 
 
 
 
 
    end
   
end