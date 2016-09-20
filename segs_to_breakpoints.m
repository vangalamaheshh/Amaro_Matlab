function [breakpoints_e breakpoints_b]=segs_to_breakpoints(e_seg,b_seg,targets,thresh,n_seg)

e_seg.diff_start=zeros(slength(e_seg),1);
e_seg.diff_end=zeros(slength(e_seg),1);

b_seg.diff_start=zeros(slength(b_seg),1);
b_seg.diff_end=zeros(slength(b_seg),1);

for i=1:slength(e_seg)
    if i>1&&i<slength(e_seg)
        
        if abs(e_seg.corrected_total_cn(i)-e_seg.corrected_total_cn(i-1))>thresh && isequal(e_seg.Chromosome{i},e_seg.Chromosome{i-1})
            e_seg.diff_start(i,1)=1;
        end
        
        if abs(e_seg.corrected_total_cn(i)-e_seg.corrected_total_cn(i+1))>thresh && isequal(e_seg.Chromosome{i},e_seg.Chromosome{i+1})
            e_seg.diff_end(i,1)=1;
        end
    
    elseif i==1
    if abs(e_seg.corrected_total_cn(i)-e_seg.corrected_total_cn(i+1))>thresh
        e_seg.diff_end(i,1)=1;
    end
    elseif i==slength(e_seg)
    if abs(e_seg.corrected_total_cn(i)-e_seg.corrected_total_cn(i-1))>thresh
        e_seg.diff_start(i,1)=1;
    end
    
    end

end


for i=1:slength(b_seg)
    if i>1&&i<slength(b_seg)
        
    if abs(b_seg.corrected_total_cn(i)-b_seg.corrected_total_cn(i-1))>thresh && isequal(b_seg.Chromosome{i},b_seg.Chromosome{i-1})
        b_seg.diff_start(i,1)=1;
    end
    if abs(b_seg.corrected_total_cn(i)-b_seg.corrected_total_cn(i+1))>thresh && isequal(b_seg.Chromosome{i},b_seg.Chromosome{i+1})
        b_seg.diff_end(i,1)=1;
    end
    
    elseif i==1
    if abs(b_seg.corrected_total_cn(i)-b_seg.corrected_total_cn(i+1))>thresh
        b_seg.diff_end(i,1)=1;
    end
    elseif i==slength(b_seg)
    if abs(b_seg.corrected_total_cn(i)-b_seg.corrected_total_cn(i-1))>thresh
        b_seg.diff_start(i,1)=1;
    end
    
    end

end



e_seg.dCN=[0;diff(e_seg.corrected_total_cn)];
b_seg.dCN=[0;diff(b_seg.corrected_total_cn)];

breakpoints_e.pos=[e_seg.xs(e_seg.diff_start==1)];%;e_seg.xe(e_seg.diff_end==1)];
breakpoints_b.pos=[b_seg.xs(b_seg.diff_start==1)];%;b_seg.xe(b_seg.diff_end==1)];

for i=1:slength(breakpoints_b)
    breakpoints_b.Chromosome{i,1}=b_seg.Chromosome{ismember(b_seg.xs,breakpoints_b.pos(i))};
    breakpoints_b.Breakpoint_Position{i,1}=b_seg.Startbp{ismember(b_seg.xs,breakpoints_b.pos(i))};
end
for i=1:slength(breakpoints_e)
    breakpoints_e.Chromosome{i,1}=e_seg.Chromosome{ismember(e_seg.xs,breakpoints_e.pos(i))};
    breakpoints_e.Breakpoint_Position{i,1}=e_seg.Startbp{ismember(e_seg.xs,breakpoints_e.pos(i))};
end

breakpoints_n.pos=n_seg.xb;
breakpoints_e.dCN=e_seg.dCN(e_seg.diff_start==1);
breakpoints_b.dCN=b_seg.dCN(b_seg.diff_start==1);

for i=1:slength(breakpoints_e)

    
    [v s]=min(abs(targets.xmid-breakpoints_e.pos(i)));
    s=max([1 s-2]);
    breakpoints_e.target_range_start(i,1)=targets.xs(s);
    
    [v e]=min(abs(targets.xmid-breakpoints_e.pos(i)));
    e=min([slength(targets) e+2]);
    breakpoints_e.target_range_end(i,1)=targets.xe(e);
    
end


for i=1:slength(breakpoints_n)

    
    [v s]=min(abs(targets.xmid-breakpoints_n.pos(i)));
    s=max([1 s-2]);
    breakpoints_n.target_range_start(i,1)=targets.xs(s);
    
    [v e]=min(abs(targets.xmid-breakpoints_n.pos(i)));
    e=min([slength(targets) e+2]);
    breakpoints_n.target_range_end(i,1)=targets.xe(e);
    
end
% overlap=zeros(slength(breakpoints_e),1);
% for i=1:slength(breakpoints_e)
%    if overlap(i)==0
%     o=find(breakpoints_e.target_range_start<breakpoints_e.pos(i)&breakpoints_e.target_range_end>breakpoints_e.pos(i)&breakpoints_e.pos~=breakpoints_e.pos(i));
%         if ~isempty(overlap)
%          overlap(o)=1;
%         end
%    end
% end
% breakpoints_e=reorder_struct(breakpoints_e,overlap==0);

for i=1:slength(breakpoints_b)

    
    [v s]=min(abs(targets.xmid-breakpoints_b.pos(i)));
    s=max([1 s-2]);
    breakpoints_b.target_range_start(i,1)=targets.xs(s);
    
    [v e]=min(abs(targets.xmid-breakpoints_b.pos(i)));
    e=min([slength(targets) e+2]);
    breakpoints_b.target_range_end(i,1)=targets.xe(e);
    
    
end




for i=1:slength(breakpoints_n)
    if ~isempty(find(breakpoints_e.target_range_start<breakpoints_n.pos(i) & breakpoints_n.pos(i)<breakpoints_e.target_range_end))
        
        breakpoints_e=reorder_struct(breakpoints_e,~(breakpoints_e.target_range_start<breakpoints_n.pos(i) & breakpoints_n.pos(i)<breakpoints_e.target_range_end));
        
        
    end
end

for i=1:slength(breakpoints_n)
    if ~isempty(find(breakpoints_b.target_range_start<breakpoints_n.pos(i) & breakpoints_n.pos(i)<breakpoints_b.target_range_end))        
        breakpoints_b=reorder_struct(breakpoints_b,~(breakpoints_b.target_range_start<breakpoints_n.pos(i) & breakpoints_n.pos(i)<breakpoints_b.target_range_end));
        
        
    end
end




% 
% overlap=zeros(slength(breakpoints_b),1);
% for i=1:slength(breakpoints_b)
%    if overlap(i)==0
%     o=find(breakpoints_b.target_range_start<breakpoints_b.pos(i)&breakpoints_b.target_range_end>breakpoints_b.pos(i)&breakpoints_b.pos~=breakpoints_b.pos(i));
%         if ~isempty(overlap)
%          overlap(o)=1;
%         end
%    end
% end
% breakpoints_b=reorder_struct(breakpoints_b,overlap==0);
end