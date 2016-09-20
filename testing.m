for j=1:slength(e_seg)
    if e_seg.Segment_Mean(j)<-1 || e_seg.Segment_Mean(j)>1
    [distance_s, index_s]=min(abs(b_seg.xs-e_seg.xs(j)));
    [distance_e, index_e]=min(abs(b_seg.xe-e_seg.xe(j)));
    
    if distance_s<distance_e
        index=index_s;
    else
        index=index_e;
    end
    
    s_diff=abs(b_seg.xs(index)-e_seg.xs(j));
    e_diff=abs(b_seg.xe(index)-e_seg.xe(j));
    segment_overlap=overlap(b_seg.xs(index),b_seg.xe(index),e_seg.xs(j),e_seg.xe(j));
    overlap_b=segment_overlap/(b_seg.xe(index)-b_seg.xs(index));
    overlap_e=segment_overlap/(e_seg.xe(j)-e_seg.xs(j));
    if overlap_b>.5 && overlap_e>.5 && ((b_seg.Segment_Mean(index)<-.5 && e_seg.Segment_Mean(j)<-.5)||(b_seg.Segment_Mean(index)>.5&&e_seg.Segment_Mean(j)>.5))
        samples.shared_segments(ismember(samples.individual_id,individuals{i}))=samples.shared_segments(ismember(samples.individual_id,individuals{i}))+1;
    end
    end



end