function plot_ziggs(orig_dat,Qzig,broad_level_p,broad_level_q,chr_start,centromere,chr_end)
  
  plot(1:length(orig_dat),orig_dat);
  hold on
  
  line(1:centromere-chr_start,repmat(broad_level_p,1,centromere-chr_start), ...
       'Color','red');
  
  line(centromere+1-chr_start:chr_end-chr_start,repmat(broad_level_q,1,chr_end-centromere), ...
       'Color','green');
  
  for j=1:size(Qzig,1)
    cur_seg = Qzig(j,:);
    line((cur_seg(2)-chr_start):(cur_seg(3)-chr_start),repmat(cur_seg(7), ...
                                                      1,cur_seg(3)- ...
                                                      cur_seg(2)+1), ...
         'Color','magenta','LineStyle','--');
    if cur_seg(6)<cur_seg(7)
      line(repmat(mean([cur_seg(2)-chr_start cur_seg(3)-chr_start]),1,length(cur_seg(6):.001: ...
                                                        cur_seg(7))), ...
           cur_seg(6):.001:cur_seg(7),'Color','cyan','LineStyle',':','Marker','+');
    else
      line(repmat(mean([cur_seg(2)-chr_start cur_seg(3)-chr_start]),1,length(cur_seg(7):.001: ...
                                                        cur_seg(6))), ...
           cur_seg(7):.001:cur_seg(6),'Color','cyan','LineStyle',':', ...
           'Marker','x');
    end
  end
  
  hold off
 
  
  
  
