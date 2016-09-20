function [CL,regs,pvs]=after_segmentation(C,should_log,should_use_X,is_hmm,idx,smooth_sz,peak_sz,min_delta,ext, ...
                                          pv_thresh,segments_file,use_new,less_memory,ts,subtract_median,clear_segments,cnv_file,cyto)

if ~exist('segments_file','var') || isempty(segments_file)
  segments_file='output';
end

if ~exist('ext','var') || isempty(ext)
  ext=[];
end

if ~exist('should_log','var') || isempty(should_log)
  should_log=1;
end

if ~exist('should_use_X','var') || isempty(should_use_X)
  should_use_X=1;
end

if ~exist('pv_thresh','var') || isempty(pv_thresh)
  pv_thresh=0.25;
end

if ~exist('min_delta','var') || isempty(min_delta)
  min_delta=0;
end

if ~exist('is_hmm','var') || isempty(is_hmm)
  is_hmm=0;
end

if ~exist('smooth_sz','var') || isempty(smooth_sz)
  smooth_sz=3;
end

if ~exist('peak_sz','var') || isempty(peak_sz)
  peak_sz=5;
end

if ~exist('use_new','var') || isempty(use_new)
  use_new=0;
end

if ~exist('less_memory','var') || isempty(less_memory)
  less_memory=0;
end

if ~exist('ts','var') || isempty(ts)
  ts=[0.3 0.3];
end

if ~exist('subtract_median','var') || isempty(subtract_median)
  subtract_median=0;
end

if ~exist('cnv_file','var')
  cnv_file=[];
end

% keyboard
if less_memory
  n_parts=5;
  parts=get_parts(1:size(C.dat,2),n_parts);
  for i=1:n_parts    
    CLs{i}=preprocess_segmented_data(reorder_D_cols(C,parts{i}),segments_file,should_log,is_hmm,smooth_sz, ...
                                     peak_sz,min_delta,less_memory,parts{i},subtract_median,cnv_file);
    CLs{i}=rmfield_if_exists(CLs{i},{'sm1','sm2','sm2j','sm3','cbs_fixed','cbs','smooth','raw'});
  end
  clear C;
  CL=unite_Ds(CLs,'cols');
else
  CL=preprocess_segmented_data(C,segments_file,should_log,is_hmm,smooth_sz,peak_sz,min_delta,less_memory,[],subtract_median,cnv_file);
end

if exist('clear_segments','var') && ~isempty(clear_segments)
  for i=1:length(clear_segments)
    CL.dat(clear_segments{i}(1):clear_segments{i}(2),:)=repmat(CL.dat(clear_segments{i}(1),:), ...
                                                      diff(clear_segments{i})+1,1);
    disp(['cleared segement ' num2str(clear_segments{i})]);
  end
end

save(['CL_' num2str(size(CL.dat,2)) ext '.mat'],'CL','-V7.3');

if exist('idx','var') && ~isempty(idx)
  CL=reorder_D_cols(CL,idx);
end

if ~should_use_X
  CL=reorder_D_rows(CL,find(CL.chrn~=23));
else
  disp('Are you sure you have normalized against females?');
end

if (0)
  CL.dat(CL.dat<-0.30001)=-0.30001;
  CL.dat(CL.dat>0.30001)=0.30001;
end
keyboard
if ~use_new
  t1s=-0.3;
  t2s=0.3;
  %GG (06/03/22) was CL.smooth but CL.dat=CL.smooth 
  [h_amp,h_del,naamp,nadel]=glioma_perm(CL.dat,0, ...
                                        struct('method','log',...
                                               't1s',t1s,...
                                               't2s',t2s));
  
  
  [thst_amp,thst_del]=exact_permutations(CL,t1s,t2s);

  [regs,pvs]=generate_regs(CL,['peaks' ext],0.3,h_amp,h_del,naamp,nadel,thst_amp,thst_del,pv_thresh,501,0);
else
  score_type=struct('method','nxa','amp_thresh',ts(1),'del_thresh',-ts(2),'res',0.001);
%  score_type=struct('method','nxa','amp_thresh',ts(1),'del_thresh',-ts(2),...
%                    'res',0.0001,'max_segment_size',2000);
  [q,p,d,ads]=snp_score_permutations(CL,score_type,-1);
  save(['stats' ext '.mat'],'q','p','d','ads','score_type');
 
  qv_thresh=0.25;
  for k=1:2
    score_thresh(k)=min(ads{k}(find(q{k}<=qv_thresh)));
  end
  regs=generate_regs_by_peel_off(CL,ads,d,q,score_type,score_thresh,501,struct('method','origloo'));
  
  pvs=q;
   plot_snp_score([ num2str(size(CL.dat,2)) ext],CL,q,ads,d,0.25,1,1,'vert',cyto); % 0
end


