function regs=get_bootstrap_wide_peaks(regs,b,ci)

for k=1:2
  for i=1:length(regs{k})
    p_st=regs{k}(i).peak_wide_st;
    p_en=regs{k}(i).peak_wide_en;
    m_st=max(find(b{k}(1:(p_st-1))==0))+1;
    if isempty(m_st)
      m_st=1;
    end
    m_en=min(find(b{k}((p_en+1):end)==0))+p_en-1;
    if isempty(m_en)
      m_en=length(b{k});
    end
    max_val=max(b{k}(m_st:m_en));
    disp([ k i p_st p_en m_st m_en max_val]);
   if max_val==0
      disp([ k i p_st p_en m_st m_en max_val]);
      regs{k}(i).b_wide_peak_st=p_st;
      regs{k}(i).b_wide_peak_en=p_en;
    else
      thresh=floor((1-ci)*max_val);
      regs{k}(i).b_wide_peak_st=min(find(b{k}(m_st:m_en)>=thresh))+m_st-1;
      regs{k}(i).b_wide_peak_en=max(find(b{k}(m_st:m_en)>=thresh))+m_st-1;
    end
  end
end
