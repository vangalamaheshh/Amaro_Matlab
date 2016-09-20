function alpha = calculate_alpha(focals,amp_thresh,del_thresh,max_amp,max_del,pseudocount,plot_hists)
%CALCULATE_ALPHA estimate dependence of event frequency on amplitude from data
%    
%   ALPHA = calculate_alpha(FOCALS,T_AMP,T_DEL,MAX_AMP,MAX_DEL,PDEUDOCOUNT,PLOT_HISTS);
%
% Returns the exponential factors 'alpha' derived from fitting an
% exponential curve to the amplitude of focal events. ALPHA is a
% vector of the amplification alpha followed by the deletion alpha.
% 
% FOCALS is a structure containing copy number data matrices for
% amplification (.amp), deletion (.del), amplification over
% deletion (.aod) and deletion over amplification (.doa).
%
% T_AMP and T_DEL are amplification and deletion cutoff values.
%
% MAX_AMP and MAX_DEL are the extreme amplification and deletion
% ranges to be used in the estimation. 
% 
% PSEUDOCOUNT is a percent of the total amplitude added to the
% score distribution to eliminate zeros.
%
% PLOT_HISTS is a flag that causes figures of alpha vs data to be
% generated. 


  % Check inputs, set defaults
      
  if ~exist('max_amp','var') || isempty(max_amp)
    max_amp = 2;
  end
  
  if ~exist('max_del','var') || isempty(max_del)
    max_del = 1;
  end
  
  if ~exist('pseudocount','var') || isempty(pseudocount)
    pseudocount = 0;
  end
  
  if ~exist('plot_hists','var') || isempty(plot_hists)
    plot_hists = 0;
  end
  
  verbose('Computing alpha from data!',20);
  
  % Calculate histograms
  
  focalAmps = focals.amp + focals.aod;
  focalDels = focals.del + focals.doa;
 
  hxA = amp_thresh:(max_amp-amp_thresh)/50:max_amp;
  hxD = del_thresh:(max_del-del_thresh)/50:max_del;

  yy = histc(focalAmps,hxA);
  yy2 = histc(focalDels,hxD);
  
  hc = sum(yy,2);
  hc2 = sum(yy2,2);
  
  % Add pseudocounts (pseudocount % of total) 
  hc = hc+repmat(pseudocount/100*sum(hc)/length(hc),size(hc));
  hc2 = hc2+repmat(pseudocount/100*sum(hc2)/length(hc2),size(hc2));
  
  % Normalize distributions
  hc = hc/sum(hc);
  hc2 = hc2/sum(hc2);

  % Compute alpha
  [alpha(1) beta(1)] = compute_alpha(hc,hxA);
  [alpha(2) beta(2)] = compute_alpha(hc2,hxD);
  
  if any(isnan(alpha)) | any(isinf(alpha))
    error(['Alpha unbounded or not computable.  Should consider adding ' ...
           'pseudocounts or using different threshold parameters to make computation more robust!'])
  end

  if plot_hists
    figure()
    subplot(2,1,1)
    semilogy(hxA,hc);
    hold on
    line(hxA,beta(1)*exp(-alpha(1)*hxA),'Color','red');
    title('Amps')
    
    subplot(2,1,2)
    semilogy(hxD,hc2);
    hold on
    line(hxD,beta(2)*exp(-alpha(2)*hxD),'Color','red');
    title('Dels')
  end

  verbose(['alpha = [' num2str(alpha(1)) ' ' num2str(alpha(2)) ']'],30);
