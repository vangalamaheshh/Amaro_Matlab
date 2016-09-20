function [pwr,app_pwr]=freq_test_power(p,alpha,p0,n,tail);
% [pwr,app_pwr]=freq_test_power(p,alpha,p0,n,tail);
%  calculate the power of rejecting the null hypothesis H0:P=p0
%  given n trials and assuming the true P is p

if ~exist('tail','var')
  tail=-1;
end
% pwr=freq_test_power(0.4,0.01,0.5,1:500);
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TESTS
%
if nargin==1
  switch p
   case 1
    N=50;
    p=0.2;
    
    x=0;
    C=5;
    for i=0:C
      x=x+binopdf(i,N,p);
    end
    x
    
    n1=2*(C+1);
    n2=2*(N-(C+1)+1);
    1-fcdf(n2*p/(n1*(1-p)),n1,n2)
    
    binocdf(C,N,p)

   case 2
    tic
    p=0.4; % 5/61 True (EGFR and Rb)
%    p=0.08; % 5/61 True (EGFR and Rb)
    p0=0.262; % 39/61 (EGFR) * 25/61 (Rb);
    alpha=0.01;
    n=1:500;
    [pwr,app_pwr]=freq_test_power(p,alpha,p0,n,1);
    toc
    keyboard
    figure(1); clf;
    plot(pwr,'b-'); hold on
    plot(app_pwr,'r-');
    set(gca,'FontSize',10);
    title(['Power(Sample size) for p0=' num2str(p0) ' (null hyp.); \alpha=' ...
           num2str(alpha) '; p=' num2str(p) ' (true)'],'FontSize',18);
    xlabel('N','FontSize',13);
    ylabel('Power(N)','FontSize',13);
    legend({'Exact','Gaussian Approx.'},'Location','SouthEast');
    print_pdf(['power_vs_N_' num2str(p0) '_' num2str(alpha) '_' ...
               num2str(p)]);
    
   case 3
    N=[50 100 200 400 800];
    f=0.1:0.1:0.9;
    lambda=[0.1 0.2 0.4 0.8];
    alpha=0.01;
    pwr=zeros(length(f),length(f),length(N),length(lambda));
    for li=1:length(lambda)
      for ni=1:length(N)
        for fi=1:length(f)
          for fj=1:fi
            pwr(fi,fj,ni,li)=power(lambda(li)*f(fi)*f(fj),alpha,N(ni),f(fi)*f(fj));
          end
        end
        pwr(:,:,ni,li)=pwr(:,:,ni,li)+triu(pwr(:,:,ni,li)',1);
      end
    end
    
    for ni=1:length(N)
      figure(1); clf;
      for li=1:length(lambda)
        subplot(2,2,li);
%        meshc(pwr(:,:,ni,li));
        imagesc(pwr(:,:,ni,li));      
        caxis([0 1]);
        colorbar;
        title(['Power(rejecting independence at \alpha=' num2str(alpha) ...
               ') ' 10 'for N=' num2str(N(ni)) '; \lambda=' ...
               num2str(lambda(li))]);
        xlabel('P_a');
        ylabel('P_b');
        set(gca,'XTick',1:length(f),'XTickLabel',f);
        set(gca,'YTick',1:length(f),'YTickLabel',f);
      end
      print_pdf(['power_' num2str(N(ni))]);
      disp('hit a key');
    end
    approx_pwr=[];
    return
    
   case 4
    p0=0.16;
    p=0.016;
    alpha=0.01;
    power_thresh=0.8;
    x0=fzero(@(x) norminv(power_thresh,exp(x)*p,sqrt(exp(x)*p*(1- ...
                                                      p)))-norminv(alpha,exp(x)*p0,sqrt(exp(x)*p0*(1-p0))),log(50));
    x=fzero(@(x) find_N(p,alpha,x,p0,power_thresh),x0);
    power(p,alpha,floor(exp(x)),p0)
    keyboard
   
   case 5
    f=0.1:0.1:0.9;
    lambda=[0.1 0.2 0.4 0.8];
    alpha=0.01;
    power_thresh=0.8;
    pwr_n=zeros(length(f),length(f),length(lambda));
    for li=1:length(lambda)
      for fi=1:length(f)
        for fj=1:fi
          p0=f(fi)*f(fj);
          p=lambda(li)*p0;
          pwr_n(fi,fj,li)=sample_size(p,alpha,p0,power_thresh);
        end
      end
      pwr_n(:,:,li)=pwr_n(:,:,li)+triu(pwr_n(:,:,li)',1);
    end
    pwr=pwr_n;
    keyboard
    
    figure(1); clf;
    for li=1:length(lambda)
      subplot(2,2,li);
      %        meshc(pwr(:,:,ni,li));
      x=pwr(:,:,li);
      x(x<1)=1;
      imagesc(log10(x));
      caxis([0.6 4]);
      clh=colorbar;
      set(clh,'YTickLabel',num2str(floor(10.^(get(clh,'YTick'))')));
      title(['N for rejecting independence at \alpha=' num2str(alpha) ...
             'with power>0.8 ' 10 '\lambda=' ...
             num2str(lambda(li))]);
      xlabel('P_a');
      ylabel('P_b');
      set(gca,'XTick',1:length(f),'XTickLabel',f);
      set(gca,'YTick',1:length(f),'YTickLabel',f);
    end
    print_pdf(['N_pwr_0.8']);
    approx_pwr=[];
    return

   case 7
    
  end
end

for i=1:length(n)
  ni=n(i);
  pwr(i)=power(p,alpha,ni,p0,tail);
  app_pwr(i)=approx_power(p,alpha,ni,p0,tail);
end

function n=sample_size(p,alpha,p0,power_thresh)
x0=fzero(@(x) norminv(power_thresh,exp(x)*p,...
                      sqrt(exp(x)*p*(1-p)))-...
         norminv(alpha,exp(x)*p0,...
                 sqrt(exp(x)*p0*(1-p0))),log(50));
x=fzero(@(x) find_N(p,alpha,x,p0,power_thresh),x0);
n=floor(exp(x));

function res=find_N(p,alpha,lgN,p0,power_thresh)
res=power(p,alpha,floor(exp(lgN)),p0)-power_thresh;

function pwr=power(p,alpha,N,p0,tail)
if ~exist('tail','var')
  tail=-1;
end
if tail <0
  CL=binoinv(alpha,N,p0); % return the smallest integer, x, for which binocdf(x,N,p0)>=alpha 
                          % since we want the largest integer for which
                          % binocdf(x,N,p0)<=alpha we added the following
  if binocdf(CL,N,p0)>alpha
    CL=CL-1;
  end
  if CL<0
    pwr=0;
  else
    pwr=binocdf(CL,N,p);
  end
elseif tail>0
%  keyboard
  CU=binoinv(1-alpha,N,p0); % return the smallest integer, x, for which binocdf(x,N,p0)>=1-alpha 
                            % since we want the largest integer for which
                            % binocdf(x,N,p0)<=alpha we add
  if binocdf(CU,N,p0)<1-alpha
    CU=CU+1;
  end
  if CU>N
    pwr=0;
  else
    pwr=1-binocdf(CU,N,p);
  end
else % two sided
  CL=binoinv(alpha/2,N,p0); % return the smallest integer, x, for which binocdf(x,N,p0)>=alpha 
                          % since we want the largest integer for which
                          % binocdf(x,N,p0)<=alpha we added the following
  if binocdf(CL,N,p0)>alpha/2
    CL=CL-1;
  end
  if CL<0
    pwr_left=0;
  else
    pwr_left=binocdf(CL,N,p);
  end
  
  CU=binoinv(1-alpha/2,N,p0); % return the smallest integer, x, for which binocdf(x,N,p0)>=1-alpha 
                            % since we want the largest integer for which
                            % binocdf(x,N,p0)<=alpha we add
  if binocdf(CU,N,p0)<1-alpha/2
    CU=CU+1;
  end
  if CU>N
    pwr_right=0;
  else
    pwr_right=1-binocdf(CU,N,p);
  end
  
  pwr=pwr_left+pwr_right;
end



function pwr=approx_power(p,alpha,N,p0,tail)
if ~exist('tail','var')
  tail=-1;
end
if tail<0
  CL=norminv(alpha,N*p0,sqrt(N*p0*(1-p0)));
  if CL<0
    pwr=0;
  else
    pwr=normcdf(CL,N*p,sqrt(N*p*(1-p)));
  end
elseif tail>0
  CU=norminv(1-alpha,N*p0,sqrt(N*p0*(1-p0)));
  if CU>N
    pwr=0;
  else
    pwr=1-normcdf(CU,N*p,sqrt(N*p*(1-p)));
  end
else
  CL=norminv(alpha/2,N*p0,sqrt(N*p0*(1-p0)));
  if CL<0
    pwr_left=0;
  else
    pwr_left=normcdf(CL,N*p,sqrt(N*p*(1-p)));
  end
  CU=norminv(1-alpha/2,N*p0,sqrt(N*p0*(1-p0)));
  if CU>N
    pwr_right=0;
  else
    pwr_right=1-normcdf(CU,N*p,sqrt(N*p*(1-p)));
  end
  pwr=pwr_left+pwr_right;
end

function x=binocdf1(C,N,p)
n1=2*(C+1);
n2=2*(N-(C+1)+1);
x=1-fcdf(n2*p/(n1*(1-p)),n1,n2);
