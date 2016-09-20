function [throw_out,CL501]=histogram_qc(CL501,hist_fname,calc_smooth,showfig,show_extra,no_log,delta,sig)
% FIXME: don't popup figures

if ~exist('showfig','var')
    showfig=0;
end

if ~exist('show_extra','var')
    show_extra=1;
end

if ~exist('no_log','var')
    no_log=0;
end

if ~exist('sig','var')
    sig=0.05;
end

if ~exist('delta','var')
    delta=0.01;
end

if ischar(CL501.sdesc)
    CL501.sdesc=cellstr(CL501.sdesc);
end

if exist(hist_fname,'file')
    verbose(['Appending ' hist_fname],10);
end

if min(min(CL501.dat))< 0
    islogtransformed = 1;
else
    islogtransformed = 0;
end


if exist('calc_smooth','var') && calc_smooth
    CL=CL501;
    if ~islogtransformed % not log-transformed
        CL = preproc_log2trans(CL,1,1);
    else
       verbose('data are already log-transformed',10);
    end
    verbose('smoothing copy number',10);
    CL501=smooth_copy_number(CL,501,0,'mean'); % 'median');
end

%delta=0.01; % 0.001
%sig=0.05; % 0.05;
verbose(['delta=' num2str(delta) '; sig=' num2str(sig)],10);
if no_log
    delta=2^(delta+1)-2;
    sig=2^(sig+1)-2;
    %  CL501.smooth=2.^(CL501.smooth+1);
    CL501 = preproc_invlog2trans(CL501,'smooth');
end
vx=-5*sig:delta:5*sig;
v=exp(-(vx).^2/(2*sig.^2));
v=v./sum(v);
throw_out=[];

if showfig
    figure(1);
    set(gcf,'Visible','off');
    clf;
end

I=size(CL501.smooth,2);
J=1;
CL501=add_chrn(CL501);
r=find(CL501.chrn~=23); % no X chromosome;
for j=1:J %8
    for i=1:I %12
        if (i+(j-1)*I <= size(CL501.smooth,2))
            if no_log
                e=histc(CL501.smooth(r,i+(j-1)*I),0:delta:8);
                %      e(3006)=e(3005); % remove artificial peak at 0
                ec=conv(e,v);
                p=find_peaks(ec,delta); % delta was 0.01
                xc=(-5*sig:delta:(8+5*sig));
            else
                e=histc(CL501.smooth(r,i+(j-1)*I),-3:delta:3);
                %      e(3006)=e(3005); % remove artificial peak at 0
                ec=conv(e,v);
                p=find_peaks(ec,delta); % delta was 0.01
                xc=(-3-5*sig):delta:(3+5*sig);
            end
            if (showfig)
                clf;
                set(gcf,'Visible','off');
                %        subplot(I,J,i+(j-1)*I);
                if no_log
                    bh=bar(0:delta:8,e); hold on
                else
                    bh=bar(-3:delta:3,e); hold on
                end
                set(bh,'EdgeColor',[0.8 0.8 0.8],'FaceColor',[0.8 0.8 0.8])
                plot(xc(1:min(length(xc),length(ec))),ec(1:min(length(xc),length(ec)))); hold on
                ax=axis;
                if no_log
                    axis([0 4 ax(3:4)]);
                else
                    axis([(-3-5*sig) (3+5*sig) ax(3:4)]);
                end
                ax=axis;
                if show_extra
                    cah=gca;
                    pos=get(gca,'Position');
                    ah=axes('position',[pos(1)+0.01*pos(3) pos(2)+0.69*pos(4) pos(3)*0.4 pos(4)*0.3]);

                    plot(xc,ec);
                    set(gcf,'Visible','off');
                    axis([-0.6 0.6 ax(3:4)]);
                    set(gca,'XTick',[-0.3 0 0.3],'XTickLabel',num2str([-0.3;0;0.3]));
                    axes(cah);
                    for k=1:length(p)
                        line([xc(p(k)) xc(p(k))],[ax(3) ax(4)],'Color','r');
                    end
                end
                th=title(CL501.sdesc{i+(j-1)*I});

                set(th,'interpreter','none');
            end
            if length(p)==1
                if (showfig)
                    line([ax(1) ax(2)],[ax(4) ax(3)],'Color','k');
                end
                throw_out=[throw_out i+(j-1)*I];
                verbose(CL501.sdesc{i+(j-1)*I},10);
            end
            if showfig && show_extra
                axes(ah);
            end
        end
        if (showfig)
            print('-f1','-dpsc','-append',hist_fname); % hists.ps
        end
    end
end

%
%if (showfig)
%  print_D('histograms',{{'png','-r180'}},0);
%  % pause
%end
