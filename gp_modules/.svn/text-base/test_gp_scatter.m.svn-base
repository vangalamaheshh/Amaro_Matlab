close all
clear all
cd /xchip/projects/analytic_workshop/GaddyHO

D=read_mit_gct_file('all_aml_mll.mas5.mad300.gct');
D=read_mit_cls_file(D,'all_aml_mll.cls',1);
c=hsv2rgb([ repmat((0:(1/9):1)',3,1) [ ones(20,1); 0.5*ones(10,1)] ...
            [ones(10,1); 0.5*ones(10,1); ones(10,1);]]);
c=c(mod((1:7:(size(c,1)*7))-1,size(c,1))+1,:);
D=add_supmark(D,c);
D.supmark(1).marker=cellstr(repmat('o',size(D.supmark(1).marker,1),1));

close all
figure(1); clf;
scatter_plot_D(D.dat(1:3,:)',D,1,1)
grid on

[typeacc,typedesc,D1,range,non_empty]=decollapse_supdat(D,1);
legend(deblank(cellstr(D1.supacc(range,:))));

HOpath='/xchip/projects/analytic_workshop/GaddyHO/';
gp_dim_red([ HOpath 'all_aml_mll.mas5.mad300.gct'],...
           [ HOpath 'all_aml_mll.mas5.mad300.pca3.gct'],...
           'PCA','3','cols')
gp_scatter_plot([ HOpath 'all_aml_mll.mas5.mad300.pca3.gct'],...
                '1','1:3',[ HOpath 'all_aml_mll.cls']);
