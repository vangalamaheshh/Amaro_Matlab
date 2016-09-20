function G = temp_mutsig2_gather(dirname,write_results_mat_flag)

if ~exist('write_results_mat_flag','var'), write_results_mat_flag=true; end

fprintf('Gathering output...');
system(['cd ' dirname ';grep REPORT sc*/*.out > results.txt']);
fprintf('\n');

G=[]; G.line = load_lines([dirname '/results.txt']);
if isempty(G.line), fprintf('NO RESULTS!\n'); return; end

G = parse_in(G,'line','REPORT\]?\s+(\S+)\s+(.*)$',{'gene','details'});


if ~isempty(grep('eff',G.details)) % report_type==2 new version
%scatter.0000000001/job.764644.out:[REPORT]  TP53        nmuts 58   len 1311  nperm 1000     (2722 perm/sec)  CONS k 0     p 0.000000 eff 1.85  CLUST k 0     p 0.000000 eff 6.23  JOINT k 0     p 0.000000

  G = parse_in(G,'details',['nmuts\s+(\d+)\s+len\s+(\d+)\s+nperm\s+(\d+).*'...
                       '\s+CONS\s+k\s+(\d+)\s+p\s+(\S+)\s+eff\s+(\S+)'...
                      '\s+CLUST\s+k\s+(\d+)\s+p\s+(\S+)\s+eff\s+(\S+)'...
                      '\s+JOINT\s+k\s+(\d+)\s+p\s+(\S+)$'],...
               {'nmuts','len','nperm',...
                'kcons','pcons','effcons',...
                'kclust','pclust','effclust',...
                'kjoint','pjoint'},1:11);

elseif ~isempty(grep('me/ma/mo/mx',G.details))   % report_type==2 old version
%  G = parse_in(G,'line',['REPORT\]?\s+(\S+)\s+nmuts\s+(\d+)\s+len\s+(\d+)\s+nperm\s+(\d+)'...
%                      '\s+CONS\s+k\s+(\d+)\s+.*'...
%                      '\s+CLUST\s+k\s+(\d+)\s+.*'...
%                      '\s+JOINT\s+k\s+(\d+)\s+.*'],...
%               {'gene','nmuts','len','nperm','kcons','kclust','kjoint'},2:7);

%scatter.0000000001/job.764644.out:[REPORT]  AHNAK       nmuts 15   len 17801 nperm 1000     CONS k 898   p 0.90 e (me/ma/mo/mx) -0.56/-0.57/-0.01/-2.07  CLUST k 6     p 0.01 e 0.20/0.19/0.26/-0.10  JOINT k 34    p 0.03 e -0.37/-0.37/1.02/-2.00

  G = parse_in(G,'details',['nmuts\s+(\d+)\s+len\s+(\d+)\s+nperm\s+(\d+).*'...
                      '\s+CONS\s+k\s+(\d+)\s+p\s+(\S+)\s+e\s+\(me/ma/mo/mx\)\s+(\S+)/(\S+)/(\S+)/(\S+)'...
                      '\s+CLUST\s+k\s+(\d+)\s+p\s+(\S+)\s+e\s+(\S+)/(\S+)/(\S+)/(\S+)'...
                      '\s+JOINT\s+k\s+(\d+)\s+p\s+(\S+)\s+e\s+(\S+)/(\S+)/(\S+)/(\S+)$'],...
               {'nmuts','len','nperm',...
                'kcons','pcons','econs_me','econs_ma','econs_mo','econs_mx',...
                'kclust','pclust','eclust_me','eclust_ma','eclust_mo','eclust_mx',...
                'kjoint','pjoint','ejoint_me','ejoint_ma','ejoint_mo','ejoint_mx'},1:21);
elseif ~isempty(grep('kcons',G.details))
  G = parse_in(G,'details',['nmuts\s+(\d+)\s+len\s+(\d+)\s+kcons\s+(\d+)'...
                      '\s+kclust\s+(\d+)\s+kjoint\s+(\d+)\s+nperm\s+(\d+)\s+'],...
               {'nmuts','len','kcons','kclust','kjoint','nperm'},1:6);
else
  G = parse_in(G,'details','(\S+) (\S+) (\S+) (\S+) (\S+)$',{'kcons','kclust','kjoint','nperm'},1:4);
  G.kcons=G.kcons-1; G.kclust=G.kclust-1;  % (old versions needed this)
end

G = rmfield(G,{'line','details'});

G = sort_struct(G,'nperm');
[u ui uj] = unique(G.gene,'last');
G = reorder_struct(G,ui);

flds={'cons','clust','joint'}; for i=1:length(flds), fld=flds{i};
  [G.(['p' fld]) G.(['ci' fld])] = calc_pval_and_ci_ratio(G.(['k' fld]),G.nperm);
end

G = sort_struct(G,'pjoint');


if write_results_mat_flag
  save([dirname '/results.mat'],'G');
end



