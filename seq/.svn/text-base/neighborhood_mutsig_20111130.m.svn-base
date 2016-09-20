function G = neighborhood_mutsig(G,V,P)
% based on /xchip/cga1/lawrence/mut/analysis/20110909_pancan/run.m

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'vidx',1:slength(V));
P = impose_default_value(P,'genes_to_analyze',G.name);
P = impose_default_value(P,'use_categories',false);
P = impose_default_value(P,'max_neighbors',50);
P = impose_default_value(P,'theta_max',1.5);
P = impose_default_value(P,'theta_halt',0.5);
P = impose_default_value(P,'qual_min',0.2);
P = impose_default_value(P,'incl_silent',true);
P = impose_default_value(P,'incl_flank',false);
P = impose_default_value(P,'scale_for_negative_selection',true);
P = impose_default_value(P,'pval_calc_method','weighted',{'weighted','weighted_finely','hyge2','binomial'});
P = impose_default_value(P,'p_mention_threshold',0.0001);
P = impose_default_value(P,'display_mode','summary',{'detail','mention','summary','final'});
P = impose_default_value(P,'show_texthist',false);

% check input data
ng=slength(G);nv=slength(V);
if ~(ng>=1 && nv>=1), error('problem with G or V'); end

demand_fields(G,{'name'});
flds = {'nnon','Nnon'};
if P.incl_silent, flds = [flds 'nsil' 'Nsil']; end
if P.incl_flank, flds = [flds 'nflank' 'Nflank']; end
demand_fields(G,flds);
sz = size(getfield(G,flds{1}));
if sz(1)~=ng, error('problem with G'); end
if length(sz)>2, error('3D matrix not appropriate for n and N fields'); end
for i=2:length(flds)
  szi = size(getfield(G,flds{i}));
  if length(sz)~=length(szi) || ~all(szi==sz)
    error('problem with n and N fields in G struct');
  end
end
if sz(2)>1
  if ~P.use_categories
    fprintf('Using total category only\n');
    P = impose_default_value(P,'total_category','*required*');
    for i=1:length(flds)
      tmp = getfield(G,flds{i});
      tmp = tmp(:,P.total_category);
      G = setfield(G,flds{i},tmp);
    end
    pcat = [];
    npcat = [];
    totcat = P.total_category;
  else
    fprintf('Using categories\n');
    P = impose_default_value(P,'point_categories','*required*');
    P = impose_default_value(P,'non_point_categories','*required*');
    P = impose_default_value(P,'total_category','*required*');
    pcat = P.point_categories;
    npcat = P.non_point_categories;
    totcat = P.total_category;
  end
else
  if P.use_categories
    error('No category information provided (i.e. n and N columns)');
  else
    fprintf('No categories\n');
    pcat = [];
    npcat = [];
    totcat = 1;
  end
end

vidx = P.vidx;
gta = P.genes_to_analyze;
max_neighbors = P.max_neighbors;
theta_max = P.theta_max;
theta_halt = P.theta_halt;
qual_min = P.qual_min;
incl_silent = P.incl_silent;
preport = P.p_mention_threshold;

% convert covariate raw values to Z-scores
ng=slength(G);nv=slength(V);
Z = nan(ng,nv);
for vi=1:nv
  if length(V.val{vi})~=ng, error('V.val must match genelist in G'); end
  V.missing{vi,1} = (isnan(V.val{vi}) | isinf(V.val{vi}));
  mn = mean(V.val{vi}(~V.missing{vi}));
  sd = std(V.val{vi}(~V.missing{vi}));
  Z(:,vi) = (V.val{vi}-mn)./sd;
  if ~isempty(V.missing{vi})
    Z(V.missing{vi},vi) = 0;
  end
end

globalrate_non = sum(G.nnon)/sum(G.Nnon); G.Fnon = (G.nnon./G.Nnon)/globalrate_non;
if P.incl_silent, globalrate_sil = sum(G.nsil)/sum(G.Nsil); G.Fsil = (G.nsil./G.Nsil)/globalrate_sil; end
if P.incl_flank, globalrate_flank = sum(G.nflank)/sum(G.Nflank); G.Fflank = (G.nflank./G.Nflank)/globalrate_flank; end

Go=G;
if ischar(gta), gta={gta}; end
if islogical(gta), gta=find(gta); end
if isnumeric(gta)
  idx = gta;
  if any(idx)<1 || any(idx)>ng, error('P.genes_to_analyze out of bounds'); end
elseif iscellstr(gta)
  idx = listmap(gta,Go.name);
else
  error('unknown format for P.genes_to_analyze');
end
G = reorder_struct(Go,idx(~isnan(idx)));

header = sprintf('%-3s %-8s %-4s %-9s %-4s %-9s','ng','neighb','nsil','Nsil','nnon','Nnon');
if P.incl_flank, header = [header sprintf(' %-4s %-9s','nflk','Nflk')]; end
header = [header sprintf(' %-5s %-5s %-4s %-9s %-5s %-5s %-5s','dist2','rms','nfit','Nfit','Fmle','Fthe','qual')];
if P.show_texthist, header = [header sprintf(' %-3s %-3s %-3s %-3s %-3s %-3s %-3s %-3s %-3s','0.1','0.5','1.0','2.0','3.0','5.0','10','20','50')]; end
header = [header sprintf(' %-7s %-7s','p','bon')];

if strcmpi(P.display_mode,'mention'), fprintf('%-10s  %s\n','gene',header); end

z=nan(slength(G),1); G.nnei=z;G.nfit=z;G.Nfit=z;G.Fmle=z;G.Fthe=z;G.qual=z;G.p_cat=z;G.bon=z;G.effect_cat=z;

for gi=1:slength(G), gio = find(strcmp(Go.name,G.name{gi})); c=1; Fobs = G.Fnon(gi,c);

  if strcmpi(P.display_mode,'detail')
    fprintf(['\n__( %d / %d)__[ %s ]' repmat('_',1,46) '[ Fobs = %4.2f ]' repmat('_',1,58) '\n'],gi,slength(G),G.name{gi},Fobs);
    fprintf('%s\n',header);
  end

  % calculate distances from this gene
  dist2 = sum(bsxfun(@minus,Z(:,vidx),Z(gio,vidx)).^2,2);
  [tmp,ord] = sort(dist2); ord= [gio;ord(ord~=gio)];

  % for "total" category and for each individual point-mutation category (if using categories)
  for cati = 1:length(pcat)+1
    if cati==1
      c = totcat;
    else
      c = pcat(cati-1);
    end

    % for each neighborhood
    ntot=0; Ntot=0; sumdist2=0; old_theta=inf;
    for ni=1:max_neighbors, gidx = ord(ni); sumdist2=sumdist2+dist2(gidx);

      % add up mutations
      for muttype=1:3
        if muttype==1   % silemt
          if ~incl_silent,continue;end;
          n=Go.nsil(gidx,c); N=Go.Nsil(gidx,c);
          if P.scale_for_negative_selection
            scale = globalrate_non/globalrate_sil;
            %        n=n*scale;
            N=N/scale;
          end
        elseif muttype==2  % nonsilent
          if ni==1,continue;end;
          n=Go.nnon(gidx,c);N=Go.Nnon(gidx,c);
        elseif muttype==3   % flanking
          if ~P.incl_flank,continue;end;
          n=Go.nflank(gidx,c);N=Go.Nflank(gidx,c);
          if P.scale_for_negative_selection
            scale = globalrate_non/globalrate_flank;
            %        n=n*scale;
            N=N/scale;
          end
        end
        if (N>n),ntot=ntot+n;Ntot=Ntot+N; end
      end

      % statistics of neighborhood
      [mle ci] = binofit(round(ntot),round(Ntot)); theta=(ci(2)-ci(1))/mle; rms=sqrt(sumdist2/ni);
      if ntot==0 && Ntot==0, nfit=sum(Go.nnon(:,c)); Nfit=sum(Go.Nnon(:,c)); ff2sz=1000; else nfit=ntot; Nfit=Ntot; ff2sz=20; end
      qual = fisher_exact_test(Go.nsil(gio,c),Go.Nsil(gio,c),ntot,Ntot);

      % stopping criterion #1 (exclude this neighbor)
      %    stop if theta is already below theta_max and quality has just dropped below qual_min
      if ni>1 && (incl_silent==true && old_theta<=theta_max && qual<qual_min)
        if cati==1  % print final outcome (only if we're working on the "total" category)
          if strcmpi(P.display_mode,'mention') && G.p_cat(gi,c)<preport
            fprintf('%-10s  %s',G.name{gi},old_str);
          elseif strcmpi(P.display_mode,'summary') && G.p_cat(gi,c)<preport
            fprintf('\n'); pr(sort_struct(reorder_struct(G,G.p_cat(:,c)<preport),{'p','effect'},[1 -1]));
          end
        end
        break
      end

      % calculate p-value
      if P.show_texthist
        ff=[0.1 0.5 1 2 3 5 10 20 50]; lk = binopdf(round(nfit),round(Nfit),globalrate_non*ff);
      end
      nobs = G.nnon(gi,c); Nobs = G.Nnon(gi,c);
      if strcmp(P.pval_calc_method,'weighted_finely')
        ff2sz=ff2sz*10;
      end
      if strcmp(P.pval_calc_method,'weighted') || strcmp(P.pval_calc_method,'weighted_finely')
        ff2=geoseries(0.05,100,ff2sz); lk2 = binopdf(round(nfit),round(Nfit),globalrate_non*ff2);
        pdist = 1-binocdf(nobs-1,Nobs,globalrate_non*ff2); pval = sum(pdist.*(lk2/sum(lk2)));
      elseif strcmp(P.pval_calc_method,'hyge2')
        pval = 1-hyge2cdf(nobs-1,Nobs,nfit,Nfit);
      elseif strcmp(P.pval_calc_method,'binomial')
        pval = 1-binocdf(nobs-1,Nobs,nfit/Nfit);
        %      pval = 1-binocdf(nobs-1,Nobs,1.0e-5);
      end
      if Nobs==0, pval=1; end
      if pval>1, pval=1; end
      if pval<0, pval=0; end
      
      % Bonferroni correction
      bon = min(1,ng*pval);

      if cati==1  % (print stats if we're working on the "total" category)
        str = sprintf('%-3.0f %-8s %-4.0f %-9.0f %-4.0f %-9.0f',...
                      ni,Go.name{gidx}(1:min(8,length(Go.name{gidx}))),Go.nsil(gidx),Go.Nsil(gidx),Go.nnon(gidx),Go.Nnon(gidx));
        if P.incl_flank, str = [str sprintf(' %-4.0f %-9.0f',Go.nflank(gidx),Go.Nflank(gidx))]; end
        str = [str sprintf(' %-5.2f %-5.2f %-4.0f %-9.0f %-5.2f %-5.2f %-5.2f',...
                           dist2(gidx),rms,nfit,Nfit,mle/globalrate_non,theta,qual)];
        if P.show_texthist, str = [str sprintf([' %-3.0f %-3.0f %-3.0f %-3.0f %-3.0f %-3.0f %-3.0f %-3.0f %-3.0f'],-log10(lk))]; end
        str = [str sprintf(' %-7s %-7s\n',format_number(pval,1,7),format_number(bon,1,7))];
        if strcmpi(P.display_mode,'detail'), fprintf('%s',str); end
      end
     
      % update gene's statistics
      G.nnei(gi,c)=ni-1; G.nfit(gi,c) = round(ntot); G.Nfit(gi,c) = round(Ntot);
      G.Fmle(gi,c)=mle/globalrate_non; G.Fthe(gi,c)=theta; G.qual(gi,c)=qual; G.p_cat(gi,c)=pval; G.bon(gi,c)=bon;
      G.effect_cat(gi,c) = G.Fnon(gi,c)./G.Fmle(gi,c);
      old_str=str;old_theta=theta;

      % stopping criterion #2 (include this neighbor)
      %     stop if we have achieved a theta as low as theta_halt
      %     or if we've reached max_neighbors
      if (ni==max_neighbors || theta<=theta_halt)
        if cati==1  % print final outcome (only if we're working on the "total" category)
          if strcmpi(P.display_mode,'mention') && G.p_cat(gi,c)<preport
            fprintf('%-10s  %s',G.name{gi},str);
          elseif strcmpi(P.display_mode,'summary') && G.p_cat(gi,c)<preport
            fprintf('\n'); pr(sort_struct(reorder_struct(G,G.p_cat(:,c)<preport),{'p','effect'},[1 -1]));
          end
        end
        break
      end
    end % next neighborhood size
  end % next category
end % next gene

c = totcat;
G.p = G.p_cat(:,c);
G.q = calc_fdr_value(G.p);
G.effect = G.effect_cat(:,c);

G=sort_struct(G,{'p','effect'},[1 -1]);
G.rank=(1:slength(G))';

if strcmpi(P.display_mode,'final')
  tmp = G;
  G = rmfield_if_exist(G,{'nobs','Nobs','nintron','Nintron'});
  pr(G,1:min(slength(G),40))
end
