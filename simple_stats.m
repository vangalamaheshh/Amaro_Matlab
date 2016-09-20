function r=simple_stats(x)

x=as_column(x);
r.min=min(x);
r.max=max(x);
r.mean=mean(x);
r.std=std(x);
r.mad=mad(x);
r.q=quantile(x,[ .025 0.1 .25 .50 .75 0.9 .975]);

disp(['[ min mean max]=' num2str([r.min r.mean r.max],'%f ')]);
disp(['[ std mad]=' num2str([r.std r.mad],'%f ')]);
disp(['[ quantiles:  .025 0.1 .25 .50 .75 0.9 .975]=' num2str(r.q,'%f ')]);;

