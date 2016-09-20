function fh_test1(bamname,outname)

fprintf(['bamname = ' bamname '\n']);

save_textfile(['<h3>Test Output</h3><p>bamname = ' bamname '<hr><p><img src="plot.png"><hr>'],'report.html');

save_lines({'Testing.';'';'bamname:';['    ' bamname];'';'End test.'},outname);

a = rand(200,3);
figure(1);
scatter(a(:,1),a(:,2),[],a,'o','filled');
title('Results Plot','fontsize',20);
text(0.5,0.5,bamname,'horizontalalign','center','interpreter','none');
xlabel('Velocity (furlongs/fornight)');
ylabel('Luminosity (Norwegian standard units)');
print('-dpng','-r150','plot.png');
print('-dpng','-r300','plot300.png');
%close(1);

%pause(2);

close all   % <-------------- crucial!  allows xvnc to close!


