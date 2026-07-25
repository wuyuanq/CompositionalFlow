
fh = figure();

x=[1 2 3 4];
y=[1 1; 1.82 1.86; 3.12 3.31; 4.81 5.40];
bar(x,y);
set(gca,'XTickLabel',{'64','128','256','512'},'fontsize',16)
xlabel('The number of processors','fontsize',16);
ylabel('Speedup','fontsize',16);
hold on

spu = [1,2,4,8];
plot(x, spu, 'k-*', 'linewidth', 2);
legend('480*480','640*640','ideal','location','northwest');

saveas(fh, 'bar.fig');