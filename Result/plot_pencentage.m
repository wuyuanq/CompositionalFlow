

fh = figure();

subplot(1,2,1);
x=[1 2 4 5];
y=[8607.3 986.8 227.5; 15120.2 1390.7 355.8; 5.11 988.3 148.8; 9.1 1391.7 257.4];
bar(x,y,'stacked');
set(gca,'XTickLabel',{'480 no S','640 no S','480 with S','640 with S'},'fontsize',16)
legend('Flash time','Solver time','Other time');
ylabel('Run time(s)','fontsize',16);

subplot(1,2,2);
x=[1 2 4 5];
y=[87.6 10.0 2.4; 89.6 8.2 2.2; 0.45 86.5 13.05; 0.55 83.9 15.55];
bar(x,y,'stacked');
set(gca,'XTickLabel',{'480 no S','640 no S','480 with S','640 with S'},'fontsize',16)
ylabel('Percentage(%)','fontsize',16);

saveas(fh, 'percentage.fig');