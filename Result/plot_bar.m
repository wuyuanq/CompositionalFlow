
fh = figure();

subplot(1,2,1);
x=[1 2 4 5 7 8 10 11];
y=[19000 4.9 380.1; 9.6 4.9 380.5; 9850 3.8 285.2; 5.6 3.8 285.6; 5200 3.3 245.7; 2.8 3.3 245.9; 2756 3.0 233; 1.4 3.0 232.6];
bar(x,y,'stacked');
set(gca,'XTickLabel',{'128','','256','','512','','1024',''},'fontsize',16)
legend('Flash time','Solver time','Commnication and other time');
xlabel('The number of processors','fontsize',16);
ylabel('Time(s)','fontsize',16);

subplot(1,2,2);
x=[1 2 4 5 7 8 10 11];
y=[33200 6.9 662.7; 16.4 6.9 662.7; 17000 4.8 495.8; 8.4 4.8 495.8; 8890 3.8 417.8; 4.4 3.8 417.8; 4600 3.4 391.3; 2.3 3.4 391.3];
bar(x,y,'stacked');
set(gca,'XTickLabel',{'128','','256','','512','','1024',''},'fontsize',16)
legend('Flash time','Solver time','Commnication and other time');
xlabel('The number of processors','fontsize',16);
ylabel('Time(s)','fontsize',16);

saveas(fh, 'bar.fig');