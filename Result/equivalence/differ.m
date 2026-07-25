
fh = figure();

x = 1:0.2*365*24*4;

y1 = load('40/soln_cmp2PhFlw_aveDiff1.txt');
y2 = load('40/soln_cmp2PhFlw_aveDiff2.txt');
plot(x, y1, 'k-', 'linewidth', 2);
hold on
plot(x, y2, 'k--', 'linewidth', 2);
hold on

y1 = load('80/soln_cmp2PhFlw_aveDiff1.txt');
y2 = load('80/soln_cmp2PhFlw_aveDiff2.txt');
plot(x, y1, 'r-', 'linewidth', 2);
hold on
plot(x, y2, 'r--', 'linewidth', 2);
hold on

y1 = load('120/soln_cmp2PhFlw_aveDiff1.txt');
y2 = load('120/soln_cmp2PhFlw_aveDiff2.txt');
plot(x, y1, 'g-', 'linewidth', 2);
hold on
plot(x, y2, 'g--', 'linewidth', 2);
hold on

y1 = load('160/soln_cmp2PhFlw_aveDiff1.txt');
y2 = load('160/soln_cmp2PhFlw_aveDiff2.txt');
plot(x, y1, 'b-', 'linewidth', 2);
hold on
plot(x, y2, 'b--', 'linewidth', 2);
hold on

y1 = load('200/soln_cmp2PhFlw_aveDiff1.txt');
y2 = load('200/soln_cmp2PhFlw_aveDiff2.txt');
plot(x, y1, 'c-', 'linewidth', 2);
hold on
plot(x, y2, 'c--', 'linewidth', 2);
hold on

xlabel('Time steps','fontsize',16);
ylabel('L1-norm','fontsize',16);
legend('40*40 1','40*40 2', ...
'80*80 1','80*80 2', ...
'120*120 1','120*120 2', ...
'160*160 1','160*160 2', ...
'200*200 1','200*200 2','location','northeast');

saveas(fh, 'diff.fig');
