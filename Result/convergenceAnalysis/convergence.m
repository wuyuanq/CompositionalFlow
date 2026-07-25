
timestep = 876;

% load data
base = load('640/soln_cmp2PhFlw_moleHistory.txt');

error = zeros(5,1);
x = zeros(5,1);
for n = 1 : 5
    filename = [num2str(40*(n+1)),'/soln_cmp2PhFlw_moleHistory.txt'];
    data = load(filename);

    for k = 2 : timestep 
        error(n) = error(n) + abs(base(k)-data(k))/abs(base(k));
    end
    error(n) = error(n)/timestep;
    x(n) = (40*(n+1))^2;
end
error
fh = figure();

loglog(x,error,'-s')
xlabel('Grid Size','fontsize',16);
ylabel('Error','fontsize',16);

saveas(fh, 'convergence.fig');




