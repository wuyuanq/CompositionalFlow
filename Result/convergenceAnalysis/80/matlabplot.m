path('/Users/wuy/Dropbox/Research/CompositionalFlow_matlab/2D', path);
model.Nc =  2;
Lx =     4.00;
Ly =     4.00;
timeEnd =   315360.0;
m4sat = 2;
nx =    80;
ny =    80;
nt =        876;
model.nx = nx;
model.ny = ny;
model.nt = nt;
model.xs = (0:nx)*Lx/nx;
model.ys = (0:ny)*Ly/ny;
model.ts = (0:nt)*timeEnd/nt;
fpwtxt = 'soln_cmp2PhFlw_Pw_raw.txt';
fuxtxt = 'soln_cmp2PhFlw_Ux_raw.txt';
fuytxt = 'soln_cmp2PhFlw_Uy_raw.txt';
fmftxt = [];
for m = 1 : model.Nc
    fk = ['soln_cmp2PhFlw_x', num2str(m), '_raw.txt'];
    fmftxt = [fmftxt; fk];
end
fswtxt = 'soln_cmp2PhFlw_Sw_raw.txt';
fxitxt = 'soln_cmp2PhFlw_c_raw.txt';
fmhtxt = 'soln_cmp2PhFlw_moleHistory.txt';
fmrtxt = 'soln_cmp2PhFlw_moleRatioHistory.txt';
ftimetxt = 'soln_cmp2PhFlw_time.txt';
model.soludoc = 'matlabplots';
if(exist(model.soludoc,'dir') ~= 7)
    system(['mkdir ' model.soludoc]);
end
RST_plot( model, fpwtxt, fuxtxt, fuytxt, fmftxt, fswtxt, fxitxt, fmhtxt, fmrtxt);
rmpath('/Users/wuy/Dropbox/Research/CompositionalFlow_matlab/2D');
