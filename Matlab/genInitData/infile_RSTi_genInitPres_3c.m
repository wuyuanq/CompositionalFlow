
path('/Users/wuy/Research/CompositionalFlow_matlab', path);
path('/Users/wuy/Research/CompositionalFlow_matlab/2D', path);

% Description of multicompnent fluid system: 
model.Nc = 3;
model.T = 220.0;      % K      

Lx = 0.1;    
Ly = 4;             % meters
timeEnd = 0.1*365*24*3600;  
m4sat = 2;  
 
nx = 1;  
ny = 40;  
nt = 0.1*365*24;
model.nx = nx;   
model.ny = ny; 
model.nt = nt;
model.xs = (0:nx)*Lx/nx; 
model.ys = (0:ny)*Ly/ny; 
model.ts = (0:nt)*timeEnd/nt; 

model.gravX = 0.0; 
model.gravY = -9.807; % unit: m/s^2

K_const = 9.869233*1.D-15;   
model.Kxx = zeros(nx,ny);  
model.Kxx(1:end, 1:end) = K_const; 
model.Kyy = model.Kxx; 

% Se = (Sw - Srw) / (1 - Srw - Srn).
% but for this case, Srw=0; Srn=0;  thus Se = Sw.
% k_rw = Se^m;  k_rn = (1-Se)^m;     % quadratic
model.kr_W = @(satW) (min(satW, 1) ).^m4sat;
model.kr_N = @(satW) (1-min(satW, 1)).^m4sat;

model.poro = zeros(nx,ny);  
model.poro(1:end, 1:end) = 0.2;

model.capP = 0.0;            % no capillary pressure
model.capPDeri = 0.0;        % not needed here  
   
model.src = zeros(model.Nc,nx,ny);   

model.isDiriX = zeros(2, ny);  
model.isDiriY = zeros(nx, 2);    

basePres = 2.0*1.D6;   % Pa 
model.PwBdryX = zeros(2, ny);
model.PwBdryY = zeros(nx, 2); 

model.zBdryX{1} = zeros(2, ny);
model.zBdryY{1} = zeros(nx, 2);
model.zBdryX{2} = zeros(2, ny);
model.zBdryY{2} = zeros(nx, 2);
model.UwBdryX = zeros(2, ny);
model.UwBdryY = zeros(nx, 2);
model.UnBdryX = zeros(2, ny);
model.UnBdryY = zeros(nx, 2);

model.PwInit = zeros(nx,ny);  
model.PwInit(1:end, 1:end) = basePres; 

model.zInit = zeros(model.Nc,nx,ny);   
model.zInit(2, 1:end, 1:end) = 0.3;
model.zInit(3, 1:end, 1:end) = 0.7;

model.ct = zeros(model.Nc, 1); % critical temperature
model.cp = zeros(model.Nc, 1); % critical pressure
model.af = zeros(model.Nc, 1); % acentric factor
model.mw = zeros(model.Nc, 1); % molar weight 
model.cv = zeros(model.Nc, 1); % critical volume
model.psatA = zeros(model.Nc, 1); % parameters to compute saturation pressures
model.psatB = zeros(model.Nc, 1);
model.psatC = zeros(model.Nc, 1);

% 1 is CO2, 2 is CH4 and 3 is C3H8
model.ct(1) = 304.4;
model.ct(2) = 190;
model.ct(3) = 370;

model.cp(1) = 7.4*1.D6;
model.cp(2) = 4.6*1.D6;
model.cp(3) = 4.2*1.D6;

model.af(1) = 0.23;
model.af(2) = 0.01;
model.af(3) = 0.15;

model.mw(1) = 0.044;
model.mw(2) = 0.016;
model.mw(3) = 0.044;

model.cv(1) = 0.0021;
model.cv(2) = 0.0062;
model.cv(3) = 0.0045;

model.psatA(1) = 6.81228;
model.psatA(2) = 6.69561;
model.psatA(3) = 6.82973;

model.psatB(1) = 1301.679;
model.psatB(2) = 405.420;
model.psatB(3) = 813.2;

model.psatC(1) = 269.506;
model.psatC(2) = 267.777;
model.psatC(3) = 248;

model.delta(:,:) = 0.0;
model.delta(1,2) = 0.15;
model.delta(2,1) = model.delta(1,2);
model.delta(1,3) = 0.1239;
model.delta(3,1) = model.delta(1,3);
model.delta(2,3) = 0.036;
model.delta(3,2) = model.delta(2,3);

model.soludoc = 'initData_3c';             

RST_compositionalTwoPhaseFlow(model)   % output to file soln_2PhFlw_Sw_RSTo.m etc 

rmpath('/Users/wuy/Research/CompositionalFlow_matlab');
rmpath('/Users/wuy/Research/CompositionalFlow_matlab/2D');

