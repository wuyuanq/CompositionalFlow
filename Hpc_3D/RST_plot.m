% This is the RST_plot() function which uses the model in the input file 
% and the solution files output by the function RST_compositionalFlow() as 
% input, and generates a series of images to demonstrate the results.

% Input parameters:
% model: the model
% fpwtxt: the pressure file in the txt form
% fuxtxt: the x-dirction velocity file in the txt form
% fuytxt: the y-dirction velocity file in the txt form
% fuztxt: the z-dirction velocity file in the txt form
% fmftxt: the mole fraction files in the txt form
% fswtxt: the saturation file in the txt form
% fxitxt: the mole density file in the txt form
% fmhtxt: the mole fraction history file of component 1 in the txt form
% fmrtxt: the mole ratio history file of desired components to component 1

% Author: Yuanqing Wu. Email: wuyuanq@gmail.com
% Last edited on October 23rd, 2013

function RST_plot( model, fpwtxt, fuxtxt, fuytxt, fuztxt, fmftxt, fswtxt, fxitxt, fmhtxt, fmrtxt )

    % simplify the symbols of the model
    Nc = model.Nc;
    nx = model.nx;
    ny = model.ny;
    nz = model.nz;
    xs = model.xs;
    ys = model.ys;
    zs = model.zs;
    nt = model.nt;
    ts = model.ts;
    soludoc = model.soludoc;
    
    % load pressures
    pw = zeros(nz, nx, ny);
    temp = load(fpwtxt);
    l = 1;
    for i = 1 : nz
        for k = 1 : ny
            for j = 1 : nx
                pw(i,j,k) = temp(l);
                l = l + 1;
            end
        end
    end
   
    % load x-direction velocities
    velx = zeros(nz, nx+1, ny);
    temp = load(fuxtxt);
    l = 1;
    for i = 1 : nz
        for k = 1 : ny  
            for j = 1 : nx+1
                velx(i,j,k) = temp(l);
                l = l + 1;
            end
        end
    end
    
    % load y-direction velocities
    vely = zeros(nz, nx, ny+1);
    temp = load(fuytxt);
    l = 1;
    for i = 1 : nz
        for k = 1 : ny+1
            for j = 1 : nx
                vely(i,j,k) = temp(l);
                l = l + 1;
            end
        end
    end
    
    % load z-direction velocities
    velz = zeros(nz+1, nx, ny);
    temp = load(fuztxt);
    l = 1;
    for i = 1 : nz+1
        for k = 1 : ny  
            for j = 1 : nx
                velz(i,j,k) = temp(l);
                l = l + 1;
            end
        end
    end
    
    % move the velocities on the faces into the center of the cell
    xvec = zeros(nz,nx,ny);
    yvec = zeros(nz,nx,ny);
    zvec = zeros(nz,nx,ny);
    for i = 1 : nz 
        for k = 1 : ny
            for j = 1 : nx
                xvec(i,j,k) = (velx(i,j,k)+velx(i, j+1, k)) / 2;
                yvec(i,j,k) = (vely(i,j,k)+vely(i, j, k+1)) / 2;
                zvec(i,j,k) = (velz(i,j,k)+velz(i+1, j, k)) / 2;
            end
        end
    end
    
    % compute the modulus of the velocity 
    vmodulus = zeros(nz,nx,ny);
    for i = 1 : nz 
        for k = 1 : ny   
            for j = 1 : nx
                vmodulus(i,j,k) = sqrt(xvec(i,j,k)^2+yvec(i,j,k)^2+zvec(i,j,k)^2);
            end
        end
    end
    
    % load mole fractions
    x = zeros(Nc, nz, nx, ny);    
    cdfmftxt = cellstr(fmftxt);    
    for m = 1 : Nc
        temp = load(char(cdfmftxt(m)));
        l = 1;
        for i = 1 : nz
            for k = 1 : ny
                for j = 1 : nx 
                    x(m, i, j, k) = temp(l);
                    l = l + 1;
                end
            end
        end
    end
    
    % load saturation of wetting phase
    sw = zeros(nz, nx, ny);
    temp = load(fswtxt);
    l = 1;
    for i = 1 : nz
        for k = 1 : ny   
            for j = 1 : nx
                sw(i,j,k) = temp(l);
                l = l + 1;
            end 
        end
    end
    
    % load molar density
    xi = zeros(nz, nx, ny);
    temp = load(fxitxt);
    l = 1;
    for i = 1 : nz
        for k = 1 : ny 
            for j = 1 : nx 
                xi(i,j,k) = temp(l);
                l = l + 1;
            end
        end
    end
    
    % load mole fraction history of component 1
    mfh = zeros(nt+1, 1);
    temp = load(fmhtxt);
    l = 1;
    for i = 2 : nt+1
        mfh(i) = temp(l);
        l = l + 1;
    end
    
    % load mole ratio history
    mfr = zeros(nt+1, 1);
    temp = load(fmrtxt);
    l = 1;
    for i = 2 : nt+1
        mfr(i) = temp(l);
        l = l + 1;
    end
    
    % construct a grid
    xcenter = zeros(nz,1);
    ycenter = zeros(nx,1);
    zcenter = zeros(ny,1);
    for i = 1 : nz
        xcenter(i) = (zs(i) + zs(i+1))/2;
    end
    for i = 1 : nx
        ycenter(i) = (xs(i) + xs(i+1))/2;
    end
    for i = 1 : ny
        zcenter(i) = (ys(i) + ys(i+1))/2;
    end
    
    % begin to draw the images
    fh = figure();
    h = title('Pressure field in the reservoir (unit: Pa)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('Z(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('X(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    h = zlabel('Y(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    axis equal
    hold on;
    
    [X,Y,Z] = meshgrid(xcenter, ycenter, zcenter); 
    xslice = [xcenter(1),xcenter(nz)]; 
    yslice = [ycenter(1),ycenter(nx)];  
    zslice = [zcenter(1),zcenter(ny)];
    h = slice(X, Y, Z, permute(pw,[2 1 3]), xslice, yslice, zslice);
    set(h,'EdgeColor','none','FaceColor','interp'); 
    t = colorbar;
    set(get(t,'title'), 'string', 'pres', 'Fontsize', 12);
    hold off;
    
    saveas(fh, [soludoc, '/Figure 1_pressure.fig']);

    fh = figure();
    h = title('Velocity component in x direction in the reservoir (unit: m/s)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('Z(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('X(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    h = zlabel('Y(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    axis equal
    hold on;
    
    [X,Y,Z] = meshgrid(xcenter, ycenter, zcenter); 
    xslice = [xcenter(1),xcenter(nz)]; 
    yslice = [ycenter(1),ycenter(nx)];  
    zslice = [zcenter(1),zcenter(ny)];
    h = slice(X, Y, Z, permute(xvec,[2 1 3]), xslice, yslice, zslice);
    set(h,'EdgeColor','none','FaceColor','interp'); 
    t = colorbar;
    set(get(t,'title'), 'string', 'velX', 'Fontsize', 12);
    hold off;
    
    saveas(fh, [soludoc, '/Figure 2_ux.fig']);
    
    fh = figure();
    h = title('Velocity component in y direction in the reservoir (unit: m/s)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('Z(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('X(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    h = zlabel('Y(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    axis equal
    hold on;
    
    [X,Y,Z] = meshgrid(xcenter, ycenter, zcenter); 
    xslice = [xcenter(1),xcenter(nz)]; 
    yslice = [ycenter(1),ycenter(nx)];  
    zslice = [zcenter(1),zcenter(ny)];
    h = slice(X, Y, Z, permute(yvec,[2 1 3]), xslice, yslice, zslice);
    set(h,'EdgeColor','none','FaceColor','interp'); 
    t = colorbar;
    set(get(t,'title'), 'string', 'velY', 'Fontsize', 12);
    hold off;
    
    saveas(fh, [soludoc, '/Figure 3_uy.fig']);
    
    fh = figure();
    h = title('Velocity component in z direction in the reservoir (unit: m/s)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('Z(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('X(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    h = zlabel('Y(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    axis equal
    hold on;
    
    [X,Y,Z] = meshgrid(xcenter, ycenter, zcenter); 
    xslice = [xcenter(1),xcenter(nz)]; 
    yslice = [ycenter(1),ycenter(nx)];  
    zslice = [zcenter(1),zcenter(ny)];
    h = slice(X, Y, Z, permute(zvec,[2 1 3]), xslice, yslice, zslice);
    set(h,'EdgeColor','none','FaceColor','interp'); 
    t = colorbar;
    set(get(t,'title'), 'string', 'velZ', 'Fontsize', 12);
    hold off;
    
    saveas(fh, [soludoc, '/Figure 4_uz.fig']);
 
    fh = figure();
    h = title('Modulus of velocity in the reservoir (unit: m/s)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('Z(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('X(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    h = zlabel('Y(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    axis equal
    hold on;
    
    [X,Y,Z] = meshgrid(xcenter, ycenter, zcenter); 
    xslice = [xcenter(1),xcenter(nz)]; 
    yslice = [ycenter(1),ycenter(nx)];  
    zslice = [zcenter(1),zcenter(ny)];
    h = slice(X, Y, Z, permute(vmodulus,[2 1 3]), xslice, yslice, zslice);
    set(h,'EdgeColor','none','FaceColor','interp'); 
    t = colorbar;
    set(get(t,'title'), 'string', '|v|', 'Fontsize', 12);
    hold off;
    
    saveas(fh, [soludoc, '/Figure 5_umodulus.fig']);

    fh = figure();
    h = title('Streamlines in the reservoir');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('Z(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('X(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    h = zlabel('Y(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    axis equal
    hold on;
    
    [X,Y,Z] = meshgrid(xcenter, ycenter, zcenter); 
    %xslice = [xcenter(1),xcenter(nz)]; 
    %yslice = [ycenter(1),ycenter(nx)];  
    %zslice = [zcenter(1),zcenter(ny)];
    xslice = [xcenter(nz)]; 
    yslice = [ycenter(1)];  
    zslice = [zcenter(ny)];
    h = streamslice(X, Y, Z, permute(xvec,[2 1 3]), permute(yvec,[2 1 3]), permute(zvec,[2 1 3]), xslice, yslice, zslice, 2);
    set(h, 'color', 'black');
    hold off;
    
    saveas(fh, [soludoc, '/Figure 6_streamline.fig']);
    
    for m = 1 : Nc
        fh = figure();
        strtitle = ['Fraction mole field of component ', num2str(m), ' in the reservoir'];
        h = title(strtitle);
        set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        h = xlabel('Z(m)');
        set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        h = ylabel('X(m)');
        set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
        h = zlabel('Y(m)');
        set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        axis equal
        hold on;
    
        [X,Y,Z] = meshgrid(xcenter, ycenter, zcenter); 
        xslice = [xcenter(1),xcenter(nz)]; 
        yslice = [ycenter(1),ycenter(nx)];  
        zslice = [zcenter(1),zcenter(ny)];
        x3d = zeros(nz,nx,ny);
        for i = 1 : nz       
            for j = 1 : nx
                for k = 1 : ny
                    x3d(i,j,k) = x(m,i,j,k);
                end
            end
        end
        h = slice(X, Y, Z, permute(x3d,[2 1 3]), xslice, yslice, zslice);
        set(h,'EdgeColor','none','FaceColor','interp'); 
        t = colorbar;
        set(get(t,'title'), 'string', 'fraction', 'Fontsize', 12);
        hold off;
        
        fstr = [soludoc, '/Figure 7-', num2str(m), '_molefraction.fig'];
        saveas(fh,fstr);
    end
    
    fh = figure();
    h = title('Saturation field in the reservoir');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('Z(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('X(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    h = zlabel('Y(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    axis equal
    hold on;
    
    [X,Y,Z] = meshgrid(xcenter, ycenter, zcenter); 
    xslice = [xcenter(1),xcenter(nz)]; 
    yslice = [ycenter(1),ycenter(nx)];  
    zslice = [zcenter(1),zcenter(ny)];
    h = slice(X, Y, Z, permute(sw,[2 1 3]), xslice, yslice, zslice);
    set(h,'EdgeColor','none','FaceColor','interp'); 
    t = colorbar;
    set(get(t,'title'), 'string', 'satur', 'Fontsize', 12);
    hold off;
    
    saveas(fh, [soludoc, '/Figure 8_saturation.fig']);
    
    fh = figure();
    h = title('Molar density field in the reservoir (unit: mol/m^3)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('Z(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('X(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    h = zlabel('Y(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    axis equal
    hold on;
    
    [X,Y,Z] = meshgrid(xcenter, ycenter, zcenter); 
    xslice = [xcenter(1),xcenter(nz)]; 
    yslice = [ycenter(1),ycenter(nx)];  
    zslice = [zcenter(1),zcenter(ny)];
    h = slice(X, Y, Z, permute(xi,[2 1 3]), xslice, yslice, zslice);
    set(h,'EdgeColor','none','FaceColor','interp'); 
    t = colorbar;
    set(get(t,'title'), 'string', 'xi', 'Fontsize', 12);
    hold off;
    
    saveas(fh, [soludoc, '/Figure 9_molar density.fig']);
    
    
    fh = figure();
    h = title('Issue mole fraction history of Component 2 in the reservoir');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('Time(s)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Mole fraction(%)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    hold on;
    
    plot(ts, mfh*100);
    hold off;
    
    saveas(fh, [soludoc, '/Figure 10_issuemolehistory.fig']);
    
    fh = figure();
    h = title('Mole ratio history in the reservoir');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('Time(s)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Mole ratio');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    hold on;
    
    plot(ts, mfr);
    hold off;
    
    saveas(fh, [soludoc, '/Figure 11_moleratio.fig']);
    
end

