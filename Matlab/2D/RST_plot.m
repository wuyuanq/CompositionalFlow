% This is the RST_plot() function which uses the model in the input file 
% and the solution files output by the function RST_compositionalTwoPhaseFlow()
% as input, and generates a series of images to demonstrate the results.

% Input parameters:
% model: the model
% fpwtxt: the pressure file in the txt form
% fuxtxt: the x-dirction velocity file in the txt form
% fuytxt: the y-dirction velocity file in the txt form
% fmftxt: the mole fraction files in the txt form
% fswtxt: the saturation of wetting phase file in the txt form
% fxitxt: the molar density in the txt form
% fmhtxt: the mole fraction history file of component 1 in the txt form
% fmrtxt: the mole ratio history file of desired components to component 1
% ftimetxt: the time when finishing each iteration

% Author: Yuanqing Wu. Email: wuyuanq@gmail.com
% Last edited on January 23rd, 2014

function RST_plot( model, fpwtxt, fuxtxt, fuytxt, fmftxt, fswtxt, fxitxt, fmhtxt, fmrtxt)%, ftimetxt )

    % simplify the symbols of the model
    Nc = model.Nc;
    nx = model.nx;
    ny = model.ny;
    nt = model.nt;
    xs = model.xs;
    ys = model.ys;
    ts = model.ts;
    soludoc = model.soludoc;
    
    % load pressures
    pw = zeros(ny, nx);
    temp = load(fpwtxt);
    k = 1;
    for j = 1 : nx
        for i = 1 : ny 
            pw(i,j) = temp(k);
            k = k + 1;
        end
    end
    
    % load x-direction velocities
    velx = zeros(ny, nx+1);
    temp = load(fuxtxt);
    k = 1;
    for j = 1 : nx+1
        for i = 1 : ny  
            velx(i,j) = temp(k);
            k = k + 1;
        end
    end
    
    % load y-direction velocities
    vely = zeros(ny+1, nx);
    temp = load(fuytxt);
    k = 1;
    for j = 1 : nx
        for i = 1 : ny+1
            vely(i,j) = temp(k);
            k = k + 1;
        end
    end
    
    % move the velocities on the edges into the center of the cell
    xvec = zeros(ny, nx);
    yvec = zeros(ny, nx);
    for j = 1 : nx
        for i = 1 : ny  
            xvec(i,j) = (velx(i,j)+velx(i,j+1))/2;
            yvec(i,j) = (vely(i,j)+vely(i+1,j))/2;
        end
    end
    
    % compute the modulus of the velocity 
    vmodulus = zeros(ny, nx);
    for j = 1 : nx
        for i = 1 : ny   
            vmodulus(i,j) = sqrt(xvec(i,j)^2+yvec(i,j)^2);
        end
    end
    
    % load mole fractions
    x = zeros(Nc, ny, nx);    
    cdfmftxt = cellstr(fmftxt);    
    for m = 1 : Nc
        temp = load(char(cdfmftxt(m)));
        k = 1;
        for j = 1 : nx
            for i = 1 : ny 
                x(m, i, j) = temp(k);
                k = k + 1;
            end
        end
    end
    
    % load saturation of wetting phase
    sw = zeros(ny, nx);
    temp = load(fswtxt);
    k = 1;
    for j = 1 : nx
        for i = 1 : ny
            sw(i,j) = temp(k);
            k = k + 1;
        end
    end
    
    % load molar density
    xi = zeros(ny, nx);
    temp = load(fxitxt);
    k = 1;
    for j = 1 : nx
        for i = 1 : ny
            xi(i,j) = temp(k);
            k = k + 1;
        end
    end
    
    % load mole fraction history of component 2
    mfh = zeros(nt+1, 1);
    temp = load(fmhtxt);
    k = 1;
    for i = 2 : nt+1
        mfh(i) = temp(k);
        k = k + 1;
    end
    mfh(1) = 0;
    
    % load mole ratio history
    mfr = zeros(nt+1, 1);
    temp = load(fmrtxt);
    k = 1;
    for i = 2 : nt+1
        mfr(i) = temp(k);
        k = k + 1;
    end
    mfr(1) = mfr(2);
    
    % load iteration time
    %time = zeros(nt+1, 1);
    %temp = load(ftimetxt);
    %k = 1;
    %for i = 2 : nt+1
    %    time(i) = temp(k);
    %    k = k + 1;
    %end
    %time(1) = 0;
    
    % construct a grid
    xcenter = zeros(nx, 1);
    ycenter = zeros(ny, 1);
    for i = 1 : nx
        xcenter(i) = (xs(i) + xs(i+1))/2;
    end
    for i = 1 : ny
        ycenter(i) = (ys(i) + ys(i+1))/2;
    end
    
    % begin to draw the images
    fh = figure();
    h = title('Pressure (Pa)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('X(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Y(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    axis equal
    hold on;
    
    [X,Y] = meshgrid(xcenter, ycenter); 
    contourf(X, Y, pw, 200, 'linecolor', 'none');
    t = colorbar;
    set(get(t,'title'), 'string', 'pres', 'Fontsize', 12);
    hold off;
    
    saveas(fh, [soludoc, '/Figure 1_pressure.fig']);

    fh = figure();
    h = title('Darcy velocity (m/s) (x-component)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('X(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Y(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    axis equal
    hold on;
    
    [X,Y] = meshgrid(xs, ycenter); 
    contourf(X, Y, velx, 200, 'linecolor', 'none');
    t = colorbar;
    set(get(t,'title'), 'string', 'velX', 'Fontsize', 12);
    hold off;
    
    saveas(fh, [soludoc, '/Figure 2_ux.fig']);
    
    fh = figure();
    h = title('Darcy velocity (m/s) (y-component)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('X(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Y(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    axis equal
    hold on;
    
    [X,Y] = meshgrid(xcenter, ys); 
    contourf(X, Y, vely, 200, 'linecolor', 'none');
    t = colorbar;
    set(get(t,'title'), 'string', 'velY', 'Fontsize', 12);
    hold off;
    
    saveas(fh, [soludoc, '/Figure 3_uy.fig']);
 
    fh = figure();
    h = title('Modulus of Darcy velocity (m/s)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('X(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Y(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    axis equal
    hold on;
    
    [X,Y] = meshgrid(xcenter, ycenter); 
    contourf(X, Y, vmodulus, 200, 'linecolor', 'none');
    t = colorbar;
    set(get(t,'title'), 'string', '|v|', 'Fontsize', 12);
    hold off;
    
    saveas(fh, [soludoc, '/Figure 4_umodulus.fig']);

    fh = figure();
    h = title('Streamlines of Darcy velocity');
    %h = title('Streamlines of Darcy velocity at the end of 3 year');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('X(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Y(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    axis equal
    hold on;
    
    [X,Y] = meshgrid(xcenter, ycenter);
    h = streamslice(X, Y, xvec, yvec, 2);
    set(h, 'color', 'black');
    hold off;
    
    saveas(fh, [soludoc, '/Figure 5_streamline.fig']);
    
    for m = 1 : Nc
        fh = figure();
        strtitle = ['Mole fraction of component ', num2str(m)];
        %strtitle = ['Mole fraction of propane at the end of 3 year'];
        h = title(strtitle);
        set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        h = xlabel('X(m)');
        set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        h = ylabel('Y(m)');
        set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
        axis equal
        hold on;
    
        [X,Y] = meshgrid(xcenter, ycenter); 
        x2d = zeros(ny, nx); % because contourf can only draw 2-D matrix, we change x from 3-D to 2-D.
        x2d(1:end, 1:end) = x(m, 1:end, 1:end);
        contourf(X, Y, x2d, 200, 'linecolor', 'none');
        t = colorbar;
        set(get(t,'title'), 'string', 'fraction', 'Fontsize', 12);
        hold off;
        
        fstr = [soludoc, '/Figure 6-', num2str(m), '_molefraction.fig'];
        saveas(fh,fstr);
    end
    
    fh = figure();
    h = title('Saturation of the oil phase');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('X(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Y(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    axis equal
    hold on;
    
    [X,Y] = meshgrid(xcenter, ycenter); 
    contourf(X, Y, sw, 200, 'linecolor', 'none');
    t = colorbar;
    set(get(t,'title'), 'string', 'saturation', 'Fontsize', 12);
    hold off;
    
    saveas(fh, [soludoc, '/Figure 7_saturation.fig']);
    
    fh = figure();
    h = title('Molar concentration of mixture (mol/m^2)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('X(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Y(m)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    axis equal
    hold on;
    
    [X,Y] = meshgrid(xcenter, ycenter); 
    contourf(X, Y, xi, 200, 'linecolor', 'none');
    t = colorbar;
    set(get(t,'title'), 'string', 'molar density', 'Fontsize', 12);
    hold off;
    
    saveas(fh, [soludoc, '/Figure 8_molarconcentration.fig']);
    
    fh = figure();
    h = title('Recovery rate of Component 3');
    %h = title('Issue mole fraction history of propane in 3 year');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('Time(s)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Recovery rate(%)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    hold on;
    
    plot(ts, mfh*100);
    hold off;
    
    saveas(fh, [soludoc, '/Figure 9_issuemolehistory.fig']);
    
    fh = figure();
    h = title('Mole ratio history');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = xlabel('Time(s)');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    h = ylabel('Mole ratio');
    set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    hold on;
    
    plot(ts, mfr);
    hold off;
    
    saveas(fh, [soludoc, '/Figure 10_moleratio.fig']);
    
    %fh = figure();
    %h = title('Iteration time from a CCFD simulation of compositional flow');
    %set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    %h = xlabel('Iteration');
    %set(h, 'fontsize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    %h = ylabel('Time(s)');
    %set(h, 'fontsize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    %hold on;
    
    %plot(ts, time);
    %hold off;
    
    %saveas(fh, [soludoc, '/Figure 11_time.fig']);
    
end

