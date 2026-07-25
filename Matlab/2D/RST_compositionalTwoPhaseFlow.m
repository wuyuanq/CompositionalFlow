% This is the RST_twoPhaseFlow() function which can be called by the
% input files. It computes the pressures, saturations of wetting phase 
% and then gives out the velocities. After that, it outputs the results 
% to a series of solution files with the help of function
% RST_writeFile() and then calls the RST_plot() function to draw the 
% images of the results;

% Input parameters:
% model: the model

% Author: Yuanqing Wu. Email: wuyuanq@gmail.com
% Last edited on January 17th, 2014

function RST_compositionalTwoPhaseFlow( model )

    % simplify the symbols of the model
    Nc = model.Nc;
    nx = model.nx;
    ny = model.ny;
    nt = model.nt;
    xs = model.xs;
    ys = model.ys;
    ts = model.ts;
    Kxx = model.Kxx;
    Kyy = model.Kyy;
    srcm = model.src; 
    PwBdryX = model.PwBdryX;
    PwBdryY = model.PwBdryY;
    PwInit = model.PwInit;
    UwBdryX = model.UwBdryX;
    UwBdryY = model.UwBdryY;
    UnBdryX = model.UnBdryX;
    UnBdryY = model.UnBdryY;
    zInit = model.zInit;
    soludoc = model.soludoc;
    
    % create the new document to store the results
    if(exist(soludoc, 'dir') == 0)
        mkdir(soludoc);
    end
    
    % the file to store the issue mole fraction history
    fmhtxt = [soludoc, '/soln_cmp2PhFlw_moleHistory.txt'];
    fmhtxtid = fopen(fmhtxt, 'w');
    
    % the file to store the left mole ratio history
    fmrtxt = [soludoc, '/soln_cmp2PhFlw_moleRatioHistory.txt'];
    fmrtxtid = fopen(fmrtxt, 'w');
  
    % rectify the directions of the velocities
    UwBdryX(1, 1:end) = -UwBdryX(1, 1:end);
    UwBdryY(1:end, 1) = -UwBdryY(1:end, 1);
    UnBdryX(1, 1:end) = -UnBdryX(1, 1:end);
    UnBdryY(1:end, 1) = -UnBdryY(1:end, 1);
    
    % define pressures of wetting phase and initialize them
    Pw = zeros(ny+2, nx+2);
    Pw(1, 2:nx+1) = PwBdryY(1:end, 1);
    Pw(ny+2, 2:nx+1) = PwBdryY(1:end, 2);
    Pw(2:ny+1, 1) = PwBdryX(1, 1:end);
    Pw(2:ny+1, nx+2) = PwBdryX(2, 1:end);
    Pw(2:end-1,2:end-1) = PwInit(1:end, 1:end)';
    
    % define the velocities of two phases and initialize them
    Uwx = zeros(ny, nx+1);
    Uwy = zeros(ny+1, nx);
    Uwx(1:end, 1) = UwBdryX(1, 1:end);
    Uwx(1:end, nx+1) = UwBdryX(2, 1:end);
    Uwy(1, 1:end) = UwBdryY(1:end, 1);
    Uwy(ny+1, 1:end) = UwBdryY(1:end, 2);
    Unx = zeros(ny, nx+1);
    Uny = zeros(ny+1, nx);
    Unx(1:end, 1) = UnBdryX(1, 1:end);
    Unx(1:end, nx+1) = UnBdryX(2, 1:end);
    Uny(1, 1:end) = UnBdryY(1:end, 1);
    Uny(ny+1, 1:end) = UnBdryY(1:end, 2);
        
    % define the saturation of wetting phase in the cells
    Sw = zeros(ny, nx);
    
    % define the mobilities on the bars
    lambdawx = zeros(ny, nx+1);
    lambdawy = zeros(ny+1, nx);
    lambdanx = zeros(ny, nx+1);
    lambdany = zeros(ny+1, nx);
    
    % define the permeabilities on the bars
    Kxxbar = zeros(ny, nx+1);
    Kyybar = zeros(ny+1, nx);
    
    % initialize Kxxbar using harmonic weighting method
    Kxxbar(1:ny,1) = Kxx(1,1:ny);
    Kxxbar(1:ny,nx+1) = Kxx(nx, 1:ny);
    for i = 1 : ny
        for j = 2 : nx  
            ltotal = xs(j+1) - xs(j-1);
            lleft = xs(j) - xs(j-1);
            lright = xs(j+1) - xs(j);
            Kxxbar(i,j) = ltotal / (lleft/Kxx(j-1,i)+lright/Kxx(j,i));
        end
    end
    
    % initialize Kyybar using harmonic weighting method
    Kyybar(1,1:nx) = Kyy(1:nx,1);
    Kyybar(ny+1,1:nx) = Kyy(1:nx, ny);
    for i = 2 : ny
        for j = 1 : nx        
            ltotal = ys(i+1) - ys(i-1);
            lup = ys(i+1) - ys(i);
            ldown = ys(i) - ys(i-1);
            Kyybar(i,j) = ltotal / (ldown/Kyy(j, i-1)+lup/Kyy(j,i));
        end
    end
    
    % define the mole fractions in the cells and initialize them
    z = zeros(Nc, ny, nx);
    for m = 1 : Nc
        for i = 1 : ny
            for j = 1 : nx
                z(m,i,j) = zInit(m,j,i);
            end
        end
    end
    
    % define the mass densities of wetting phase and nonwetting phase and  
    % initialize them
    densiW = zeros(ny, nx);
    densiN = zeros(ny, nx);
    
    % define the mass density on the bars
    densiWbarx = zeros(ny, nx+1);
    densiWbary = zeros(ny+1, nx);
    densiNbarx = zeros(ny, nx+1);
    densiNbary = zeros(ny+1, nx);
    
    % define the sources of wetting and nonwetting phases and compute their
    % values
    src = zeros(Nc, ny, nx);
    for m = 1 : Nc
        for i = 1 : ny
            for j = 1 : nx
                src(m,i,j) = srcm(m,j,i);
            end
        end
    end
    
    % define the molar fractions of wetting and nonwetting phase in the
    % cells
    xW = zeros(Nc, ny, nx);
    xN = zeros(Nc, ny, nx);
    
    % define the molar fractions of wetting and nonwetting phase on the
    % bars
    xWbarx = zeros(Nc, ny, nx+1);
    xWbary = zeros(Nc, ny+1, nx);
    xNbarx = zeros(Nc, ny, nx+1);
    xNbary = zeros(Nc, ny+1, nx);
    
    % define the molar density of each phase in the cells
    xiW = zeros(ny, nx);
    xiN = zeros(ny, nx);
    xi = zeros(ny, nx); % the fluid molar density
    
    % define the molar density of each phase on the bars
    xiWbarx = zeros(ny, nx+1);
    xiWbary = zeros(ny+1, nx);
    xiNbarx = zeros(ny, nx+1);
    xiNbary = zeros(ny+1, nx);
    
    % define the partial molar volume in the cells
    v = zeros(Nc, ny, nx);
    
    % the total mixture compressibility
    Cf = zeros(ny, nx);
    
    % define the viscosity of wetting and nonwetting phase
    viscW = zeros(ny, nx);
    viscN = zeros(ny, nx);
    
    % define the left mole of each component
    leftmole = zeros(Nc, 1);
    
    firsttime = 1;
    totalmole = 0.0;
    timestep = ts(2) - ts(1);
    
    % time iteration 
    for t = 2 : nt+1
        t
        [ xW, xN, xiW, xiN, densiW, densiN, Sw, v, Cf, viscW, viscN, xi, lambdawx, lambdanx, xiWbarx, xiNbarx, densiWbarx, ...
            densiNbarx, xWbarx, xNbarx, lambdawy, lambdany, xiWbary, xiNbary, densiWbary, densiNbary, xWbary, xNbary, ...
            firsttime, totalmole ] = computeParameters( model, Pw, z, Uwx, Uwy, Unx, Uny, xW, xN, xiW, xiN, xi, ...
            densiW, densiN, Sw, v, Cf, viscW, viscN, lambdawx, lambdanx, xiWbarx, xiNbarx, densiWbarx, densiNbarx, xWbarx, ...
            xNbarx, lambdawy, lambdany, xiWbary, xiNbary, densiWbary, densiNbary, xWbary, xNbary, Kxxbar, Kyybar, ...
            firsttime, totalmole, leftmole, fmhtxtid, fmrtxtid );

        [ Pw ] = computePres( model, Pw, lambdawx, lambdanx, lambdawy, lambdany, ...
            xiWbarx, xiNbarx, xiWbary, xiNbary, densiWbarx, densiWbary, densiNbarx, densiNbary, xWbarx, xNbarx, xWbary, ...
            xNbary, v, Cf, src, timestep );
        
        [ Uwx, Uwy, Unx, Uny ] = computeVel( model, Pw, lambdawx, lambdanx, lambdawy, lambdany, ...
            densiWbarx, densiWbary, densiNbarx, densiNbary );
      
        [ z ] = computez( model, timestep, xi, z, src, Uwx, Uwy, Unx, Uny, xiWbarx, xiWbary, xiNbarx, ...
            xiNbary, xWbarx, xWbary, xNbarx, xNbary );  
        
    end
    
    % close the mole history file
    fclose(fmhtxtid);
    
    % close the mole rate history file
    fclose(fmrtxtid);
    
    % write the results to the solution files
    [ fpwtxt, fuxtxt, fuytxt, fmftxt, fswtxt, fxitxt ] = RST_writeFile( model, Pw, Uwx+Unx, Uwy+Uny, z, Sw, xi );
    
    % draw the result images
    RST_plot( model, fpwtxt, fuxtxt, fuytxt, fmftxt, fswtxt, fxitxt, fmhtxt, fmrtxt ); 
    
end

function [ xW, xN, xiW, xiN, densiW, densiN, Sw, v, Cf, viscW, viscN, xi, lambdawx, lambdanx, xiWbarx, xiNbarx, densiWbarx, ...
    densiNbarx, xWbarx, xNbarx, lambdawy, lambdany, xiWbary, xiNbary, densiWbary, densiNbary, xWbary, xNbary, ...
    firsttime, totalmole] = computeParameters( model, Pw, z, Uwx, Uwy, Unx, Uny, xW, xN, ...
    xiW, xiN, xi, densiW, densiN, Sw, v, Cf, viscW, viscN, lambdawx, lambdanx, xiWbarx, xiNbarx, densiWbarx, densiNbarx, ...
    xWbarx, xNbarx, lambdawy, lambdany, xiWbary, xiNbary, densiWbary, densiNbary, xWbary, xNbary, Kxxbar, ...
    Kyybar, firsttime, totalmole, leftmole, fmhtxtid, fmrtxtid )

    Nc = model.Nc;
    nx = model.nx;
    ny = model.ny;
    kr_W = model.kr_W;
    kr_N = model.kr_N; 
    zBdryX = model.zBdryX;
    zBdryY = model.zBdryY;
    isDiriX = model.isDiriX;
    isDiriY = model.isDiriY;
    
    % compute a series of parameters of wetting and nonwetting phase
    for i = 1 : ny
        for j = 1 : nx

            [ xW(1:Nc,i,j), xN(1:Nc,i,j), xiW(i,j), xiN(i,j), densiW(i,j), densiN(i,j), Sw(i,j), v(1:Nc,i,j), Cf(i,j), isW, isN ] ...
                = RST_flashcalculation( model, Pw(i+1,j+1), z(1:Nc,i,j) );    

            if(isW)
                [ viscW(i,j) ] = RST_viscosity( model, xW(1:Nc,i,j), xiW(i,j), Pw(i+1,j+1), 'l' );
            else
                viscW(i,j) = inf;
            end    
            if(isN)
                [ viscN(i,j) ] = RST_viscosity( model, xN(1:Nc,i,j), xiN(i,j), Pw(i+1,j+1), 'g' );
            else
                viscN(i,j) = inf;
            end
            
        end
    end   
    
    % compute the fluid mole density
    xi(1:end,1:end) = xiW(1:end,1:end).*Sw(1:end,1:end) + xiN(1:end,1:end).*(1-Sw(1:end,1:end));
    
    % compute the left mole of the components
    leftmole(1:Nc) = 0;
    for m = 1 : Nc
        for i = 1 : ny
            for j = 1 : nx
                leftmole(m) = leftmole(m) + xi(i,j)*z(m,i,j);
            end
        end 
    end
    
    if(firsttime == 1)
        firsttime = 0;
        for m = 2 : Nc
            totalmole = totalmole + leftmole(m);
        end
    end
        
    totaldesiredleftmole = 0.0;
    for m = 2 : Nc
        totaldesiredleftmole = totaldesiredleftmole + leftmole(m);
    end
        
    fprintf(fmhtxtid, '%10e\n', (totalmole-totaldesiredleftmole)/totalmole);
             
    % print the mole ratio of desired components to component 1 in the well
    fprintf(fmrtxtid, '%10e\n', totaldesiredleftmole/leftmole(1));
        
    % compute the mobilities, the fractional flow functions, the molar densities,
    % the molar fractions on the bars, using the single-point upstream weighting.
    % lambda means: mobility = permeability/viscosity
    for i = 1 : ny
        if((Uwx(i,1) > 0) || (Unx(i,1) > 0))
            [ xWbarx(1:Nc,i,1), xNbarx(1:Nc,i,1), xiWbarx(i,1), xiNbarx(i,1), densiWbarx(i,1), densiNbarx(i,1), Swtemp, ~, ~, isW, isN ] ...
                = RST_flashcalculation( model, Pw(i+1,2), zBdryX(1:Nc,1,i) );   
            if(isW)
                [ viscWtemp ] = RST_viscosity( model, xWbarx(1:Nc,i,1), xiWbarx(i,1), Pw(i+1,2), 'l' );
            end
            if(isN)
                [ viscNtemp ] = RST_viscosity( model, xNbarx(1:Nc,i,1), xiNbarx(i,1), Pw(i+1,2), 'g' );
            end
        end
        if(Uwx(i,1) > 0)    
            lambdawx(i,1) = Kxxbar(i,1)*kr_W(Swtemp)/viscWtemp; 
        elseif(Uwx(i,1) < 0)  
            lambdawx(i,1) = Kxxbar(i,1)*kr_W(Sw(i,1))/viscW(i,1);                
            xiWbarx(i,1) = xiW(i,1);                    
            densiWbarx(i,1) = densiW(i,1);              
            xWbarx(1:Nc,i,1) = xW(1:Nc,i,1);  
        else
            if(isDiriX(1,i) == 1)
                lambdawx(i,1) = Kxxbar(i,1)*kr_W(Sw(i,1))/viscW(i,1);
                xiWbarx(i,1) = xiW(i,1);
                densiWbarx(i,1) = densiW(i,1);
                xWbarx(1:Nc,i,1) = xW(1:Nc,i,1);
            end 
        end
        if(Unx(i,1) > 0)          
            lambdanx(i,1) = Kxxbar(i,1)*kr_N(Swtemp)/viscNtemp; 
        elseif(Unx(i,1) < 0) 
            lambdanx(i,1) = Kxxbar(i,1)*kr_N(Sw(i,1))/viscN(i,1);
            xiNbarx(i,1) = xiN(i,1);
            densiNbarx(i,1) = densiN(i,1);
            xNbarx(1:Nc,i,1) = xN(1:Nc,i,1);
        else
            if(isDiriX(1,i) == 1) 
                lambdanx(i,1) = Kxxbar(i,1)*kr_N(Sw(i,1))/viscN(i,1);
                xiNbarx(i,1) = xiN(i,1);
                densiNbarx(i,1) = densiN(i,1);
                xNbarx(1:Nc,i,1) = xN(1:Nc,i,1);
            end 
        end          
    end
    
    for i = 1 : ny
        if((Uwx(i,nx+1) < 0) || (Unx(i,nx+1) < 0))
            [ xWbarx(1:Nc,i,nx+1), xNbarx(1:Nc,i,nx+1), xiWbarx(i,nx+1), xiNbarx(i,nx+1), densiWbarx(i,nx+1), densiNbarx(i,nx+1), Swtemp, ~, ~, isW, isN ] ...
                = RST_flashcalculation( model, Pw(i+1,nx+1), zBdryX(1:Nc,2,i) );   
            if(isW)
                [ viscWtemp ] = RST_viscosity( model, xWbarx(1:Nc,i,nx+1), xiWbarx(i,nx+1), Pw(i+1,nx+1), 'l' );
            end
            if(isN)
                [ viscNtemp ] = RST_viscosity( model, xNbarx(1:Nc,i,nx+1), xiNbarx(i,nx+1), Pw(i+1,nx+1), 'g' );
            end
        end
        if(Uwx(i,nx+1) < 0)    
            lambdawx(i,nx+1) = Kxxbar(i,nx+1)*kr_W(Swtemp)/viscWtemp; 
        elseif(Uwx(i,nx+1) > 0)
            lambdawx(i,nx+1) = Kxxbar(i,nx+1)*kr_W(Sw(i,nx))/viscW(i,nx);                
            xiWbarx(i,nx+1) = xiW(i,nx);                    
            densiWbarx(i,nx+1) = densiW(i,nx);              
            xWbarx(1:Nc,i,nx+1) = xW(1:Nc,i,nx);  
        else
            if(isDiriX(2,i) == 1) 
                lambdawx(i,nx+1) = Kxxbar(i,nx+1)*kr_W(Sw(i,nx))/viscW(i,nx);
                xiWbarx(i,nx+1) = xiW(i,nx);
                densiWbarx(i,nx+1) = densiW(i,nx);
                xWbarx(1:Nc,i,nx+1) = xW(1:Nc,i,nx);
            end 
        end
        if(Unx(i,nx+1) < 0)          
            lambdanx(i,nx+1) = Kxxbar(i,nx+1)*kr_N(Swtemp)/viscNtemp; 
        elseif(Unx(i,nx+1) > 0)
            lambdanx(i,nx+1) = Kxxbar(i,nx+1)*kr_N(Sw(i,nx))/viscN(i,nx);
            xiNbarx(i,nx+1) = xiN(i,nx);
            densiNbarx(i,nx+1) = densiN(i,nx);
            xNbarx(1:Nc,i,nx+1) = xN(1:Nc,i,nx);
        else
            if(isDiriX(2,i) == 1)
                lambdanx(i,nx+1) = Kxxbar(i,nx+1)*kr_N(Sw(i,nx))/viscN(i,nx);
                xiNbarx(i,nx+1) = xiN(i,nx);
                densiNbarx(i,nx+1) = densiN(i,nx);
                xNbarx(1:Nc,i,nx+1) = xN(1:Nc,i,nx);
            end
        end          
    end
                
    for i = 1 : ny
        for j = 2 : nx              
            if(Uwx(i,j) > 0)
                lambdawx(i,j) = Kxxbar(i,j)*kr_W(Sw(i,j-1))/viscW(i,j-1);
                xiWbarx(i,j) = xiW(i,j-1);
                densiWbarx(i,j) = densiW(i,j-1);
                xWbarx(1:Nc,i,j) = xW(1:Nc,i,j-1);              
            else
                lambdawx(i,j) = Kxxbar(i,j)*kr_W(Sw(i,j))/viscW(i,j);
                xiWbarx(i,j) = xiW(i,j);
                densiWbarx(i,j) = densiW(i,j);
                xWbarx(1:Nc,i,j) = xW(1:Nc,i,j);     
            end
            if(Unx(i,j) > 0)
                lambdanx(i,j) = Kxxbar(i,j)*kr_N(Sw(i,j-1))/viscN(i,j-1);
                xiNbarx(i,j) = xiN(i,j-1);
                densiNbarx(i,j) = densiN(i,j-1);
                xNbarx(1:Nc,i,j) = xN(1:Nc,i,j-1);              
            else
                lambdanx(i,j) = Kxxbar(i,j)*kr_N(Sw(i,j))/viscN(i,j);
                xiNbarx(i,j) = xiN(i,j);
                densiNbarx(i,j) = densiN(i,j);
                xNbarx(1:Nc,i,j) = xN(1:Nc,i,j);         
            end
        end
    end
    
    for i = 1 : nx
        if((Uwy(1,i) > 0) || (Uny(1,i) > 0))
            [ xWbary(1:Nc,1,i), xNbary(1:Nc,1,i), xiWbary(1,i), xiNbary(1,i), densiWbary(1,i), densiNbary(1,i), Swtemp, ~, ~, isW, isN ] ...
                = RST_flashcalculation( model, Pw(2,i+1), zBdryY(1:Nc,i,1) );  
            if(isW)
                [ viscWtemp ] = RST_viscosity( model, xWbary(1:Nc,1,i), xiWbary(1,i), Pw(2,i+1), 'l' );
            end
            if(isN)
                [ viscNtemp ] = RST_viscosity( model, xNbary(1:Nc,1,i), xiNbary(1,i), Pw(2,i+1), 'g' );
            end
        end
        if(Uwy(1,i) > 0)    
            lambdawy(1,i) = Kyybar(1,i)*kr_W(Swtemp)/viscWtemp; 
        elseif(Uwy(1,i) < 0)
            lambdawy(1,i) = Kyybar(1,i)*kr_W(Sw(1,i))/viscW(1,i);                
            xiWbary(1,i) = xiW(1,i);                    
            densiWbary(1,i) = densiW(1,i);              
            xWbary(1:Nc,1,i) = xW(1:Nc,1,i);  
        else
            if(isDiriY(i,1) == 1)
                lambdawy(1,i) = Kyybar(1,i)*kr_W(Sw(1,i))/viscW(1,i);
                xiWbary(1,i) = xiW(1,i);
                densiWbary(1,i) = densiW(1,i);
                xWbary(1:Nc,1,i) = xW(1:Nc,1,i);
            end 
        end
        if(Uny(1,i) > 0)          
            lambdany(1,i) = Kyybar(1,i)*kr_N(Swtemp)/viscNtemp; 
        elseif(Uny(1,i) < 0)
            lambdany(1,i) = Kyybar(1,i)*kr_N(Sw(1,i))/viscN(1,i);
            xiNbary(1,i) = xiN(1,i);
            densiNbary(1,i) = densiN(1,i);
            xNbary(1:Nc,1,i) = xN(1:Nc,1,i);
        else
            if(isDiriY(i,1) == 1) 
                lambdany(1,i) = Kyybar(1,i)*kr_N(Sw(1,i))/viscN(1,i);
                xiNbary(1,i) = xiN(1,i);
                densiNbary(1,i) = densiN(1,i);
                xNbary(1:Nc,1,i) = xN(1:Nc,1,i);
            end 
        end          
    end
    
    for i = 1 : nx
        if((Uwy(ny+1,i) < 0) || (Uny(ny+1,i) < 0))
            [ xWbary(1:Nc,ny+1,i), xNbary(1:Nc,ny+1,i), xiWbary(ny+1,i), xiNbary(ny+1,i), densiWbary(ny+1,i), densiNbary(ny+1,i), Swtemp, ~, ~, isW, isN ] ...
                = RST_flashcalculation( model, Pw(ny+1,i+1), zBdryY(1:Nc,i,2) ); 
            if(isW)
                [ viscWtemp ] = RST_viscosity( model, xWbary(1:Nc,ny+1,i), xiWbary(ny+1,i), Pw(ny+1,i+1), 'l' );
            end
            if(isN)
                [ viscNtemp ] = RST_viscosity( model, xNbary(1:Nc,ny+1,i), xiNbary(ny+1,i), Pw(ny+1,i+1), 'g' );
            end
        end
        if(Uwy(ny+1,i) < 0)    
            lambdawy(ny+1,i) = Kyybar(ny+1,i)*kr_W(Swtemp)/viscWtemp; 
        elseif(Uwy(ny+1,i) > 0)
            lambdawy(ny+1,i) = Kyybar(ny+1,i)*kr_W(Sw(ny,i))/viscW(ny,i);                
            xiWbary(ny+1,i) = xiW(ny,i);                    
            densiWbary(ny+1,i) = densiW(ny,i);              
            xWbary(1:Nc,ny+1,i) = xW(1:Nc,ny,i);  
        else
            if(isDiriY(i,2) == 1)
                lambdawy(ny+1,i) = Kyybar(ny+1,i)*kr_W(Sw(ny,i))/viscW(ny,i);
                xiWbary(ny+1,i) = xiW(ny,i);
                densiWbary(ny+1,i) = densiW(ny,i);
                xWbary(1:Nc,ny+1,i) = xW(1:Nc,ny,i);
            end 
        end
        if(Uny(ny+1,i) < 0)          
            lambdany(ny+1,i) = Kyybar(ny+1,i)*kr_N(Swtemp)/viscNtemp; 
        elseif(Uny(ny+1,i) > 0)
            lambdany(ny+1,i) = Kyybar(ny+1,i)*kr_N(Sw(ny,i))/viscN(ny,i);
            xiNbary(ny+1,i) = xiN(ny,i);
            densiNbary(ny+1,i) = densiN(ny,i);
            xNbary(1:Nc,ny+1,i) = xN(1:Nc,ny,i);
        else
            if(isDiriY(i,2) == 1) 
                lambdany(ny+1,i) = Kyybar(ny+1,i)*kr_N(Sw(ny,i))/viscN(ny,i);
                xiNbary(ny+1,i) = xiN(ny,i);
                densiNbary(ny+1,i) = densiN(ny,i);
                xNbary(1:Nc,ny+1,i) = xN(1:Nc,ny,i);
            end 
        end          
    end
        
    for i = 2 : ny
        for j = 1 : nx                
            if(Uwy(i,j) > 0)
                lambdawy(i,j) = Kyybar(i,j)*kr_W(Sw(i-1,j))/viscW(i-1,j);
                xiWbary(i,j) = xiW(i-1,j);
                densiWbary(i,j) = densiW(i-1,j);
                xWbary(1:Nc,i,j) = xW(1:Nc,i-1,j);              
            else
                lambdawy(i,j) = Kyybar(i,j)*kr_W(Sw(i,j))/viscW(i,j);
                xiWbary(i,j) = xiW(i,j);
                densiWbary(i,j) = densiW(i,j);
                xWbary(1:Nc,i,j) = xW(1:Nc,i,j);  
            end
            if(Uny(i,j) > 0)
                lambdany(i,j) = Kyybar(i,j)*kr_N(Sw(i-1,j))/viscN(i-1,j);
                xiNbary(i,j) = xiN(i-1,j);
                densiNbary(i,j) = densiN(i-1,j);
                xNbary(1:Nc,i,j) = xN(1:Nc,i-1,j);               
            else
                lambdany(i,j) = Kyybar(i,j)*kr_N(Sw(i,j))/viscN(i,j);
                xiNbary(i,j) = xiN(i,j);                
                densiNbary(i,j) = densiN(i,j);
                xNbary(1:Nc,i,j) = xN(1:Nc,i,j);        
            end
        end
    end

end

% The function computePresandVel() is to compute the new pressures and
% velocities
function [ Pw ] = computePres( model, Pw, lambdawx, lambdanx, lambdawy, lambdany, ...
    xiWbarx, xiNbarx, xiWbary, xiNbary, densiWbarx, densiWbary, densiNbarx, densiNbary, xWbarx, xNbarx, xWbary, xNbary, ...
        v, Cf, src, timestep )

    Nc = model.Nc;
    nx = model.nx;
    ny = model.ny;
    xs = model.xs;
    ys = model.ys;
    gravX = model.gravX;
    gravY = model.gravY;
    poro = model.poro;
    isDiriX = model.isDiriX;
    isDiriY = model.isDiriY;
    UwBdryX = model.UwBdryX;
    UwBdryY = model.UwBdryY;
    UnBdryX = model.UnBdryX;
    UnBdryY = model.UnBdryY;

    % define the coefficient matrices A and b
    A = sparse(nx*ny, nx*ny);
    b = zeros(nx*ny, 1);
    
    % rectify the directions of the velocities
    UwBdryX(1, 1:end) = -UwBdryX(1, 1:end);
    UwBdryY(1:end, 1) = -UwBdryY(1:end, 1);
    UnBdryX(1, 1:end) = -UnBdryX(1, 1:end);
    UnBdryY(1:end, 1) = -UnBdryY(1:end, 1);
       
    % compute the coefficiencies of the wetting phase pressure equations        
    for i = 2 : ny+1
        for j = 2 : nx+1
                
            r = (i-2)*nx + (j-1); % r is the index of the current cell

            xedge = xs(j) - xs(j-1); % x-edge of the current cell
            yedge = ys(i) - ys(i-1); % y-edge of the current cell            
                
            % ledge means the x-edge of the left cell 
            if(j ~= 2)
                ledge = xs(j-1) - xs(j-2);
            else
                ledge = 0;
            end
            
            % redge means the x-edge of the right cell
            if(j ~= nx+1)
                redge = xs(j+1) - xs(j);    
            else
                redge = 0;
            end
            
            % uedge means the y-edge of the up cell
            if(i ~= ny+1)
                uedge = ys(i+1) - ys(i);
            else
                uedge = 0;
            end
            
            % dedge means the y-edge of the down cell
            if(i ~= 2)
                dedge = ys(i-1) - ys(i-2);
            else
                dedge = 0;
            end
                
            A(r,r) = poro(j-1,i-1)*Cf(i-1,j-1)/timestep;
            b(r,1) = sum(v(1:Nc,i-1,j-1).*src(1:Nc,i-1,j-1)) + poro(j-1,i-1)*Cf(i-1,j-1)/timestep*Pw(i,j);
                                
            % some coefficients
            cotwx2 = lambdawx(i-1,j-1)*xiWbarx(i-1,j-1)*sum(xWbarx(1:Nc,i-1,j-1).*v(1:Nc,i-1,j-1));
            cotnx2 = lambdanx(i-1,j-1)*xiNbarx(i-1,j-1)*sum(xNbarx(1:Nc,i-1,j-1).*v(1:Nc,i-1,j-1));
            cotwx1 = lambdawx(i-1,j)*xiWbarx(i-1,j)*sum(xWbarx(1:Nc,i-1,j).*v(1:Nc,i-1,j-1));
            cotnx1 = lambdanx(i-1,j)*xiNbarx(i-1,j)*sum(xNbarx(1:Nc,i-1,j).*v(1:Nc,i-1,j-1));
            cotwy1 = lambdawy(i,j-1)*xiWbary(i,j-1)*sum(xWbary(1:Nc,i,j-1).*v(1:Nc,i-1,j-1));
            cotny1 = lambdany(i,j-1)*xiNbary(i,j-1)*sum(xNbary(1:Nc,i,j-1).*v(1:Nc,i-1,j-1));
            cotwy2 = lambdawy(i-1,j-1)*xiWbary(i-1,j-1)*sum(xWbary(1:Nc,i-1,j-1).*v(1:Nc,i-1,j-1));
            cotny2 = lambdany(i-1,j-1)*xiNbary(i-1,j-1)*sum(xNbary(1:Nc,i-1,j-1).*v(1:Nc,i-1,j-1));           
         
            % compute the coefficients of the equations                                    
            up = -2*(cotwy1+cotny1)/yedge/(yedge+uedge); % up means the coefficient of the pressure of the up cell
            down = -2*(cotwy2+cotny2)/yedge/(yedge+dedge); % down means the coefficient of the pressure of the down cell
            left = -2*(cotwx2+cotnx2)/xedge/(xedge+ledge); % left means the coefficient of the pressure of the left cell
            right = -2*(cotwx1+cotnx1)/xedge/(xedge+redge); % right means the coefficient of the pressure of the right cell            
            
            if((i == ny+1) && (isDiriY(j-1,2) == 0)) % Neuman boundary
                b(r,1) = b(r,1) - UwBdryY(j-1,2)*xiWbary(i,j-1)*sum(xWbary(1:Nc,i,j-1).*v(1:Nc,i-1,j-1))/yedge ...
                    - UnBdryY(j-1,2)*xiNbary(i,j-1)*sum(xNbary(1:Nc,i,j-1).*v(1:Nc,i-1,j-1))/yedge;
            elseif((i == ny+1) && (isDiriY(j-1,2) == 1)) % Dirichlet boundary
                b(r,1) = b(r,1) - up*Pw(i+1,j) - cotwy1*densiWbary(i,j-1)*gravY/yedge - cotny1*densiNbary(i,j-1)*gravY/yedge;
                A(r,r) = A(r,r) - up;
            else
                A(r,r+nx) = up;
                A(r,r) = A(r,r) - up;
                b(r,1) = b(r,1) - cotwy1*densiWbary(i,j-1)*gravY/yedge - cotny1*densiNbary(i,j-1)*gravY/yedge;
            end
                
            if((i == 2) && (isDiriY(j-1,1) == 0))
                b(r,1) = b(r,1) + UwBdryY(j-1,1)*xiWbary(i-1,j-1)*sum(xWbary(1:Nc,i-1,j-1).*v(1:Nc,i-1,j-1))/yedge ...
                    + UnBdryY(j-1,1)*xiNbary(i-1,j-1)*sum(xNbary(1:Nc,i-1,j-1).*v(1:Nc,i-1,j-1))/yedge;
            elseif((i == 2) && (isDiriY(j-1,1) == 1))
                b(r,1) = b(r,1) - down*Pw(i-1,j) + cotwy2*densiWbary(i-1,j-1)*gravY/yedge + cotny2*densiNbary(i-1,j-1)*gravY/yedge;
                A(r,r) = A(r,r) - down;
            else
                A(r,r-nx) = down;
                A(r,r) = A(r,r) - down;
                b(r,1) = b(r,1) + cotwy2*densiWbary(i-1,j-1)*gravY/yedge + cotny2*densiNbary(i-1,j-1)*gravY/yedge;
            end
                
            if((j == 2) && (isDiriX(1,i-1) == 0))
                b(r,1) = b(r,1) + UwBdryX(1,i-1)*xiWbarx(i-1,j-1)*sum(xWbarx(1:Nc,i-1,j-1).*v(1:Nc,i-1,j-1))/xedge ...
                    + UnBdryX(1,i-1)*xiNbarx(i-1,j-1)*sum(xNbarx(1:Nc,i-1,j-1).*v(1:Nc,i-1,j-1))/xedge;
            elseif((j == 2) && (isDiriX(1,i-1) == 1))
                b(r,1) = b(r,1) - left*Pw(i,j-1) + cotwx2*densiWbarx(i-1,j-1)*gravX/xedge + cotnx2*densiNbarx(i-1,j-1)*gravX/xedge;
                A(r,r) = A(r,r) - left;
            else
                A(r,r-1) = left;
                A(r,r) = A(r,r) - left;
                b(r,1) = b(r,1) + cotwx2*densiWbarx(i-1,j-1)*gravX/xedge + cotnx2*densiNbarx(i-1,j-1)*gravX/xedge;
            end
                
            if((j == nx+1) && (isDiriX(2,i-1) == 0))
                b(r,1) = b(r,1) - UwBdryX(2,i-1)*xiWbarx(i-1,j)*sum(xWbarx(1:Nc,i-1,j).*v(1:Nc,i-1,j-1))/xedge ...
                    - UnBdryX(2,i-1)*xiNbarx(i-1,j)*sum(xNbarx(1:Nc,i-1,j).*v(1:Nc,i-1,j-1))/xedge;  
            elseif((j == nx+1) && (isDiriX(2,i-1) == 1))
                b(r,1) = b(r,1) - right*Pw(i,j+1) - cotwx1*densiWbarx(i-1,j)*gravX/xedge - cotnx1*densiNbarx(i-1,j)*gravX/xedge;
                A(r,r) = A(r,r) - right;
            else
                A(r,r+1) = right;
                A(r,r) = A(r,r) - right;
                b(r,1) = b(r,1) - cotwx1*densiWbarx(i-1,j)*gravX/xedge - cotnx1*densiNbarx(i-1,j)*gravX/xedge;
            end
        end
    end
    
    % compute the new pressures of wetting phase
    xx = A\b;
    for i = 2 : ny+1
        for j = 2 : nx+1
            Pw(i,j) = xx((i-2)*nx+(j-1),1);
        end
    end
    
end

% The function computeVel() is to compute the new velocities
function [ Uwx, Uwy, Unx, Uny ] = computeVel( model, Pw, lambdawx, lambdanx, lambdawy, lambdany, ...
    densiWbarx, densiWbary, densiNbarx, densiNbary )

    nx = model.nx;
    ny = model.ny;
    xs = model.xs;
    ys = model.ys;
    gravX = model.gravX;
    gravY = model.gravY;
    isDiriX = model.isDiriX;
    isDiriY = model.isDiriY;
    UwBdryX = model.UwBdryX;
    UwBdryY = model.UwBdryY;
    UnBdryX = model.UnBdryX;
    UnBdryY = model.UnBdryY;
    Uwx = zeros(ny, nx+1);
    Uwy = zeros(ny+1, nx);
    Unx = zeros(ny, nx+1);
    Uny = zeros(ny+1, nx);
    
    % rectify the directions of the velocities
    UwBdryX(1, 1:end) = -UwBdryX(1, 1:end);
    UwBdryY(1:end, 1) = -UwBdryY(1:end, 1);
    UnBdryX(1, 1:end) = -UnBdryX(1, 1:end);
    UnBdryY(1:end, 1) = -UnBdryY(1:end, 1);
    
    % compute the total new velocities in x direction  
    for i = 1 : ny
        if(isDiriX(1,i) == 1)
            Uwx(i,1) = -lambdawx(i,1)*((Pw(i+1,2)-Pw(i+1,1))*2/(xs(2)-xs(1)) - densiWbarx(i,1)*gravX);
            Unx(i,1) = -lambdanx(i,1)*((Pw(i+1,2)-Pw(i+1,1))*2/(xs(2)-xs(1)) - densiNbarx(i,1)*gravX);
        else
            Uwx(i,1) = UwBdryX(1,i);
            Unx(i,1) = UnBdryX(1,i);
        end
    end
    
    for i = 1 : ny
        if(isDiriX(2,i) == 1)
            Uwx(i,nx+1) = -lambdawx(i,nx+1)*((Pw(i+1,nx+2)-Pw(i+1,nx+1))*2/(xs(nx+1)-xs(nx)) - densiWbarx(i,nx+1)*gravX);
            Unx(i,nx+1) = -lambdanx(i,nx+1)*((Pw(i+1,nx+2)-Pw(i+1,nx+1))*2/(xs(nx+1)-xs(nx)) - densiNbarx(i,nx+1)*gravX);
        else
            Uwx(i,nx+1) = UwBdryX(2,i);
            Unx(i,nx+1) = UnBdryX(2,i);
        end
    end
    
    for i = 1 : ny
        for j = 2 : nx
            Uwx(i,j) = -lambdawx(i,j)*((Pw(i+1,j+1)-Pw(i+1,j))*2/(xs(j+1)-xs(j-1)) - densiWbarx(i,j)*gravX);
            Unx(i,j) = -lambdanx(i,j)*((Pw(i+1,j+1)-Pw(i+1,j))*2/(xs(j+1)-xs(j-1)) - densiNbarx(i,j)*gravX);
        end
    end
    
    % compute the total new velocities in y direction
    for j = 1 : nx
        if(isDiriY(j,1) == 1)
            Uwy(1,j) = -lambdawy(1,j)*((Pw(2,j+1)-Pw(1,j+1))*2/(ys(2)-ys(1)) - densiWbary(1,j)*gravY);
            Uny(1,j) = -lambdany(1,j)*((Pw(2,j+1)-Pw(1,j+1))*2/(ys(2)-ys(1)) - densiNbary(1,j)*gravY);
        else
            Uwy(1,j) = UwBdryY(j,1);
            Uny(1,j) = UnBdryY(j,1);
        end
    end
    
    for j = 1 : nx
        if(isDiriY(j,2) == 1)
            Uwy(ny+1,j) = -lambdawy(ny+1,j)*((Pw(ny+2,j+1)-Pw(ny+1,j+1))*2/(ys(ny+1)-ys(ny)) - densiWbary(ny+1,j)*gravY);
            Uny(ny+1,j) = -lambdany(ny+1,j)*((Pw(ny+2,j+1)-Pw(ny+1,j+1))*2/(ys(ny+1)-ys(ny)) - densiNbary(ny+1,j)*gravY);
        else
            Uwy(ny+1,j) = UwBdryY(j,2);
            Uny(ny+1,j) = UnBdryY(j,2);
        end
    end

    for i = 2 : ny
        for j = 1 : nx
            Uwy(i,j) = -lambdawy(i,j)*((Pw(i+1,j+1)-Pw(i,j+1))*2/(ys(i+1)-ys(i-1)) - densiWbary(i,j)*gravY);
            Uny(i,j) = -lambdany(i,j)*((Pw(i+1,j+1)-Pw(i,j+1))*2/(ys(i+1)-ys(i-1)) - densiNbary(i,j)*gravY);
        end
    end
    
end

% The function computezandxi() is to compute the new mole fraction of
% each component and new mole density of each cell.
function [ z ] = computez( model, timestep, xi, z, src, Uwx, Uwy, Unx, Uny, xiWbarx, ...
    xiWbary, xiNbarx, xiNbary, xWbarx, xWbary, xNbarx, xNbary )

    Nc = model.Nc;
    nx = model.nx;
    ny = model.ny;
    xs = model.xs;
    ys = model.ys;
    poro = model.poro;
    
    for i = 1 : ny
        for j = 1 : nx
            xic = zeros(Nc,1); % molar density of each component
            for m = 1 : Nc
                right = Uwx(i,j+1)*xiWbarx(i,j+1)*xWbarx(m,i,j+1) + Unx(i,j+1)*xiNbarx(i,j+1)*xNbarx(m,i,j+1);
                left = Uwx(i,j)*xiWbarx(i,j)*xWbarx(m,i,j) + Unx(i,j)*xiNbarx(i,j)*xNbarx(m,i,j);
                up = Uwy(i+1,j)*xiWbary(i+1,j)*xWbary(m,i+1,j) + Uny(i+1,j)*xiNbary(i+1,j)*xNbary(m,i+1,j);
                down = Uwy(i,j)*xiWbary(i,j)*xWbary(m,i,j) + Uny(i,j)*xiNbary(i,j)*xNbary(m,i,j);
                div = (right-left)/(xs(j+1)-xs(j)) + (up-down)/(ys(i+1)-ys(i));
                xic(m) = (src(m,i,j)-div)*timestep/poro(j,i) + z(m,i,j)*xi(i,j);
                if((xic(m)<0) || ~isreal(xic(m)))
                    xic(m)
                    error('Please tune the time step.');
                end
            end
            xinew = sum(xic(1:Nc));
            z(1:Nc,i,j) = xic(1:Nc)/xinew;
        end
    end  
    
end