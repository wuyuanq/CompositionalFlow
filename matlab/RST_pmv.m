% This is the RST_pmv() function which computes the partial molar volume of
% each component.

% Input parameters:
% model: the model
% comp: the component serial number
% P: the pressure
% z: the mole fraction of each component
% Vold: the old volume in the cell
% K: equilibrium ratio

% Return values:
% v: the partial molar volume of each component, unit: m^3/mol

% Author: Yuanqing Wu. Email: wuyuanq@gmail.com
% Last edited on November 23rd, 2012

function [ v ] = RST_pmv( model, comp, P, z, Vold, K )

    % simplify the notations of model
    Nc = model.Nc;
    
    MAXTIME = 1000;
    
    if(abs(z(comp)) < 1.D-12) 
        v = 0.0;
        return
    end 

    x = zeros(Nc, 1);
    y = zeros(Nc, 1);
    znew = zeros(Nc, 1);

    incmole = z(comp)*1.D-5; % need adjust according to the case
    molenew = 1+incmole;

    for i = 1 : Nc
        znew(i) = z(i)/molenew;
    end 
    znew(comp) = (z(comp)+incmole)/molenew;

    for t = 1 : MAXTIME

        h0 = 0;
        for i = 1 : Nc
            h0 = h0 + K(i)*znew(i);
        end 
        h0 = h0 - 1;
        if(h0 <= 0) 
            for i = 1 : Nc
                x(i) = znew(i);
            end 
            yt = 0;
            for i = 1 : Nc
                yt = yt + K(i)*znew(i);
            end 
            for i = 1 : Nc
                y(i) = K(i)*znew(i)/yt;
            end 
            [ K, v, jump ] = submoduleforpmv( model, x, y, znew, P, molenew, Vold, incmole );
            if(jump)
                break
            end
            continue
        end 

        h1 = 0;
        for i = 1 : Nc
            h1 = h1 + znew(i)/K(i);
        end 
        h1 = 1 - h1;
        if(h1 >= 0) 
            for i = 1 : Nc
                y(i) = znew(i);
            end 
            xt = 0;
            for i = 1 : Nc
                xt = xt + znew(i)/K(i);
            end 
            for i = 1 : Nc
                x(i) = znew(i)/K(i)/xt;
            end 
            [ K, v, jump ] = submoduleforpmv( model, x, y, znew, P, molenew, Vold, incmole );
            if(jump)
                break
            end
            continue
        end 

        beta = 0.5;
        lp = 0;
        rp = 1;
        while(true)
            hbeta = 0;
            for i = 1 : Nc
                hbeta = hbeta + ((K(i)-1)*znew(i))/(1+beta*(K(i)-1));
            end 
            if(abs(hbeta) < 1.D-12) 
                break
            end 

            if(hbeta > 0) 
                lp = beta;
            else
                rp = beta;
            end 
            beta = (lp + rp)/2;

        end 

        for i = 1 : Nc
            x(i) = znew(i)/(1+(K(i)-1)*beta);
            y(i) = K(i)*x(i);
        end 

        [ K, v, jump ] = submoduleforpmv( model, x, y, znew, P, molenew, Vold, incmole );
        if(jump)
            break
        end
    end 
    
end

function [ K, v, jump ] = submoduleforpmv( model, x, y, znew, P, molenew, Vold, incmole )

    Nc = model.Nc;
    jump =false;
    v = 0;
    
    [ criteria, K, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, xiL, xiG, rhoL, rhoG, ~, ~ ] ...
        = substeps( model, x, y, P );

    if(rhoL < rhoG)
        [ criteria, K, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, xiL, xiG, ~, ~, ~, ~ ] ...
            = substeps( model, y, x, P );
    end 

    if(criteria < 1.D-12) 
        for i = 1 : Nc
            if(znew(i) ~= 0) 
                nn1 = i;
                break
            end 
        end 
        for i = nn1+1 : Nc
            if(znew(i) ~= 0) 
                nn2 = i;
                break
            end 
        end 

        phasemole(1) = (y(nn2)*znew(nn1)-y(nn1)*znew(nn2))/(x(nn1)*y(nn2)-x(nn2)*y(nn1))*molenew;
        phasemole(2) = (x(nn1)*znew(nn2)-x(nn2)*znew(nn1))/(x(nn1)*y(nn2)-x(nn2)*y(nn1))*molenew;
        Vnew = phasemole(1)/xiL + phasemole(2)/xiG;
        v = (Vnew-Vold)/incmole;

        jump = true;
    end
    
end
