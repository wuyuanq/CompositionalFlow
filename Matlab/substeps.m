% This is the sub function. It realizes the steps 4 to 6.

% Input parameters:
% model: the model
% x: the mole fraction of each component in the liquid phase
% y: the mole fraction of each component in the gas phase
% P: the pressure in the cell

% Return values:
% criteria: the convergent criteria
% K: the equilibrium ratio
% ZL: the compressibility of liquid phase
% ZG: the compressibility of gas phase
% am: parameter
% bm: parameter
% al: parameter
% ag: parameter
% bl: parameter
% bg: parameter
% AL: parameter
% AG: parameter
% BL: parameter
% BG: parameter
% xiL: the molar density of liquid phase, unit: mol/m^3
% xiG: the molar density of gas phase, unit: mol/m^3
% rhoL: the mass density of liguid phase, unit: kg/m^3
% rhoG: the mass density of gas phase, unit: kg/m^3
% CfL: the compressibility coefficient in liquid phase
% CfG: the compressibility coefficient in gas phase

% Author: Yuanqing Wu. Email: wuyuanq@gmail.com
% Last edited on January 23rd, 2014

function [ criteria, K, ZL, ZG, am, bm, al, ag, bl, bg, AL, AG, BL, BG, xiL, xiG, rhoL, rhoG, CfL, CfG ] ...
    = substeps( model, x, y, P )

    % simplify the notations of model
    Nc = model.Nc;   
            
    [ ZL, ZG, am, bm, al, ag, bl, bg, AL, AG, BL, BG, xiL, xiG, rhoL, ...
        rhoG, CfL, CfG, phil, phig  ] = RST_fugacitycoef( model, x, y, P );

    fl = zeros(Nc, 1);
    fg = zeros(Nc, 1);
    K = zeros(Nc, 1);

    for i = 1 : Nc
        fl(i) = x(i)*phil(i)*P;
        fg(i) = y(i)*phig(i)*P;
        K(i) = phil(i)/phig(i);
    end 

    criteria = 0;
    nn = 0; % the number of existing components
    for i = 1 : Nc
        if((x(i) ~= 0)||(y(i) ~= 0)) % gurantee that component i exists
            criteria = criteria + (log(fg(i)/fl(i)))^2;
            nn = nn + 1;
        end 
    end 
    criteria = criteria / nn;

end