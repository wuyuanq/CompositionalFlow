% The function RST_fugacitycoef is to compute the fugacity coefficients of
% liquid and gas phases

% Input parameters:
% model: the model
% x: the mole fraction of each component in the liquid phase
% y: the mole fraction of each component in the gas phase
% P: the pressure in the cell

% Return values:
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
% xiL: the molar density of liquid phase
% xiG: the molar density of gas phase
% rhoL: the mass density of liquid phase
% rhoG: the mass density of gas phase
% CfL: the compressibility factor of liquid phase
% CfG: the compressibility factor of gas phase
% phil: the fugacity coefficient of liquid phase
% phig: the fugacity coefficient of gas phase

% Author: Yuanqing Wu. Email: wuyuanq@gmail.com
% Last edited on January 23rd, 2014

function [ ZL, ZG, am, bm, al, ag, bl, bg, AL, AG, BL, BG, xiL, xiG, rhoL, ...
    rhoG, CfL, CfG, phil, phig  ] = RST_fugacitycoef( model, x, y, P )

    Nc = model.Nc;
    delta = model.delta;
    phil = zeros(Nc, 1);
    phig = zeros(Nc, 1);
        
    [ ZL, ZG, am, bm, al, ag, bl, bg, AL, AG, BL, BG, xiL, xiG, rhoL, rhoG, CfL, CfG ] ...
        = RST_PREOS( model, x, y, P );

    for i = 1 : Nc
        CL = 0;
        for j = 1 : Nc
            CL = CL + x(j)*(1-delta(i,j))*sqrt(am(i)*am(j));
        end 
        CL = CL * 2.0/al;

        CG = 0;
        for j = 1 : Nc
            CG = CG + y(j)*(1-delta(i,j))*sqrt(am(i)*am(j));
        end 
        CG = CG * 2.0/ag;

        phil(i) = exp(bm(i)/bl*(ZL-1)-log(ZL-BL) - AL/(2*sqrt(2.0)* ...
            BL)*(CL-bm(i)/bl)*log((ZL+2.414*BL)/(ZL-0.414*BL)));
        phig(i) = exp(bm(i)/bg*(ZG-1)-log(ZG-BG) - AG/(2*sqrt(2.0)* ...
            BG)*(CG-bm(i)/bg)*log((ZG+2.414*BG)/(ZG-0.414*BG)));
    end 

end