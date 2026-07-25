% The function RST_comprefac is to compute the compressibility factor of the cell

% Input parameters:
% model: the model
% x: the mole fraction of each component in the liquid phase
% y: the mole fraction of each component in the gas phase
% P: the pressure in the cell
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
% sL: the saturation

% Return values:
% Cf: the compressibility factor

% Author: Yuanqing Wu. Email: wuyuanq@gmail.com
% Last edited on January 23rd, 2014

function [ Cf ] = RST_comprefac( model, x, y, P, ZL, ZG, am, bm, al, ag, bl, bg, ...
    AL, AG, BL, BG, xiL, xiG, sL )

    R = 8.314;

    Nc = model.Nc;
    T = model.T;
    delta = model.delta;

    A = zeros(Nc, Nc);
    b = zeros(Nc, 1);

    tempsumg = zeros(Nc, 1);
    tempsuml = zeros(Nc, 1);
    paragy = zeros(Nc, 1);
    paralx = zeros(Nc, 1);
    parAGy = zeros(Nc, 1);
    parALx = zeros(Nc, 1);
    parBGy = zeros(Nc, 1);
    parBLx = zeros(Nc, 1);
    parZGy = zeros(Nc, 1);
    parZLx = zeros(Nc, 1);
    par2gy = zeros(Nc, 1);
    par2lx = zeros(Nc, 1);
    par31gy = zeros(Nc, 1);
    par31lx = zeros(Nc, 1);
    par33gy = zeros(Nc, 1);
    par33lx = zeros(Nc, 1);
    rightgy = zeros(Nc, 1);
    rightlx = zeros(Nc, 1);
    par1gy = zeros(Nc, Nc);
    par1lx = zeros(Nc, Nc);
    par32gy = zeros(Nc, Nc);
    par32lx = zeros(Nc, Nc);
    par3gy = zeros(Nc, Nc);
    par3lx = zeros(Nc, Nc);
    parfgy = zeros(Nc, Nc);
    parflx = zeros(Nc, Nc);
    parfgn = zeros(Nc, Nc);
    parfln = zeros(Nc, Nc);
    par1gp = zeros(Nc, 1);
    par1lp = zeros(Nc, 1);
    par2gp = zeros(Nc, 1);
    par2lp = zeros(Nc, 1);
    par3gp = zeros(Nc, 1);
    par3lp = zeros(Nc, 1);
    parfgp = zeros(Nc, 1);
    parflp = zeros(Nc, 1);
    parnlp = zeros(Nc, 1);
    parngp = zeros(Nc, 1);
    parAGn = zeros(Nc, 1);
    parALn = zeros(Nc, 1);
    parBGn = zeros(Nc, 1);
    parBLn = zeros(Nc, 1);
    parZGn = zeros(Nc, 1);
    parZLn = zeros(Nc, 1);

    for i = 1 : Nc
        tempsumg(i) = 0;
        for j = 1 : Nc
            tempsumg(i) = tempsumg(i) + y(j)*(1-delta(i,j))*sqrt(am(i)*am(j));
        end 
        tempsumg(i) = 2*tempsumg(i);
    end 

    for i = 1 : Nc
        tempsuml(i) = 0;
        for j = 1 : Nc
            tempsuml(i) = tempsuml(i) + x(j)*(1-delta(i,j))*sqrt(am(i)*am(j));
        end 
        tempsuml(i) = 2*tempsuml(i);
    end 

    parZGAG = (BG-ZG)/(3*ZG^2-2*ZG*(1-BG)+(AG-3*BG^2-2*BG));
    parZLAL = (BL-ZL)/(3*ZL^2-2*ZL*(1-BL)+(AL-3*BL^2-2*BL));

    parZGBG = (-ZG^2+2*(3*BG+1)*ZG+(AG-2*BG-3*BG^2))/ ...
        (3*ZG^2-2*(1-BG)*ZG+(AG-3*BG^2-2*BG));
    parZLBL = (-ZL^2+2*(3*BL+1)*ZL+(AL-2*BL-3*BL^2))/ ...
        (3*ZL^2-2*(1-BL)*ZL+(AL-3*BL^2-2*BL));

    % suppose Vf = 1, then ng = xiG*(1-sL)*Vf = xiG*(1-sL)
    ng = xiG*(1-sL);
    nl = xiL*sL;

    for i = 1 : Nc

        paragy(i) = 0;
        for j = 1 : Nc
            if(j ~= i) 
                paragy(i) = paragy(i) + y(j)*(1-delta(i,j))*sqrt(am(i)*am(j));
            end 
        end 
        paragy(i) = 2*paragy(i);

        paralx(i) = 0;
        for j = 1 : Nc
            if(j ~= i) 
                paralx(i) = paralx(i) + x(j)*(1-delta(i,j))*sqrt(am(i)*am(j));
            end 
        end 
        paralx(i) = 2*paralx(i);

        parAGy(i) = P/(R*T)^2*paragy(i);
        parALx(i) = P/(R*T)^2*paralx(i);

        parBGy(i) = P*bm(i)/R/T;
        parBLx(i) = P*bm(i)/R/T;

        parZGy(i) = parZGAG*P/(R*T)^2*paragy(i) + parZGBG*P*bm(i)/R/T;

        parZLx(i) = parZLAL*P/(R*T)^2*paralx(i) + parZLBL*P*bm(i)/R/T;

        par2gy(i) = 1/(ZG-BG)*(parZGy(i)-parBGy(i));
        par2lx(i) = 1/(ZL-BL)*(parZLx(i)-parBLx(i));

        par31gy(i) = parAGy(i)/(2*sqrt(2.0)*BG)-AG/(2*sqrt(2.0)*BG^2)*parBGy(i);
        par31lx(i) = parALx(i)/(2*sqrt(2.0)*BL)-AL/(2*sqrt(2.0)*BL^2)*parBLx(i);

        par33gy(i) = ((parZGy(i)+2.414*parBGy(i))*(ZG-0.414*BG)-(ZG+2.414*BG)* ...
            (parZGy(i)-0.414*parBGy(i)))/(ZG+2.414*BG)/(ZG-0.414*BG);
        par33lx(i) = ((parZLx(i)+2.414*parBLx(i))*(ZL-0.414*BL)-(ZL+2.414*BL)* ...
            (parZLx(i)-0.414*parBLx(i)))/(ZL+2.414*BL)/(ZL-0.414*BL);

        rightgy(i) = bm(i)/bg*(ZG-1)-log(ZG-BG) - AG/(2*sqrt(2.0)* ...
            BG)*(tempsumg(i)-bm(i)/bg)*log((ZG+2.414*BG)/(ZG-0.414*BG));
        rightlx(i) = bm(i)/bl*(ZL-1)-log(ZL-BL) - AL/(2*sqrt(2.0)* ...
            BL)*(tempsuml(i)-bm(i)/bl)*log((ZL+2.414*BL)/(ZL-0.414*BL));

    end 

    for j = 1 : Nc
        for i = 1 : Nc

            par1gy(i,j) = bm(i)/bg*parZGy(j) - bm(i)*bm(j)/bg^2*(ZG-1);
            par1lx(i,j) = bm(i)/bl*parZLx(j) - bm(i)*bm(j)/bl^2*(ZL-1);

            par32gy(i,j) = 2/ag*(1-delta(i,j))*sqrt(am(i)*am(j)) -paragy(j)/ ...
                ag^2*tempsumg(i) + bm(i)*bm(j)/bg^2;

            par32lx(i,j) = 2/al*(1-delta(i,j))*sqrt(am(i)*am(j)) -paralx(j)/ ...
                al^2*tempsuml(i) + bm(i)*bm(j)/bl^2;

            par3gy(i,j) = AG/(2*sqrt(2.0)*BG)*(tempsumg(i)/ag-bm(i)/bg)*par33gy(j) ...
                + AG/(2*sqrt(2.0)*BG)*log((ZG+2.414*BG)/(ZG-0.414*BG))* ...
                par32gy(i,j)+(tempsumg(i)/ag-bm(i)/bg)*log((ZG+2.414*BG)/ ...
                (ZG-0.414*BG))*par31gy(j);

            par3lx(i,j) = AL/(2*sqrt(2.0)*BL)*(tempsuml(i)/al-bm(i)/bl)*par33lx(j) ...
                + AL/(2*sqrt(2.0)*BL)*log((ZL+2.414*BL)/(ZL-0.414*BL))* ...
                par32lx(i,j)+(tempsuml(i)/al-bm(i)/bl)*log((ZL+2.414*BL)/ ...
                (ZL-0.414*BL))*par31lx(j);

            parfgy(i,j) = y(i)*P*exp(rightgy(i))*(par1gy(i,j)-par2gy(j)- par3gy(i,j));
            parflx(i,j) = x(i)*P*exp(rightlx(i))*(par1lx(i,j)-par2lx(j)- par3lx(i,j));

            if(i == j)
                parfgy(i,j) = parfgy(i,j) + P*exp(rightgy(i));
                parflx(i,j) = parflx(i,j) + P*exp(rightlx(i));
            end 

            parfgn(i,j) = parfgy(i,j)/ng;
            parfln(i,j) = parflx(i,j)/nl;
        end 
    end 

    parZGp = parZGAG*ag/(R*T)^2 + parZGBG*bg/R/T;
    parZLp = parZLAL*al/(R*T)^2 + parZLBL*bl/R/T;

    parAGp = ag/(R*T)^2;
    parALp = al/(R*T)^2;

    parBGp = bg/R/T;
    parBLp = bl/R/T;

    par31gp = parAGp/(2*sqrt(2.0)*BG) - AG*parBGp/(2*sqrt(2.0)*BG^2);
    par31lp = parALp/(2*sqrt(2.0)*BL) - AL*parBLp/(2*sqrt(2.0)*BL^2);

    par33gp = ((parZGp-0.414*parBGp)*(ZG+2.414*BG)-(ZG-0.414*BG)* ...
        (parZGp+2.414*parBGp))/((ZG+2.414*BG)*(ZG-0.414*BG));
    par33lp = ((parZLp-0.414*parBLp)*(ZL+2.414*BL)-(ZL-0.414*BL)* ...
        (parZLp+2.414*parBLp))/((ZL+2.414*BL)*(ZL-0.414*BL));

    for i = 1 : Nc

        par1gp(i) = bm(i)/bg*parZGp;
        par1lp(i) = bm(i)/bl*parZLp;

        par2gp(i) = 1/(ZG-BG)*(parZGp-bg/R/T);
        par2lp(i) = 1/(ZL-BL)*(parZLp-bl/R/T);

        par3gp(i) = AG/(2*sqrt(2.0)*BG)*(tempsumg(i)/ag-bm(i)/bg)*par33gp + ...
            (tempsumg(i)/ag-bm(i)/bg)*log((ZG+2.414*BG)/(ZG-0.414*BG))*par31gp;
        par3lp(i) = AL/(2*sqrt(2.0)*BL)*(tempsuml(i)/al-bm(i)/bl)*par33lp + ...
            (tempsuml(i)/al-bm(i)/bl)*log((ZL+2.414*BL)/(ZL-0.414*BL))*par31lp;

        parfgp(i) = y(i)*exp(rightgy(i)) + y(i)*P*exp(rightgy(i))*(par1gp(i)-par2gp(i)-par3gp(i));
        parflp(i) = x(i)*exp(rightlx(i)) + x(i)*P*exp(rightlx(i))*(par1lp(i)-par2lp(i)-par3lp(i));

    end 

    for i = 1 : Nc
        for j = 1 : Nc
            A(i,j) = parfgn(i,j) + parfln(i,j);
        end 
    end 

    for i = 1 : Nc
        b(i) = parfgp(i) - parflp(i);
    end 

    solx = A\b; 
    
    for i = 1 : Nc
        parnlp(i) = solx(i);
        parngp(i) = -solx(i);
    end 

    for i = 1 : Nc

        parAGn(i) = ag/(R*T)^2/parngp(i)+P/(R*T)^2/ng*paragy(i);
        parALn(i) = al/(R*T)^2/parnlp(i)+P/(R*T)^2/nl*paralx(i);

        parBGn(i) = bg/R/T/parngp(i)+P/R/T*bm(i)/ng;
        parBLn(i) = bl/R/T/parnlp(i)+P/R/T*bm(i)/nl;

        parZGn(i) = parZGAG*parAGn(i)+parZGBG*parBGn(i);
        parZLn(i) = parZLAL*parALn(i)+parZLBL*parBLn(i);

    end 

    parZGPtotaln = 0;
    for i = 1 : Nc
        if(~isnan(parZGn(i)))
            parZGPtotaln = parZGPtotaln + parZGn(i)*parngp(i);
        end 
    end 
    parZGPtotaln = parZGPtotaln + parZGp;

    parZLPtotaln = 0;
    for i = 1 : Nc
        if(~isnan(parZLn(i))) 
            parZLPtotaln = parZLPtotaln + parZLn(i)*parnlp(i);
        end 
    end 
    parZLPtotaln = parZLPtotaln + parZLp;

    parngptotaln = 0;
    for i = 1 : Nc
        parngptotaln = parngptotaln + parngp(i);
    end 

    parnlptotaln = 0;
    for i = 1 : Nc
        parnlptotaln = parnlptotaln + parnlp(i);
    end 

    parvfp = -R*T/P^2*(ZG*ng+ZL*nl)+R*T/P*(ng*parZGPtotaln+ZG*parngptotaln + ...
        nl*parZLPtotaln+ZL*parnlptotaln);

    % since Vf=1, Cf = -parvfp/Vf = -parvfp
    Cf = -parvfp;
    
end

