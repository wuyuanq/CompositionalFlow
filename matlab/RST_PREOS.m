% The function RST_PREOS is to use the PR-EOS to get the compressibility
% factors and so on of the liquid and gas phases.

% Input parameters:
% model: the model
% x: the mole fraction in liquid phase of each component
% y: the mole fraction in gas phase of each component
% P: the pressure

% Return value:
% ZL: the compressibility factor of the liquid phase
% ZG: the compressibility factor of the gas phase 
% am, bm, al, ag, bl, bg, AL, AG, BL, BG: some coefficients
% xiL: the molar density in liquid phase, unit: mol/m^3
% xiG: the molar density in gas phase, unit: mol/m^3
% rhoL: the mass density of liguid phase, unit: kg/m^3
% rhoG: the mass density of gas phase, unit: kg/m^3
% CfL: the compressibility coefficient in liquid phase
% CfG: the compressibility coefficient in gas phase

% Reference: Chap 7. in <<Reservoir Simulation: Mathematical Techniques in Oil Recovery>>

% Author: Yuanqing Wu. Email: wuyuanq@gmail.com
% Last edited on January 23rd, 2014

function [ ZL, ZG, am, bm, al, ag, bl, bg, AL, AG, BL, BG, xiL, xiG, rhoL, rhoG, CfL, CfG ] ...
    = RST_PREOS( model, x, y, P )

    % simplify the notations of model
    Nc = model.Nc;
    T = model.T;
    Tc = model.ct;
    Pc = model.cp;
    delta = model.delta;
    omega = model.af;
    
    R = 8.314; % gas constant, unit: J/(mol*K)
    
    am = zeros(Nc,1);
    bm = zeros(Nc,1);
    
    lambda = zeros(Nc, 1);
    for i = 1 : Nc
        lambda(i) = 0.37464 + 1.5423*omega(i) - 0.26992*omega(i)^2;
    end
    
    alpha = zeros(Nc, 1);
    for i = 1 : Nc
        alpha(i) = (1+lambda(i)*(1-sqrt(T/Tc(i))))^2;
    end
    
    for i = 1 : Nc
        am(i) = 0.45724*alpha(i)*R^2*Tc(i)^2/Pc(i); 
        bm(i) = 0.077796*R*Tc(i)/Pc(i);
    end
           
    al = 0;
    for i = 1 : Nc
        for j = 1 : Nc
            al = al + x(i)*x(j)*(1-delta(i,j))*sqrt(am(i)*am(j));
        end
    end
        
    ag = 0;
    for i = 1 : Nc
        for j = 1 : Nc
            ag = ag + y(i)*y(j)*(1-delta(i,j))*sqrt(am(i)*am(j));
        end
    end
       
    bl = 0;
    for i = 1 : Nc
        bl = bl + x(i)*bm(i);
    end
        
    bg = 0;
    for i = 1 : Nc
        bg = bg + y(i)*bm(i);
    end
        
    [ AL, BL, ZL, xiL, rhoL, CfL ] = computeLiquidPhase( model, x, P, al, bl );
        
    [ AG, BG, ZG, xiG, rhoG, CfG ] = computeGasPhase( model, y, P, ag, bg );
    
end

function [ AL, BL, ZL, xiL, rhoL, CfL ] = computeLiquidPhase(model, x, P, al, bl )

    % simplify the notations of model
    Nc = model.Nc;
    T = model.T;
    Tc = model.ct;
    mw = model.mw;
    omega = model.af;
    
    R = 8.314; % gas constant, unit: J/(mol*K)
    
    AL = al*P/(R*T)^2;
    BL = bl*P/(R*T);

    Z = roots([1 -(1-BL) (AL-3*BL^2-2*BL) -(AL*BL-BL^2-BL^3)]); % compressibility factor   

    ZRP = [];
    for i = 1 : 3
        if(isreal(Z(i)) && (Z(i)>0) && (Z(i)>BL))
            ZRP = [ZRP Z(i)];   
        end
    end
    ZL = min(ZRP);
    if(size(ZL) == 0)
        ZL
        error('ZL doesn''t have a reasonable value.')
    end
        
    xiL_bs = P/(R*T*ZL); % the molar density in liquid phase, before shifting
    
    % shift the liquid volume and correct the molar density in liquid phase
    % reference: <<GENERALIZED LIQUID VOLUME SHIFTS FOR THE PENGROBINSON
    % EQUATION OF STATE FOR C1 TO C8 HYDROCARBONS>>
    c = zeros(Nc,1); % volume correct item for each component
    C1toC2 = zeros(Nc,1);
    for i = 1 : Nc
        C1toC2(i) = 110.07*omega(i)^4 - 83.807*omega(i)^3 + 18.926*omega(i)^2 - 1.6348*omega(i) - 0.0066;
    end
    C2 = 2.013645*1.D-3; % unit: m^3/kg
    C3 = 0.89;
    for i = 1 : Nc
        c(i) = C1toC2(i)*C2 + C2*(T/Tc(i)-C3)^2;
    end
    ctotal = 0;
    for i = 1 : Nc 
        ctotal = ctotal + x(i)*c(i)*mw(i); % unit: m^3/mol
    end
    xiL = 1/(1/xiL_bs + ctotal); % unit: m^3/mol
    
    % the mass density
    rhoL = 0;
    for m = 1 : Nc
        rhoL = rhoL + x(m)*mw(m);
    end
    rhoL = rhoL*xiL;
    
    % compute the compressibility coefficient in liquid phase 
    deri_AL_p = al/(R*T)^2;
    deri_BL_p = bl/(R*T);
    deri_ZL_p = -(deri_BL_p*ZL^2+(deri_AL_p-2*(1+3*BL)*deri_BL_p)*ZL - ...
        (deri_AL_p*BL+(AL-2*BL-3*BL^2)*deri_BL_p))/(3*ZL^2-2*(1-BL)*ZL+(AL-2*BL-3*BL^2));
    deri_xiL_p_bs = 1/(R*T*ZL) - P/(R*T*ZL^2)*deri_ZL_p;
    deri_xiL_p = (1/(1+ctotal*xiL_bs)^2)*deri_xiL_p_bs;
    CfL = deri_xiL_p/xiL;

end

function [AG, BG, ZG, xiG, rhoG, CfG ] = computeGasPhase(model, y, P, ag, bg )

    % simplify the notations of model
    Nc = model.Nc;
    T = model.T;
    mw = model.mw;
    
    R = 8.314; % gas constant, unit: J/(mol*K)

    % compute the compressibility factor in gas phase
    AG = ag*P/(R*T)^2;
    BG = bg*P/(R*T);
    
    Z = roots([1 -(1-BG) (AG-3*BG^2-2*BG) -(AG*BG-BG^2-BG^3)]);
    
    ZR = [];
    for i = 1 : 3
        if isreal(Z(i))
            ZR = [ZR Z(i)];   
        end
    end   
    ZG = max(ZR); 
    if(size(ZG) == 0)
        ZG
        error('ZG doesn''t have a reasonable value. ')
    end
      
    xiG = P/(R*T*ZG); % the molar density in gas phase
    
    % the mass density
    rhoG = 0;
    for m = 1 : Nc
        rhoG = rhoG + y(m)*mw(m);
    end
    rhoG = rhoG*xiG;
    
    % compute the compressibility coefficient in gas phase 
    deri_AG_p = ag/(R*T)^2;
    deri_BG_p = bg/(R*T);
    deri_ZG_p = -(deri_BG_p*ZG^2+(deri_AG_p-2*(1+3*BG)*deri_BG_p)*ZG - ...
        (deri_AG_p*BG+(AG-2*BG-3*BG^2)*deri_BG_p))/(3*ZG^2-2*(1-BG)*ZG+(AG-2*BG-3*BG^2));
    CfG = 1/P - 1/ZG*deri_ZG_p;
    
end
