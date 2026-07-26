% The function RST_viscosity is to compute the viscosity of the cell

% Input parameters:
% model: the model
% x: the mole fraction, it is an array of size (Nc,1)
% xi: the molar density of the mixture, unit: mol/m^3
% p: the pressure, unit: Pa
% phase: the phase of the fluid

% Return value:
% mu: the viscosity of the cell, unit: Pa*s 

% Reference: Chap 7. in <<Reservoir Simulation: Mathematical Techniques in Oil Recovery>>

% Author: Yuanqing Wu. Email: wuyuanq@gmail.com
% Last edited on November 23rd, 2012

function [ mu ] = RST_viscosity( model, x, xi, p, phase )

    % simplify the notations of model
    Nc = model.Nc;
    T = model.T;
    ct = model.ct;
    cp = model.cp;
    cv = model.cv;
    mw = model.mw;
       
    cp = cp*9.8692*1.D-6; % translate the unit from Pa to atm. 1 Pa = 9.8692*10^(-6) atm 
    mw = mw*1000; % translate the unit from kg/mol to g/mol  
       
    if(phase == 'l') % computation of the liquid viscosity at low pressure      
    
        CI = zeros(Nc,1);
         for i = 1 : Nc
            CI(i) = 1/(ct(i)^(1.D0/6.D0)/sqrt(mw(i))/cp(i)^0.666666666666D0);
         end
 
         RO = xi*1.d-3;
         
         XMIU = zeros(Nc,1);
         for i = 1 : Nc
            TR = T/ct(i);
            if (TR <= 1.5D0) 
                XMIU(i)=34.D-5*CI(i)*TR^.94D0;
            else
                XMIU(i)=17.78D-5*CI(i)*(4.58D0*TR-1.67D0)^0.625D0;
            end
         end
        
        S = 0.D0;
        S1 = 0.D0;
        GML = 0.D0;

        for i = 1 : Nc
            S = S + x(i)*XMIU(i)*sqrt(mw(i));
            S1 = S1 + x(i)*sqrt(mw(i));
            GML = GML + x(i)*mw(i);
        end    
        RO = RO*GML;  
        VS = S/S1;
        SV = 0.D0;

        isheavy = false;
        if(~isheavy)
            for i = 1 : Nc
                SV = SV + x(i)*cv(i)*mw(i);  
            end             
        else
            for i = 1 : 6
                SV = SV + x(i)*cv(i);
            end

            Z7PLUS = 0.D0;
            for i = 7 : Nc
                Z7PLUS = Z7PLUS + x(i);
            end
            MW7PLUS = 430.7; 
            GM7PLUS = 0.988; % AS GIVEN IN THE HEAVY OIL CASE
            V7PLUS = 21.573 + 0.015122*MW7PLUS - 27.656*GM7PLUS + 0.070615*MW7PLUS*GM7PLUS;
            V7PLUS = V7PLUS*0.3048^3/0.453592; % cub-ft/lb-mol -> cub-m/kg-mol

            SV = SV + Z7PLUS*V7PLUS;
        end
     
        ROC = 1.D0/SV;
        if(~isheavy) 
            ROR = RO/ROC/GML;
        else
            ROR = RO/ROC;
        end

        ST = 0.D0;
        SM = 0.D0;
        SPI = 0.D0;
                 
        for i = 1 : Nc
            SM = SM + x(i)*mw(i);
            ST = ST + x(i)*ct(i);
            SPI = SPI + x(i)*cp(i);
        end

        CE = ST^(1.D0/6.D0)/sqrt(SM)/SPI^(2.D0/3.D0);

        % AUGUST 16, 2007: If ROR is greater than 10, F is suspiciously large ...
        F = .1023D0+.023364D0*ROR+.058533D0*ROR^2- 0.040758D0*ROR^3  +.0093324D0*ROR^4; 

        COR = (F^4-1.D-4)/CE;
        mu = VS+COR;
    
    elseif(phase == 'g') % computation of the gas viscosity at low pressure 
        
        A = zeros(16,1);
        B = zeros(4,1);
        
        A(1) = -2.4621182;
        A(2) = 2.97054714;
        A(3) = -0.286264054;
        A(4) = 8.05420522*10^(-3);
        A(5) = 2.80860949;
        A(6) = -3.49803305;
        A(7) = 0.36037302;
        A(8) = -0.0104432413;
        A(9) = -0.793385684;
        A(10) = 1.39643306;
        A(11) = -0.149144925;
        A(12) = 4.41015512*10^(-3);
        A(13) = 0.0839387178;
        A(14) = -0.186408848;
        A(15) = 0.0203367881;
        A(16) = -6.09579263*10^(-4);

        B(1) = 4.0;
        B(2) = 4.0;
        B(3) = 4.0;
        B(4) = 4.0;
      
        PPC = 0.0;
        TPC = 0.0;
        GML = 0.0;

        for m = 1 : Nc
            PPC = PPC + x(m)*cp(m);
            TPC = TPC + x(m)*ct(m);
            GML = GML + x(m)*mw(m);
        end

        VS1 = (7.43+0.0133*GML)*(1.8*T)^1.5/(1.8*T+75.4+13.9*GML)*1.0*10^(-4);
        PPR = p*10^(-5)/PPC;
        TPR = T/TPC;
        for i = 1 : 4
            j = 4*i;
            B(i) = A(j-3)+(A(j-2)+(A(j-1)+A(j)*PPR)*PPR)*PPR;
        end

        CE = B(1)+(B(2)+(B(3)+B(4)*TPR)*TPR)*TPR;     
        VISR = exp(CE)/TPR;
        mu = VISR*VS1;

    end
    
    % translate the unit of viscosity from cp to Pa*s
    mu = mu * 0.001; % 1 cp = 0.001 Pa*s
    
end

