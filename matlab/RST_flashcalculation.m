
% The function RST_flashcalculation is to use the PR-EOS and successive
% substitution technique to do flash calculation.

% Input parameters:
% model: the model
% P: the pressure in the cell
% z: the mole fraction of each component in the cell

% Return values:
% x: the mole fraction in liquid phase of each component
% y: the mole fraction in gas phase of each component
% xiL: the molar density in liquid phase, unit: mol/m^3
% xiG: the molar density in gas phase, unit: mol/m^3
% rhoL: the mass density of liguid phase, unit: kg/m^3
% rhoG: the mass density of gas phase, unit: kg/m^3
% sL: the saturation of wetting phase
% v: the partial molar volume, unit: m^3/mol
% Cf: the total mixture compressibility
% isW: is there wetting phase?
% isN: is there nonwetting phase?

% Reference: Chap 7. in <<Reservoir Simulation: Mathematical Techniques in Oil Recovery>>
%            Chap 4. in <<Thermodynamics of Hydrocarbon Reservoirs>>

% Author: Yuanqing Wu. Email: wuyuanq@gmail.com
% Last edited on January 5th, 2014      

function [ x, y, xiL, xiG, rhoL, rhoG, sL, v, Cf, isW, isN ] = RST_flashcalculation( model, P, z )
    
    % simplify the notations of model
    Nc = model.Nc;
    T = model.T;
    Tc = model.ct;
    Pc = model.cp;
    omega = model.af;
    psatA = model.psatA;
    psatB = model.psatB;
    psatC = model.psatC;
    delta = model.delta;
    
    R = 8.314;
    MAXTIME = 1000;
    
    x = zeros(Nc, 1); % mole fraction in liquid phase of each component
    y = zeros(Nc, 1); % mole fraction in gas phase of each component
    
    v = zeros(Nc, 1); % the partial molar volume
    
    % decide whether it is pure substance
    ispure = false;
    purenum = 0;
    for m = 1 : Nc
        if(abs(z(m)-1.D0) < 1.D-12)
            ispure = true;
            purenum = m;
            break;
        end
    end
    
    % if it is pure substance
    if(ispure)
        
        [ ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, xiL, xiG, rhoL, rhoG, CfL, CfG ] ...
            = RST_PREOS( model, z, z, P );
    
        if(T > Tc(purenum))
            x = 0;
            y = z;
            xiL = 0;
            rhoL = 0;
            Cf = CfG;
            v(purenum) = 1/xiG;
            sL = 0.0;
            isW = false;
            isN = true;
        else
            psatmix = exp(psatA(purenum)-psatB(purenum)/((T-273)+psatC(purenum)));
            psatmix = psatmix*133.322;
            if(P > psatmix)
                x = z;
                y = 0;
                xiG = 0;
                rhoG = 0;
                Cf = CfL;
                v(purenum) = 1/xiL;
                sL = 1.0;
                isN = false;
                isW = true;
            else
                x = 0;
                y = z;
                xiL = 0;
                rhoL = 0;
                Cf = CfG;
                v(purenum) = 1/xiG;
                sL = 0.0;
                isW = false;
                isN = true;
            end
        end
        return
    end
    
    % step 1: guess the initial values of K(i) at the fixed temperature and pressure 
    K = zeros(Nc, 1);
    for i = 1 : Nc
        K(i) = exp(5.37*(1+omega(i))*(1-Tc(i)/T)+log(Pc(i)/P));
    end
    
    ztest1f = z;

    Kstab1ftry = zeros(Nc+8, Nc);

    Kstab1ftry(1,:) = K;
    Kstab1ftry(2,:) = 1.d0./K;
    Kstab1ftry(3,:) = K.^(1.d0/3.d0);
    Kstab1ftry(4,:) = 1.d0./K.^(1.d0/3.d0);
    for i = 1 : Nc
        for j = 1 : Nc
            if (i == j)
                Kstab1ftry(i+4,j) = 0.9d0/ztest1f(j);
            else
                Kstab1ftry(i+4,j) = 0.1d0/(Nc-1)/ztest1f(j);
            end
        end 
    end 
    Kstab1ftry(Nc+5:Nc+8,:) = 0.d0;
    size1f_try = Nc+4;

    [ SINGLE, K_stab_1F ] = RST_stability( model, P, z, ztest1f, Kstab1ftry, size1f_try );
    
    if(SINGLE)

        [ ZL, ZG, am, bm, al, ag, bl, bg, ~, ~, ~, ~, xiL, xiG, rhoL, rhoG, CfL, CfG] ...
            = RST_PREOS( model, z, z, P );

        tempsum = zeros(Nc, 1);

        for i = 1 : Nc
            for j = 1 : Nc
                tempsum(i) = tempsum(i) + z(j)*(1-delta(i,j))*sqrt(am(i)*am(j));
            end 
            tempsum(i) = 2*tempsum(i);
        end 

        mt = 0;
        for m = 1 : Nc
            mt = mt + Tc(m)*z(m);
        end 

        if(T > mt)
            x = 0;
            y = z;
            xiL = 0;
            rhoL = 0;
            sL = 0;
            Cf = CfG;

            molevg = ZG*R*T/P;
            for m = 1 : Nc
                v(m) = (R*T/(molevg-bg)*(1+bm(m)/(molevg-bg))-(tempsum(m)-(2*ag*bm(m)*(molevg-bg)) ...
                    /(molevg*(molevg+bg)+bg*(molevg-bg)))/(molevg*(molevg+bg)+bg*(molevg-bg)))/(R*T ...
                    /(molevg-bg)^2-2*ag*(molevg+bg)/(molevg*(molevg+bg)+bg*(molevg-bg))^2);
                if(abs(z(m)) < 1.D-12) 
                    v(m) = 0.0;
                end 
            end 

            isW = false;
            isN = true;
        else
            %compute the saturation vapor pressure of each component at temperature T, Antoine equation
            psat = zeros(Nc, 1);
            for m = 1 : Nc
                psat(m) = exp(psatA(m)-psatB(m)/((T-273)+psatC(m)));
                psat(m) = psat(m)*133.322;
            end 
            psatmix = 0;
            for m = 1 : Nc
                psatmix = psatmix + psat(m)*z(m);
            end 
            if(P > psatmix)
                x = z;
                y = 0;
                xiG = 0;
                rhoG = 0;
                sL = 1;
                Cf = CfL;

                molevl = ZL*R*T/P;
                for m = 1 : Nc
                    v(m) = (R*T/(molevl-bl)*(1+bm(m)/(molevl-bl))-(tempsum(m)-(2*al*bm(m)* ...
                        (molevl-bl))/(molevl*(molevl+bl)+bl*(molevl-bl)))/(molevl*(molevl+bl)+bl* ...
                        (molevl-bl)))/(R*T/(molevl-bl)^2-2*al*(molevl+bl)/(molevl*(molevl+bl)+ ...
                        bl*(molevl-bl))^2);
                    if(abs(z(m)) < 1.D-12) 
                        v(m) = 0.0;
                    end 
                end 

                isN = false;
                isW = true;
            else
                x = 0;
                y = z;
                xiL = 0;
                rhoL = 0;
                sL = 0;
                Cf = CfG;

                molevg = ZG*R*T/P;
                for m = 1 : Nc
                    v(m) = (R*T/(molevg-bg)*(1+bm(m)/(molevg-bg))-(tempsum(m)-(2*ag*bm(m)* ...
                        (molevg-bg))/(molevg*(molevg+bg)+bg*(molevg-bg)))/(molevg*(molevg+bg)+bg* ...
                        (molevg-bg)))/(R*T/(molevg-bg)^2-2*ag*(molevg+bg)/(molevg*(molevg+bg)+ ...
                        bg*(molevg-bg))^2);
                    if(abs(z(m)) < 1.D-12) 
                        v(m) = 0.0;
                    end 
                end 

                isW = false;
                isN = true;
            end 
        end 

    else

        K = K_stab_1F(1,:);

        for t = 1 : MAXTIME

            h0 = 0;
            for i = 1 : Nc
                h0 = h0 + K(i)*z(i);
            end 
            h0 = h0 - 1;
            h1 = 0;
            for i = 1 : Nc
                h1 = h1 + z(i)/K(i);
            end 
            h1 = 1 - h1;

            if(h0 <= 0) 
                x = z;
                yt = 0;
                for i = 1 : Nc
                    yt = yt + K(i)*z(i);
                end 
                y = K.*z/yt;
                [ K, xiL, xiG, rhoL, rhoG, sL, v, Cf, isW, isN, jump ] = submodule( model, x, y, z, P );
                if(jump)
                    break
                end
                continue
            end 
                
            if(h1 >= 0)
                y = z;
                xt = 0;
                for i = 1 : Nc
                    xt = xt + z(i)/K(i);
                end 
                x = z./K/xt;
                [ K, xiL, xiG, rhoL, rhoG, sL, v, Cf, isW, isN, jump ] = submodule( model, x, y, z, P );
                if(jump)
                    break
                end
                continue
            end 

            beta = 0.5;
            lp = 0.0;
            rp = 1.0;
            while(true)
                hbeta = 0;
                for i = 1 : Nc
                    hbeta = hbeta + ((K(i)-1)*z(i))/(1+beta*(K(i)-1));
                end 
                if(abs(hbeta) < 1.D-12)
                    break
                end 

                if(hbeta > 0) 
                    lp = beta;
                else
                    rp = beta;
                end 
                beta = (lp + rp)/2.0;
            end 

            for i = 1 : Nc
                x(i) = z(i)/(1+(K(i)-1)*beta);
                y(i) = K(i)*x(i);
            end 

            [ K, xiL, xiG, rhoL, rhoG, sL, v, Cf, isW, isN, jump ] = submodule( model, x, y, z, P );
            if(jump)
                break
            end
        end 
    end 
end

function [ K, xiL, xiG, rhoL, rhoG, sL, v, Cf, isW, isN, jump ] = submodule( model, x, y, z, P )

    Nc = model.Nc;
    jump =false;
    sL = 0;
    v = zeros(Nc, 1);
    Cf = 0;
    isW = false;
    isN = false;
    
    [ criteria, K, ZL, ZG, am, bm, al, ag, bl, bg, AL, AG, BL, BG, xiL, xiG, rhoL, rhoG, CfL, CfG ] ...
        = substeps( model, x, y, P );

    if(rhoL <= rhoG) 
        [ criteria, K, ZL, ZG, am, bm, al, ag, bl, bg, AL, AG, BL, BG, xiL, xiG, rhoL, rhoG, CfL, CfG ] ...
            = substeps( model, y, x, P );
    end

    if(criteria < 1.D-12)

        for i = 1 : Nc
            if(z(i) ~= 0)
                nn1 = i;
                break
            end 
        end 
        for i = nn1+1 : Nc
            if(z(i) ~= 0)
                nn2 = i;
                break
            end 
        end 

        sL = (y(nn2)*z(nn1)-y(nn1)*z(nn2))/(y(nn2)*z(nn1)-y(nn1)*z(nn2)+(x(nn1)*z(nn2)-x(nn2)*z(nn1))*xiL/xiG);

        phasemole(1) = (y(nn2)*z(nn1)-y(nn1)*z(nn2))/(x(nn1)*y(nn2)-x(nn2)*y(nn1))*(sL*xiL+(1-sL)*xiG);
        phasemole(2) = (x(nn1)*z(nn2)-x(nn2)*z(nn1))/(x(nn1)*y(nn2)-x(nn2)*y(nn1))*(sL*xiL+(1-sL)*xiG);

        for i = 1 : Nc
            v(i) = RST_pmv( model, i, P, z, 1/(phasemole(1)+phasemole(2)), K );
        end 

        [ Cf ] = RST_comprefac( model, x, y, P, ZL, ZG, am, bm, al, ag, bl, bg, ...
            AL, AG, BL, BG, xiL, xiG, sL );

        % When sL is very small, the algorithm to compute Cf in 2 phases will return NaN,
        % so we have to process such condition.
        if(isnan(Cf)) 
            if(sL < 1.D-7) 
                Cf = CfG;
            else
                Cf = CfL;
            end 
        end 

        isW = true;
        isN = true;

        jump = true;
    end

end



