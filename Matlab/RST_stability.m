% The function RST_stability is to dicide the phase status of the fluid.

% Input parameters:
% model: the model
% Pk: the pressure
% ZI: the composition of the fluid
% ztest: paramter
% Kstabtry: the parameter
% size_try: parameter

% Return value:
% judge: boolean variable. If the fluid is single phase, it is true; 
% if the fluid is two phases, it is false.
% K_stab: the initial guess of K (the equilibrium ratio)

% Author: Yuanqing Wu. Email: wuyuanq@gmail.com
% Last edited on January 23rd, 2014

function [ judge, K_stab ] = RST_stability( model, Pk, ZI, ztest, Kstabtry, size_try )

    Nc = model.Nc;

    XE = zeros(Nc, 1);
    YE = zeros(Nc, 1);
    XOLD = zeros(Nc, 1);
    XNEW = zeros(Nc, 1);
    GOLD = zeros(Nc, 1);
    GNEW = zeros(Nc, 1);
    FUGZ = zeros(Nc, 1);
    S = zeros(Nc, 1);
    YD = zeros(Nc, 1);
    iter_stab = zeros(Nc+8, 1);
    TPD_stab = zeros(Nc+8, 1);
    K_stab = zeros(Nc+8, Nc);

    judge = false;

    size_stab = 0;
    num_judge_iter = 0;

    y = zeros(Nc, 1);
    [ ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, phig ] = RST_fugacitycoef( model, y, ZI, Pk );

    FG = log(phig);
    for i_comp = 1 : Nc
        if(ztest(i_comp) <= 0.d0)
            continue
        end 
        FUGZ(i_comp) = FG(i_comp)+log(ztest(i_comp));
    end 

    for i_stab = 1 : size_try

        for i_comp = 1 : Nc
            XE(i_comp) = ztest(i_comp) / Kstabtry(i_stab,i_comp);
        end 

        X_CAL = XE/sum(XE);

        x = zeros(Nc, 1);
        [ ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, phil, ~ ] = RST_fugacitycoef( model, X_CAL, x, Pk );

        FL = log(phil);

        for i_comp = 1 : Nc
            if((ztest(i_comp) <= 0.d0)||(XE(i_comp) <= 0.d0))
                XOLD(i_comp) = 0.d0;
                GOLD(i_comp) = 0.d0;
                XNEW(i_comp) = 0.d0;
                continue
            end 
            XOLD(i_comp) = 2.D0*sqrt(XE(i_comp));
            TM = FUGZ(i_comp)-FL(i_comp);
            GOLD(i_comp) = XOLD(i_comp)/2.D0*(log(XE(i_comp))-TM);
            XNEW(i_comp) = 2.D0*exp(TM/2.D0);
        end 

        iteration_stab = 0;

        while(iteration_stab <= 100)

            XE(1:Nc) = 0.25D0*XNEW(1:Nc).^2;
            X_CAL = XE/sum(XE);

            x = zeros(Nc, 1);
            [ ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, phil, ~ ] = RST_fugacitycoef( model, X_CAL, x, Pk );

            FL = log(phil);
            SUM0=1.D0;
            for i_comp = 1 : Nc
                if((ztest(i_comp) <= 0.d0)||(XE(i_comp) <= 0.d0))
                    GNEW(i_comp) = 0.d0;
                    S(i_comp) = 0.d0;
                    YD(i_comp) = 0.d0;
                    continue
                end 

                TM = log(XE(i_comp))+FL(i_comp)-FUGZ(i_comp);
                SUM0 = SUM0+XE(i_comp)*(TM-1.D0);
                GNEW(i_comp) = XNEW(i_comp)/2.D0*TM;
                S(i_comp) = XNEW(i_comp)-XOLD(i_comp);
                YD(i_comp) = GNEW(i_comp)-GOLD(i_comp);
            end 

            [ YE ] = update(Nc, YD, S, GNEW, YE);

            for i_comp = 1 : Nc
                GOLD(i_comp) = GNEW(i_comp);
                XOLD(i_comp) = XNEW(i_comp);
                XNEW(i_comp) = XNEW(i_comp)-YE(i_comp);
            end 

            err_stab = 0.d0;
            for i_comp = 1 : Nc
                if(abs(XNEW(i_comp)-XOLD(i_comp))>=err_stab) 
                    err_stab=abs(XNEW(i_comp)-XOLD(i_comp));
                end 
            end 

            if(err_stab < 1.d-10) 
                break
            end 

            if(err_stab>1.d30)
                break
            end 

            iteration_stab = iteration_stab+1;

        end 

        if(err_stab >= 1.d-10) 
            num_judge_iter = num_judge_iter+1;
        end 

        if(err_stab < 1.d-10) % Convergence criteria
            num_judge_comp = 0;
            for i_comp = 1 : Nc
                if (abs(FL(i_comp)-FG(i_comp)) <= 1.d-5) 
                    num_judge_comp = num_judge_comp+1;
                end 
            end 
            if(num_judge_comp < Nc) 
                size_stab = size_stab+1;

                iter_stab(size_stab) = iteration_stab;
                TPD_stab(size_stab) = SUM0;
                for i_comp = 1 : Nc
                    K_stab(size_stab,i_comp) = exp(FL(i_comp)-FG(i_comp));
                end
            end
        end

    end

    if(size_stab == 0) 

        judge = true;
        if(num_judge_iter == size_try) 
            disp('NO CONVERGENCE IN STABILITY: ASSUMED SINGLE PHASE');
        end

    elseif(size_stab >= 2)

        matrix_stab = zeros(size_stab,Nc+2);
        matrix_stab(:,1:Nc) = K_stab(1:size_stab,1:Nc);
        matrix_stab(:,Nc+1) = iter_stab(1:size_stab);
        matrix_stab(:,Nc+2) = TPD_stab(1:size_stab);
        [ matrix_stab ] = piksrt(size_stab,Nc+2,matrix_stab);
        K_stab(1:size_stab,1:Nc) = matrix_stab(:,1:Nc);
        iter_stab(1:size_stab) = round(matrix_stab(:,Nc+1));
        TPD_stab(1:size_stab) = matrix_stab(:,Nc+2);

        m_stab = 1;
        while(m_stab < size_stab)
            num_judge_comp = 0;
            for i_comp = 1 : Nc
                if (abs(K_stab(m_stab,i_comp)/K_stab(m_stab+1,i_comp)-1.d0)<=1.d-5) 
                    num_judge_comp = num_judge_comp+1;
                end 
            end 
            if (num_judge_comp == Nc)
                size_stab = size_stab-1;
                TPD_stab(m_stab+1:size_stab) = TPD_stab(m_stab+2:size_stab+1);
                iter_stab(m_stab+1:size_stab) = iter_stab(m_stab+2:size_stab+1);
                K_stab(m_stab+1:size_stab,:) = K_stab(m_stab+2:size_stab+1,:);
            else
                m_stab = m_stab+1;
            end 
        end 

        TPD_stab(size_stab+1:Nc+8) = 0.d0;
        K_stab(size_stab+1:Nc+8,:) = 0.d0;

    end

    if (TPD_stab(1) >= -1.d-10) 
        judge = true;
    end
    
end

function [ X1 ] = update(N, Y1, SS, GD, X1)

    SY = 0D0;
    HYY = 0;
    for i = 1 : N
        X1(i) = Y1(i);
        XI = X1(i);
        SY = SY+SS(i)*Y1(i);
        HYY = HYY+XI*Y1(i);
    end 

    FAC = (1D0+HYY/SY)*0.5D0;
    X1 = (X1-FAC*SS)/SY;

    XG = 0D0;
    SG0 = 0D0;
    for i = 1 : N
        XG = XG+X1(i)*GD(i);
        SG0 = SG0+SS(i)*GD(i);
    end 

    X1 = GD-XG*SS-SG0*X1;
        
end
      
% Sort 'matrix' into ascending numerical order of the last column by straight insertion method.
function [ matrix ] = piksrt(row, column, matrix)

    for j = 2 : row
        array = matrix(j,:);
        jump = 0;
        for i = j-1 : -1 : 1
            if(matrix(i,column) <= array(column))
                matrix(i+1,:) = array;
                jump = 1;
                break
            end 
            matrix(i+1,:) = matrix(i,:);
        end 
        if(jump == 0)
            i = 0;
            matrix(i+1,:) = array;
        end
    end 
    
end
        
        
        
        