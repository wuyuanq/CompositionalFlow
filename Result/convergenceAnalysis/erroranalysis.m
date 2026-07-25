
nx = 40;
ny = 40;
location = '/Users/wuy/Dropbox/Research/CompositionalFlow_results/convergenceAnalysis/';

ftxt{1} = [location,num2str(40),'/soln_cmp2PhFlw_Sw_raw.txt'];
ftxt{2} = [location,num2str(80),'/soln_cmp2PhFlw_Sw_raw.txt'];
ftxt{3} = [location,num2str(160),'/soln_cmp2PhFlw_Sw_raw.txt'];
ftxt{4} = [location,num2str(640),'/soln_cmp2PhFlw_Sw_raw.txt'];

% load fields
for n = 1 : 3
    field{n} = zeros(ny*2^(n-1), nx*2^(n-1));
    temp = load(ftxt{n});
    c = 0;
    for j = 1 : ny*2^(n-1)
        for i = 1 : nx*2^(n-1) 
            c = c + 1;
            field{n}(j,i) = temp(c);
        end
    end
end

field{4} = zeros(ny*2^4, nx*2^4);
temp = load(ftxt{4});
c = 0;
for j = 1 : ny*2^4
    for i = 1 : nx*2^4 
        c = c + 1;
        field{4}(j,i) = temp(c);
    end
end

norm = 0;
for j = 1 : ny*2^4    
    for i = 1 : nx*2^4 
        norm = norm + field{4}(j,i)^2;
    end
end
norm = sqrt(norm/(ny*2^4*nx*2^4));

% average and error
for n = 1 : 3
    convfield{n} = zeros(ny*2^(n-1), nx*2^(n-1));
    error{n} = zeros(ny*2^(n-1), nx*2^(n-1));
    aveerror = 0;
    for j = 1 : ny*2^(n-1)
        for i = 1 : nx*2^(n-1) 
            convfield{n}(j,i) = 0;
            for nj = 1 : 16/2^(n-1)
                for ni = 1 : 16/2^(n-1)
                    convfield{n}(j,i) = convfield{n}(j,i) + field{4}((j-1)*16/2^(n-1)+nj,(i-1)*16/2^(n-1)+ni);
                end
            end
            convfield{n}(j,i) = convfield{n}(j,i)/(16/2^(n-1))^2;
            error{n}(j,i) = abs(convfield{n}(j,i)-field{n}(j,i));
            aveerror = aveerror + error{n}(j,i)^2;
        end
    end
    aveerror = sqrt(aveerror/(nx*2^(n-1)*ny*2^(n-1)))/norm;
    disp(['The L2 norm error of grid ', num2str(n), ' is ', num2str(aveerror)]);
end


