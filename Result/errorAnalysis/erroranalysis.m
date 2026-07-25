
nx = 160;
ny = 160;
location = '/Users/wuy/dropbox/research/CompositionalFlow_results/errorAnalysis/';
fporo1txt = [location,'soln_cmp2PhFlw_x1_raw_160.txt'];
fporo2txt = [location,'soln_cmp2PhFlw_x1_raw_320.txt'];
fporo3txt = [location,'soln_cmp2PhFlw_x1_raw_640.txt'];

% load porosities
poro1 = zeros(ny, nx);
temp = load(fporo1txt);
k = 1;
for j = 1 : ny 
    for i = 1 : nx  
        poro1(j,i) = temp(k);
        k = k + 1;
    end
end

poro2 = zeros(ny*2, nx*2);
temp = load(fporo2txt);
k = 1;
for j = 1 : ny*2 
    for i = 1 : nx*2  
        poro2(j,i) = temp(k);
        k = k + 1;
    end
end

poro3 = zeros(ny*4, nx*4);
temp = load(fporo3txt);
k = 1;
for j = 1 : ny*4
    for i = 1 : nx*4   
        poro3(j,i) = temp(k);
        k = k + 1;
    end
end

% average
poro21 = zeros(ny, nx);
for j = 1 : ny
    for i = 1 : nx   
        poro21(j,i) = (poro2((j-1)*2+1,(i-1)*2+1)+poro2((j-1)*2+1,(i-1)*2+2)+ ...
            poro2((j-1)*2+2,(i-1)*2+1)+poro2((j-1)*2+2,(i-1)*2+2))/4;
    end
end

poro31 = zeros(ny, nx);
for j = 1 : ny
    for i = 1 : nx   
        poro31(j,i) = (poro3((j-1)*4+1,(i-1)*4+1)+poro3((j-1)*4+1,(i-1)*4+2)+ ...
            poro3((j-1)*4+1,(i-1)*4+3)+poro3((j-1)*4+1,(i-1)*4+4)+ ...
            poro3((j-1)*4+2,(i-1)*4+1)+poro3((j-1)*4+2,(i-1)*4+2)+ ...
            poro3((j-1)*4+2,(i-1)*4+3)+poro3((j-1)*4+2,(i-1)*4+4)+ ...
            poro3((j-1)*4+3,(i-1)*4+1)+poro3((j-1)*4+3,(i-1)*4+2)+ ...
            poro3((j-1)*4+3,(i-1)*4+3)+poro3((j-1)*4+3,(i-1)*4+4)+ ...
            poro3((j-1)*4+4,(i-1)*4+1)+poro3((j-1)*4+4,(i-1)*4+2)+ ...
            poro3((j-1)*4+4,(i-1)*4+3)+poro3((j-1)*4+4,(i-1)*4+4))/16;
    end
end

% calculate errors
error23 = zeros(ny, nx);
for j = 1 : ny
    for i = 1 : nx 
        error23(j,i) = abs(poro21(j,i)-poro31(j,i));
    end
end

aveerr23 = 0;
for j = 1 : ny
    for i = 1 : nx 
        aveerr23 = aveerr23 + error23(j,i)^2;
    end
end
aveerr23 = sqrt(aveerr23/(nx*ny));

error13 = zeros(ny, nx);
for j = 1 : ny
    for i = 1 : nx 
        error13(j,i) = abs(poro31(j,i)-poro1(j,i));
    end
end

aveerr13 = 0;
for j = 1 : ny
    for i = 1 : nx 
        aveerr13 = aveerr13 + error13(j,i)^2;
    end
end
aveerr13 = sqrt(aveerr13/(nx*ny));

p = log2(aveerr13/aveerr23-1);

p
aveerr13
aveerr23






