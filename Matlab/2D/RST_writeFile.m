% This is the RST_writeFile() function which outputs the 
% results of RST_compositionalTwoPhaseFlow() to a series of solution files.

% Input parameters:
% model: the model
% pw: the pressure
% ux: the velocities in the x dirction
% uy: the velocities in the y dirction
% x: the mole fraction
% Sw: the saturation of wetting phase
% xi: the molar density

% Return value:
% fpwtxt: the pressure file of wetting phase in the txt form
% fuxtxt: the x-dirction velocity file in the txt form
% fuytxt: the y-dirction velocity file in the txt form
% fmftxt: the mole fraction files in the txt form
% fswtxt: the saturation of wetting phase file in the txt form
% fxitxt: the molar density in the txt form

% Author: Yuanqing Wu. Email: wuyuanq@gmail.com
% Last edited on October 23rd, 2012

function [ fpwtxt, fuxtxt, fuytxt, fmftxt, fswtxt, fxitxt ] = RST_writeFile( model, pw, ux, uy, x, sw, xi )

    % simplify the symbols of the model
    Nc = model.Nc;
    nx = model.nx;
    ny = model.ny;
    nt = model.nt;
    xs = model.xs;
    ys = model.ys;
    ts = model.ts;
    soludoc = model.soludoc;
    
    % create the new document to store the results
    if(exist(soludoc, 'dir') == 0)
        mkdir(soludoc);
    end
    
    % define the names of the solution files
    fpw = [soludoc, '/soln_cmp2PhFlw_Pw_RSTo.m'];
    fux = [soludoc, '/soln_cmp2PhFlw_Ux_RSTo.m'];
    fuy = [soludoc, '/soln_cmp2PhFlw_Uy_RSTo.m'];
    fsw = [soludoc, '/soln_cmp2PhFlw_Sw_RSTo.m'];
    fxi = [soludoc, '/soln_cmp2PhFlw_c_RSTo.m'];
    fpwtxt = [soludoc, '/soln_cmp2PhFlw_Pw_raw.txt'];
    fuxtxt = [soludoc, '/soln_cmp2PhFlw_Ux_raw.txt'];
    fuytxt = [soludoc, '/soln_cmp2PhFlw_Uy_raw.txt'];
    fswtxt = [soludoc, '/soln_cmp2PhFlw_Sw_raw.txt'];
    fxitxt = [soludoc, '/soln_cmp2PhFlw_c_raw.txt'];
    fmf = [];
    for m = 1 : Nc
        fk = [soludoc, '/soln_cmp2PhFlw_x', num2str(m), '_RSTo.m'];
        fmf = [fmf; fk];
    end
    cdfmf = cellstr(fmf);
    fmftxt = [];
    for m = 1 : Nc
        fk = [soludoc, '/soln_cmp2PhFlw_x', num2str(m), '_raw.txt'];
        fmftxt = [fmftxt; fk];
    end
    cdfmftxt = cellstr(fmftxt);
    
    fid = fopen(fpwtxt, 'w');
    for j = 2 : nx+1
        for i = 2 : ny+1    
            fprintf(fid, '%10e\n', pw(i,j));
        end
    end
    fclose(fid);
    
    fid = fopen(fuxtxt, 'w');
    for j = 1 : nx+1
        for i = 1 : ny
            fprintf(fid, '%10e\n', ux(i,j));
        end
    end
    fclose(fid);
    
    fid = fopen(fuytxt, 'w');
    for j = 1 : nx
        for i = 1 : ny+1  
            fprintf(fid, '%10e\n', uy(i,j));
        end
    end
    fclose(fid); 
    
    for m = 1 : Nc
        fid = fopen(char(cdfmftxt(m)), 'w');
        for j = 1 : nx
            for i = 1 : ny  
                fprintf(fid, '%10e\n', x(m,i,j));
            end
        end
        fclose(fid); 
    end
    
    fid = fopen(fswtxt, 'w');
    for j = 1 : nx
        for i = 1 : ny   
            fprintf(fid, '%10e\n', sw(i,j));
        end
    end
    fclose(fid);
    
    fid = fopen(fxitxt, 'w');
    for j = 1 : nx
        for i = 1 : ny  
            fprintf(fid, '%10e\n', xi(i,j));
        end
    end
    fclose(fid);
    
    fid = fopen(fpw, 'w');
    fprintf(fid, 'Pw.type = ''cell-centered_data'';\n');
    fprintf(fid, 'Pw.simTime =      %.4f    ;\n', ts(nt+1));
    fprintf(fid, 'Pw.mesh.type = ''rectangular_mesh'';\n');
    fprintf(fid, 'Pw.mesh.xs = [');
    for i = 1 : nx
        fprintf(fid, '%.4f,\t', xs(i));
    end
    i = i + 1;
    fprintf(fid, '%.4f\t', xs(i));
    fprintf(fid, '];\n');
    fprintf(fid, 'Pw.mesh.ys = [');
    for i = 1 : ny
        fprintf(fid, '%.4f,\t', ys(i));
    end
    i = i + 1;
    fprintf(fid, '%.4f\t', ys(i));
    fprintf(fid, '];\n');
    fprintf(fid, 'Pw.array = zeros(%d, %d);\n', nx, ny); 
    fprintf(fid, 'Pw.array(1:%d, 1:%d) = reshape( load(''%s''), [%d, %d] );', nx, ny, fpwtxt, nx, ny); 
    fclose(fid);

    fid = fopen(fux, 'w');
    fprintf(fid, 'Ux.type = ''face-centered_data'';\n');
    fprintf(fid, 'Ux.simTime =      %.4f    ;\n', ts(nt+1));
    fprintf(fid, 'Ux.mesh.type = ''rectangular_mesh'';\n');
    fprintf(fid, 'Ux.mesh.xs = [');
    for i = 1 : nx
        fprintf(fid,'%.4f,\t', xs(i));
    end
    i = i + 1;
    fprintf(fid, '%.4f\t', xs(i));
    fprintf(fid, '];\n');
    fprintf(fid, 'Ux.mesh.ys = [');
    for i = 1 : ny
        fprintf(fid, '%.4f,\t', ys(i));
    end
    i = i + 1;
    fprintf(fid, '%.4f\t', ys(i));
    fprintf(fid, '];\n');
    fprintf(fid, 'Ux.array = zeros(%d,%d);\n', nx+1, ny);
    fprintf(fid, 'Ux.array(1:%d, 1:%d) = reshape( load(''%s''), [%d, %d] );\n', nx+1, ny, fuxtxt, nx+1, ny); 
    fclose(fid);
    
    fid = fopen(fuy, 'w');
    fprintf(fid, 'Uy.type = ''face-centered_data'';\n');
    fprintf(fid, 'Uy.simTime =      %.4f    ;\n', ts(nt+1));
    fprintf(fid, 'Uy.mesh.type = ''rectangular_mesh'';\n');
    fprintf(fid, 'Uy.mesh.xs = [');
    for i = 1 : nx
        fprintf(fid, '%.4f,\t', xs(i));
    end
    i = i + 1;
    fprintf(fid, '%.4f\t', xs(i));
    fprintf(fid, '];\n');
    fprintf(fid, 'Uy.mesh.ys = [');
    for i = 1 : ny
        fprintf(fid, '%.4f,\t', ys(i));
    end
    i = i + 1;
    fprintf(fid, '%.4f\t', ys(i));
    fprintf(fid, '];\n');
    fprintf(fid, 'Uy.array = zeros(%d,%d);\n', nx, ny+1);
    fprintf(fid, 'Uy.array(1:%d, 1:%d) = reshape( load(''%s''), [%d, %d] );\n', nx, ny+1, fuytxt, nx, ny+1); 
    fclose(fid);

    for m = 1 : Nc
        fid = fopen(char(cdfmf(m)), 'w');
        fprintf(fid, 'x%d.type = ''cell-centered_data'';\n',m);
        fprintf(fid, 'x%d.simTime =      %.4f    ;\n', m, ts(nt+1));
        fprintf(fid, 'x%d.mesh.type = ''rectangular_mesh'';\n', m);
        fprintf(fid, 'x%d.mesh.xs = [', m);
        for i = 1 : nx
            fprintf(fid, '%.4f,\t', xs(i));
        end
        i = i + 1;
        fprintf(fid, '%.4f\t', xs(i));
        fprintf(fid, '];\n');
        fprintf(fid, 'x%d.mesh.ys = [', m);
        for i = 1 : ny
            fprintf(fid, '%.4f,\t', ys(i));
        end
        i = i + 1;
        fprintf(fid, '%.4f\t', ys(i));
        fprintf(fid, '];\n');
        fprintf(fid, 'x%d.array = zeros(%d,%d);\n', m, nx, ny);
        fprintf(fid, 'x%d.array(1:%d, 1:%d) = reshape( load(''%s''), [%d, %d] );\n', m, nx, ny, char(cdfmftxt(m)), nx, ny); 
        fclose(fid);
    end
    
    fid = fopen(fsw, 'w');
    fprintf(fid, 'Sw.type = ''cell-centered_data'';\n');
    fprintf(fid, 'Sw.simTime =      %.4f    ;\n', ts(nt+1));
    fprintf(fid, 'Sw.mesh.type = ''rectangular_mesh'';\n');
    fprintf(fid, 'Sw.mesh.xs = [');
    for i = 1 : nx
        fprintf(fid, '%.4f,\t', xs(i));
    end
    i = i + 1;
    fprintf(fid, '%.4f\t', xs(i));
    fprintf(fid, '];\n');
    fprintf(fid, 'Sw.mesh.ys = [');
    for i = 1 : ny
        fprintf(fid, '%.4f,\t', ys(i));
    end
    i = i + 1;
    fprintf(fid, '%.4f\t', ys(i));
    fprintf(fid, '];\n');
    fprintf(fid, 'Sw.array = zeros(%d, %d);\n', nx, ny); 
    fprintf(fid, 'Sw.array(1:%d, 1:%d) = reshape( load(''%s''), [%d, %d] );', nx, ny, fswtxt, nx, ny); 
    fclose(fid);
    
    fid = fopen(fxi, 'w');
    fprintf(fid, 'c.type = ''cell-centered_data'';\n');
    fprintf(fid, 'c.simTime =      %.4f    ;\n', ts(nt+1));
    fprintf(fid, 'c.mesh.type = ''rectangular_mesh'';\n');
    fprintf(fid, 'c.mesh.xs = [');
    for i = 1 : nx
        fprintf(fid, '%.4f,\t', xs(i));
    end
    i = i + 1;
    fprintf(fid, '%.4f\t', xs(i));
    fprintf(fid, '];\n');
    fprintf(fid, 'c.mesh.ys = [');
    for i = 1 : ny
        fprintf(fid, '%.4f,\t', ys(i));
    end
    i = i + 1;
    fprintf(fid, '%.4f\t', ys(i));
    fprintf(fid, '];\n');
    fprintf(fid, 'c.array = zeros(%d, %d);\n', nx, ny); 
    fprintf(fid, 'c.array(1:%d, 1:%d) = reshape( load(''%s''), [%d, %d] );', nx, ny, fxitxt, nx, ny); 
    fclose(fid);
    
end

