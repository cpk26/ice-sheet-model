function [vv] = ism_cslip_dct(vv,pp,gg,oo )
%% Two-dimensional discrete fourier cosine transform
% Improve to fast fourier cosine transform?
% Inputs:
%   vv      struct containing initial solution variables
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv2     struct containing new solution variables

%Sw

n_x = vv.n_x;               %Number of coefficients (x/y dirs)
n_y = vv.n_y;   
Lx = gg.Lx;                 %Length of Domain (x/y dimensions)
Ly = gg.Ly;
FF = sqrt(vv.C(:));         %Function to transform

acoeff = zeros(n_x+1,n_y+1);

%% Calculate Coefficients of DCT

for j=0:n_x 
    for k = 0:n_y
        
        if j == 0 && k == 0
            acoeff(j+1,k+1) = (1/(Lx*Ly)) * sum(FF)*gg.dx*gg.dy;
            
        elseif k == 0
            AA = 2/(Lx*Ly);
            BB = sum(FF.*cos(pi*j*gg.xx(:)/Lx))*gg.dx*gg.dy;
            acoeff(j+1,k+1) = AA*BB;
            
        elseif j == 0
            AA = 2/(Lx*Ly);
            BB = sum(FF.*cos(pi*k*gg.yy(:)/Ly))*gg.dx*gg.dy;
            acoeff(j+1,k+1) = AA*BB;
            
        else
            AA = 4/(Lx*Ly);
            BB = sum(FF.*cos(pi*j*gg.xx(:)/Lx).*cos(pi*k*gg.yy(:)/Ly))*gg.dx*gg.dy;
            acoeff(j+1,k+1) = AA*BB;
            
        end      
    end
end
    
vv.acoeff = acoeff(:);
end
