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
x = linspace(0,1,gg.nI); x(1) = []; x(end) = []; %switch to generic domain
y = linspace(1,0,gg.nJ); y(1) = []; y(end) = [];
[xx,yy] = meshgrid(x,y); 
dx = abs(xx(1,2) - xx(1,1));
dy = abs(yy(1,1) - yy(2,1));
xx = xx(:); yy = yy(:);
FF = sum(sqrt(vv.C(:)));         %Function to transform

acoeff = zeros(n_x+1,n_y+1);

%% Calculate Coefficients of DCT

for j=0:n_x 
    for k = 0:n_y
        
        if j == 0 && k == 0
            acoeff(j+1,k+1) = sum(FF)*dx*dy;
            
        elseif k == 0
            acoeff(j+1,k+1) = 2*sum(FF.*cos(pi*j*xx))*dx*dy;
            
        elseif j == 0
            acoeff(j+1,k+1) = 2*sum(FF.*cos(pi*k*yy))*dx*dy;
            
        else
            acoeff(j+1,k+1) = 4*sum(FF.*cos(pi*j*xx).*cos(pi*k*yy))*dx*dy;
            
        end      
    end
end
    
vv.acoeff = acoeff(:);
end
