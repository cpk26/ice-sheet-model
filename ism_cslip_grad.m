function [vv] = ism_cslip_grad(vv, pp, gg, oo)
%% Calculate gradients of alpha coefficients
% Inputs:
%   n_x,n_y   number of alpha coefficients in the x,y directions
%   vv      solution variables
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv       solution variables

n_x = vv.n_x;           %Number of coefficients (x/y dirs)
n_y = vv.n_y;
Lx = gg.Lx;             %Length of Domain (x/y dimensions)
Ly = gg.Ly;

DP = gg.S_h * (gg.c_uh*gg.S_u'*(vv.u.*vv.lambda) + gg.c_vh*gg.S_v'*(vv.v.*vv.mu));
BB = gg.S_h*sqrt(vv.C(:));

agrad = zeros((n_x+1)*(n_y+1),1);

for j=0:n_x;
for k = 0:n_y
        
aInd = (j+1) + k*n_y;
agrad(aInd) = -2*pp.c8*sum(cos(pi*j*gg.xx(:)/Lx).*cos(pi*k*gg.yy(:)/Ly).*BB.*DP*gg.dx*gg.dy);
        
       
end
end



vv.agrad = agrad;

end

