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
m = 0:1:gg.nI-1; n = 0:1:gg.nJ-1;  %switch to generic domain
[xx,yy] = meshgrid(m,n); 

%IS n defined correclty?


DP = gg.S_h * (gg.c_uh*gg.S_u'*(vv.u.*vv.lambda) + gg.c_vh*gg.S_v'*(vv.v.*vv.mu));
BB = idct2(vv.acoeff); BB = gg.S_h*BB(:);

agrad = zeros(size(vv.acoeff));

for j = 0:n_x-1;
for k = 0:n_y-1;

if k == 0; a1 = 1/sqrt(gg.nJ); else a1 = sqrt(2/gg.nJ); end
if j == 0; a2 = 1/sqrt(gg.nI); else a1 = sqrt(2/gg.nI); end

AA = a1*a2*cos(pi*(2*xx+1)*j/(2*gg.nI)).*cos(pi*(2*yy+1)*k/(2*gg.nJ));
AA = gg.S_h*AA(:);
agrad(j+1,k+1) = -2*pp.c8*sum(AA.*BB.*DP*gg.dx*gg.dy);
         
end
end



vv.agrad = agrad;

end
