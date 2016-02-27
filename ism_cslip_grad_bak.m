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


DP = (gg.c_uh*(vv.u.*vv.lambda) + gg.c_vh*(vv.v.*vv.mu));
BB = exp(idct2(vv.acoeff)); BB = gg.S_h*BB(:);

agrad = zeros(size(vv.acoeff));


%aHess = zeros(gg.nIJ,gg.nIJ);
%[kInd,jInd] = meshgrid([0:n_x-1],[0:n_y-1]);
%parfor ii = 1:(n_x)*(n_y)

for j = 0:n_y-1;
 for k = 0:n_x-1;
   
%j = jInd(ii);
%k = kInd(ii);

if j == 0; a1 = 1/sqrt(gg.nJ); else a1 = sqrt(2/gg.nJ); end
if k == 0; a2 = 1/sqrt(gg.nI); else a1 = sqrt(2/gg.nI); end

AA = a1*a2*cos(pi*(2*yy+1)*j/(2*gg.nJ)).*cos(pi*(2*xx+1)*k/(2*gg.nI));
AA = gg.S_h*AA(:);
agrad(j+1,k+1) = -pp.c8*sum(AA.*BB.*DP*gg.dx*gg.dy);

% parfor m=1:n_y                                 %Calculate Hessian
%     for n=1:n_x
%         e = zeros(size(vv.acoeff));
%         e(m,n) = 1e-5;
%         BB2 = exp(idct2(vv.acoeff + e)); BB2 = gg.S_h*BB2(:);
%         
%         ind = (j+1) + (k)*n_y;
%         aHess(ind,:) = (pp.c8*sum(AA.*BB2.*DP*gg.dx*gg.dy)-agrad(j+1,k+1))./1e-5;
%     end
% end
         
end
end

vv.agrad = agrad;

end

