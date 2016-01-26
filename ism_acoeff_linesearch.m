
function [vv2] = ism_acoeff_linesearch(vv,aa, pp, gg, oo)
%% Newton-Raphson step to update inversion alpha coefficients 
% Based on code from C.T. Kelley, from SIAM's 'Solving Nonlinear Equations
% with Newton's Method'
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv2     updated struct with new alpha coefficients

disp('Newton-Raphson Step')

if ~isfield(vv,'tau'), vv.tau = -1; end;        %Initial step coefficient

%% Paramaters
vv2 = vv;               

maxarm = 25;            %Maximum number of steps (reductions)
iarm = 0;               %Counter of number of reductions
sigma = .5;            %Initial step size reduction
alpha = 1.d-4;          %Line search parameter
armflag = 0;            %Flag for exceeding max number of reductions
n_x = vv.n_x;           %Number of coefficients for basis in x/y dirs
n_y = vv.n_y;

[mft0] = ism_inversion_misfit(vv.u,vv.v,aa,pp,gg, oo);  %Initial misfit

tau = vv.tau; taum = 1; tauc = tau;     %Step coefficient. 'm': previous step; 'c': current step
step = mft0/norm(vv.agrad(:),1) * vv.agrad / norm(vv.agrad(:),1);

%% Initial step
vv2.acoeff = vv2.acoeff - tau*step;
vv2.C = idct2(vv2.acoeff);       %Calculate new slipperiness field

[ii] = ism_sia(aa.s,aa.h,vv2.C, pp,gg,oo);  %Calculate corresponding velocities 
vv2.u = ii.u; vv2.v = ii.v;                 %SIA                                        
[vv2] = ism_sstream(vv2,aa,pp,gg,oo );      %SSA 

[mft] = ism_inversion_misfit(vv2.u,vv2.v,aa,pp,gg, oo); %Current misfit
mft0_2 = mft0*mft0; mftc_2 = mft*mft; mftm_2 = mft*mft; %Misfit squared; 'm': previous step; 'c': current step


%% Iterate step size coefficient
%mft >= (1 - alpha*tau) * mft0

while mft >=  mft0;    %stopping conditions
    fprintf('Line search iteration: %i \n',iarm+1)
    %% Apply the three point parabolic model.
    tau = sigma*tau;
%     if iarm == 0
%         tau = sigma1*tau;
%     else
%         tau = parab3p(tauc, taum, mft0_2, mftc_2, mftm_2);
%     end

    %% Update vv2.coeff; keep the books on tau.
    step = mft0/norm(vv.agrad(:),1) * vv.agrad / norm(vv.agrad(:),1);
    vv2.acoeff = vv.acoeff - tau*step; %stepping from the original acoeff.
    taum = tauc;
    tauc = tau;

    %% Calculate misfit based on current step size
    vv2.C = idct2(vv2.acoeff);      %Calculate new slipperiness field
    [ii] = ism_sia(aa.s,aa.h,vv2.C, pp,gg,oo);  %Calculate corresponding velocities 
    vv2.u = ii.u; vv2.v = ii.v;                 %SIA                                        
    [vv2] = ism_sstream(vv2,aa,pp,gg,oo );      %SSA 
    [mft] = ism_inversion_misfit(vv2.u,vv2.v,aa,pp,gg, oo); %Current misfit
    
    %% Keep the books on the function norms.
    mftm_2 = mftc_2;
    mftc_2 = mft*mft;
    iarm = iarm+1;
    if iarm > maxarm
        disp(' Armijo failure, too many reductions ');
        armflag = 1;
        return;
    end
    
    disp(iarm)
end
fprintf('Misfit: %i \n', mft)
disp('Completed')
vv2.tau = tau; %Save tau as initial starting point for subsequent line search

end   


%
function taup = parab3p(tauc, taum, ff0, ffc, ffm)
% Apply three-point safeguarded parabolic model for a line search.
%
% C. T. Kelley, April 1, 2003
%
% This code comes with no guarantee or warranty of any kind.
%
% function taup = parab3p(tauc, taum, ff0, ffc, ffm)
%
% input:
%       tauc = current steplength
%       taum = previous steplength
%       ff0 = value of \| F(x_c) \|^2
%       ffc = value of \| F(x_c + \tauc d) \|^2
%       ffm = value of \| F(x_c + \taum d) \|^2
%
% output:
%       taup = new value of tau given parabolic model
%
% internal parameters:
%       sigma0 = .1, sigma1 = .5, safeguarding bounds for the linesearch
%

%
% Set internal parameters.
%
sigma0 = .1; sigma1 = .5;
%
% Compute coefficients of interpolation polynomial.
%
% p(tau) = ff0 + (c1 tau + c2 tau^2)/d1
%
% d1 = (tauc - taum)*tauc*taum < 0
%      so, if c2 > 0 we have negative curvature and default to
%      taup = sigam1 * tau.
%
c2 = taum*(ffc-ff0)-tauc*(ffm-ff0);
if c2 >= 0
    taup = sigma1*tauc; return
end
c1 = tauc*tauc*(ffm-ff0)-taum*taum*(ffc-ff0);
taup = -c1*.5/c2;
if taup < sigma0*tauc, taup = sigma0*tauc; end
if taup > sigma1*tauc, taup = sigma1*tauc; end
end