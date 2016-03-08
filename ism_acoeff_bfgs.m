
function [vv2] = ism_acoeff_bfgs(vv,aa, pp, gg, oo)
%% Line search to update inversion alpha coefficients 
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

 
disp('BFGS line search')

if ~isfield(vv,'aHessInv'), vv.aHessInv = speye(gg.nIJ); end;    %Acoeff Hessian
if ~isfield(vv,'bfgs_iter'), vv.bfgs_iter = 1; end;    %Acoeff Hessian
if ~isfield(oo,'bfgs_ls'), oo.bfgs_ls = 2; end;    %Acoeff Hessian


%% Paramaters
vv2 = vv;               

tau = 1;                %Initial step size
sigma = .25;            %Step size reduction
alpha = 10^-4;          %Armijo condition coefficient.
maxarm = 15;            %Maximum number of steps (reductions)
iarm = 0;               %Counter of number of reductions
armflag = 0;            %Flag for exceeding max number of reductions

[mft0] = ism_inv_misfit(vv,aa,pp,gg, oo);  %Initial velocity misfit

%% Determine step 
if vv.bfgs_iter < oo.bfgs_ls 
[vv] = ism_adjoint(vv,aa,pp,gg,oo );
[vv] = ism_cslip_grad(vv, pp, gg, oo);
sgrad = -vv.agrad(:);
sgrad = reshape(sgrad,gg.nJ,gg.nI);
vv2.bfgs_iter = vv2.bfgs_iter+1;
tau = 1/(norm(sgrad(:),2));
%tau = 1;
else
disp('Calculating sgrad')
tic
sgrad = -vv.aHessInv*vv.agrad(:);
sgrad = reshape(sgrad,gg.nJ,gg.nI);
toc
disp('Done')
tau = 1;
end

sdir = sgrad./(norm(sgrad,2));


%% Initial step
vv2.acoeff = vv.acoeff + tau*sgrad;
vv2.C = exp(idct2(vv2.acoeff));       %Calculate new slipperiness field
vv2.C = vv2.C;

[vv2] = ism_sia(aa.s,aa.h,vv2.C,vv2, pp,gg,oo); %SIA                                        
[vv2] = ism_sstream(vv2,aa,pp,gg,oo );          %SSA 

[mft] = ism_inv_misfit(vv2,aa,pp,gg, oo); %Current misfit

%% Iterate step size coefficient
while mft >= mft0 + alpha * tau * (vv.agrad(:)'*sdir(:));    %stopping conditions
    fprintf('Line search iteration: %i \n',iarm+1)
    %% Step Size
    tau = sigma*tau;

    %% Update vv2.coeff; keep the books on tau.
    vv2.acoeff = vv.acoeff + tau*sgrad; %stepping from the original acoeff.

    %% Calculate misfit based on current step size
    vv2.C = exp(idct2(vv2.acoeff));      %Calculate new slipperiness field
    [vv2] = ism_sia(aa.s,aa.h,vv2.C,vv2,pp,gg,oo);  %Calculate corresponding velocities 
    [vv2] = ism_sstream(vv2,aa,pp,gg,oo );      %SSA 
    [mft] = ism_inv_misfit(vv2,aa,pp,gg, oo); %Current misfit
    
    %% Keep the books on the function norms.
    iarm = iarm+1;
    if iarm > maxarm
        disp('Line Search Step Failure, too many reductions ');
        armflag = 1;
        vv2.armflag = 1;
        break
    end
    
    disp(iarm)
end
disp('Done line search')
fprintf('Inversion Misfit: %i \n', mft)

%Update values for new acoeff
[vv2] = ism_adjoint(vv2,aa,pp,gg,oo );
[vv2] = ism_cslip_grad(vv2, pp, gg, oo);            %calculate agrad for current alpha coeff

%% Approximate new Hessian Inverse
if (vv2.bfgs_iter >= oo.bfgs_ls) & ~armflag
disp('Updating Approximated Hessian Inverse')
tic
S = vv2.acoeff(:) - vv.acoeff(:);
Y = vv2.agrad(:) - vv.agrad(:);
p = (S'*Y);
%p = 1/(Y'*S);
%vv2.aHessInv = vv.aHess + (Y*Y')/(Y'*S) - (vv.aHess*S*S'*vv.aHess)/(S'*vv.aHess*S);
%vv2.aHessInv = (speye(gg.nIJ) - p*S*Y')*vv.aHessInv*(speye(gg.nIJ) - p*S*Y') + p*S*Y';
vv2.aHessInv = vv.aHessInv +...
    ((S'*Y +Y'*vv.aHessInv*Y)*(S*S'))/((S'*Y)^2) -...
    (vv.aHessInv*(Y*S') + (S*Y')*vv.aHessInv)/(S'*Y);

toc
disp('Done')
datestr(clock)
end

end   

