
function [cst,gradN] = ism_adjAD_optWrapper(acoeff,vv,aa, pp, gg, oo)
% Inputs:
%   acoeff  alpha coefficients determining basal cslip (C)
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   cst     cst function of misfit between observed and predicted velocities
%   grad    gradient of cost function w.r.t acoeff

vv_orig = vv;
acoeff_orig = acoeff;
vv2 = vv;

%% Construct Basal Slip Field
vv2.alpha = ism_alpha_field(acoeff,vv2, pp, gg, oo);
vv2.Cb = ism_slidinglaw(vv2.alpha,vv2.uv,vv2.Cb,vv2.F2,vv2,aa,pp,gg,oo);

Cb = reshape(gg.S_h'*vv2.Cb,gg.nJ,gg.nI);

%% Initial velocites and viscosity

[vv2] = ism_sia(aa.s,aa.h,Cb,vv2, pp,gg,oo);    %SIA

if oo.hybrid,                               %Setup Hybrid Approximation                           
Cb = vv2.Cb;                                  %Basal Slipperiness
uv = vv2.uv;
nEff = ism_visc(uv,vv2,aa,pp,gg,oo);          %Initial viscosity;
nEff_lyrs = repmat(nEff,1,oo.nl+1);

F2 = ism_falpha(2,uv,nEff_lyrs,vv2,aa,pp,gg,oo );
C = Cb(:)./(1 + (pp.c13*Cb(:)).*(F2));                   %Effective Basal Slipperiness

vv2.C = C;
vv2.F2 = F2;
vv2.nEff = nEff;
vv2.nEff_lyrs = nEff_lyrs;

else
uv = vv2.uv;
nEff = ism_visc(uv,vv2,aa,pp,gg,oo);

vv2.C = vv2.Cb;
vv2.nEff = nEff;
end



%% DEISM

oo.adj_iter = 0;                                     %Depth Integrated Model
oo.pic_iter = 50; 
[vv2, ~] = ism_deism(vv2,aa,pp,gg,oo );           


cst = ism_inv_cost(vv2.uv,vv2.Cb,vv2.alpha,vv2.F1,vv2.F2,vv2,aa,pp,gg, oo);  %Current misfit

%% Cost
cst_orig = cst;
vv2_orig = vv2;

%% Return gradient if optimization routine requests it

if nargout > 1 % gradient required
    
    %% Depth Integrated Model
    %Second DEISM call, save parameters for adjoint method
    oo.adj_iter = 1;                                     %Save picard iterations
    oo.pic_iter = 2;                                        %Reduced number of picard iteration
    
    [vv2, rr] = ism_deism(vv2,aa,pp,gg,oo ); 
       
    [ gradN ] = ism_adjAD_grad2(acoeff,vv2,rr,aa,pp,gg,oo);
    
    %% For Manual Testing
%     gradN = gradN/max(abs(gradN));
%     imagesc(reshape(gg.S_h'*gradN,gg.nJ,gg.nI))
%     acoeff_orig = acoeff;
%     cst_orig = cst;
%     acoeff = acoeff_orig - gradN;

 end

end   

