
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


vv2 = vv;

%% Construct Basal Slip Field
vv2.alpha = ism_alpha_field(acoeff,vv2, pp, gg, oo);
vv2.Cb = ism_slidinglaw(vv2.alpha,vv2.uv,vv2.Cb,vv2.F2,vv2,aa,pp,gg,oo);

Cb = reshape(gg.S_h'*vv2.Cb,gg.nJ,gg.nI);

% if oo.hybrid, 
%     vv.Cb = ism_cslip_field(vv, pp, gg, oo);    %Construct basal slipperiness
%     C = vv.Cb;
% else
%     vv.C = ism_cslip_field(vv, pp, gg, oo); 
%     C = vv.C;
% end

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
oo.pic_iter = 25; 
[vv2, ~] = ism_deism(vv2,aa,pp,gg,oo );           

if oo.hybrid
F1 = ism_falpha(1,vv2.uv,vv2.nEff_lyrs,vv2,aa,pp,gg,oo );          %Calculate F alpha factors [Hybrid]
F2 = ism_falpha(2,vv2.uv,vv2.nEff_lyrs,vv2,aa,pp,gg,oo );
cst = ism_inv_cost(vv2.uv,vv2.Cb,vv2.alpha,F1,F2,vv2,aa,pp,gg, oo);  %Current misfit

else
cst = ism_inv_cost(vv2.uv,vv2.Cb,vv2.alpha,[],[],vv2,aa,pp,gg, oo);  %Current misfit
end

%% Cost


%% Return gradient if optimization routine requests it

if nargout > 1 % gradient required
    
    %% Depth Integrated Model
    %Second DEISM call, save parameters for adjoint method
    oo.adj_iter = 1;                                     %Save picard iterations
    oo.pic_iter = 5;                                        %Reduced number of picard iteration
    
    vv2_orig = vv2;
    [vv2, rr] = ism_deism(vv2,aa,pp,gg,oo ); 
    
    oo.pic_iter = 1;
    
    if oo.hybrid, 
    C = vv2.Cb;                                                      %Use the correct Basal Drag
    F1 = ism_falpha(1,vv2.uv,vv2.nEff_lyrs,vv2,aa,pp,gg,oo );          %Update F1/F2 w/ new velocities/viscosities [Hybrid]
    F2 = ism_falpha(2,vv2.uv,vv2.nEff_lyrs,vv2,aa,pp,gg,oo );
    else C = vv2.Cb; F1 = zeros(gg.nha,1); F2 = F1; end;
    
    %% Adjoint method with automatic differentiation
    rr.runalpha = [];
    
    % Adjoint variable Uf*
    U_adi = struct('f', vv2.uv, 'dU',ones(gg.nua+gg.nva,1));
    UAD = ism_inv_cost_ADu(U_adi,C,vv2.alpha,F1,F2,vv2,aa,pp,gg,oo);
    rr.adjU = UAD.dU;
    
    %%Determine adjoint state of Cb_final* 
    C_adi = struct('f', C, 'dC',ones(gg.nha,1));
    UAC = ism_inv_cost_ADc(vv2.uv,C_adi,vv2.alpha,F1,F2,vv2,aa,pp,gg,oo);
    rr.adjC = UAC.dC;
    
    %Determine adjoint state of Alpha_final* 
    %Part 1
    C_alpha = ism_slidinglaw_dalpha([],vv2.uv,rr.Cbn(:,end-1),rr.F2n(:,end),vv2,aa,pp,gg,oo);
    C_alpha = spdiags(C_alpha,0,gg.nha,gg.nha);
    rr.adjalpha1 = C_alpha'*rr.adjC;
    
    %Part 2
    alpha_adi = struct('f', vv2.alpha, 'dalpha',ones(gg.nha,1));
    AAD = ism_inv_cost_ADa(vv2.uv,C,alpha_adi,F1,F2,vv2,aa,pp,gg,oo);
    rr.adjalpha2 = AAD.dalpha;
    
    %Combine
    rr.runalpha =  rr.adjalpha1 + rr.adjalpha2;
    
    %Determine adjoint variable Alpha_initial*
    rr = ism_adjAD_main(vv2,rr,aa,pp,gg,oo );
    
    %Gradient of cst function w.r.t acoeff
    gradN = (rr.runalpha.*exp(acoeff(:)));
    
    %% For Manual Testing
%     gradN = gradN/max(abs(gradN));
%     imagesc(reshape(gg.S_h'*gradN,gg.nJ,gg.nI))
%     acoeff_orig = acoeff;
%     cst_orig = cst;
%     acoeff = acoeff_orig - gradN;

end

end   

