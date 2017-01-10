function [duv] = adi_deism_backbone_deriv(uv, C,nEff,aa,pp,gg,oo )
%% Shallow Stream Model 
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:


[LHS] = ism_deism_fieldeq(C.f,nEff.f, aa,pp,gg,oo);              %Field Equations
[LHS_C] = ism_deism_fieldeq(zeros(gg.nha,1),nEff.dacoeff, aa,pp,gg,oo);              %Field Equations
[LHS_nEff] = ism_deism_fieldeq(C.dacoeff,zeros(gg.nha,1), aa,pp,gg,oo);              %Field Equations

dLHS_dacoeff = LHS_C + LHS_nEff;

duv = LHS\(-dLHS_dacoeff*uv);



end

