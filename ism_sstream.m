function [vv2, FE1, FE2] = ism_sstream(vv,aa,pp,gg,oo )
%% Shallow Stream Model 
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv2     struct containing new solution variables
%   J       Jacobian matrix

%% Variables (Non-Dimensionalized)
s = aa.s; s_diag = spdiags(s(:),0,gg.nIJ,gg.nIJ);      %Topography
h = aa.h; h_diag = spdiags(s(:),0,gg.nIJ,gg.nIJ);
%Use gradient instead of gg.nddx/y 
%since periodic BC conditions do not apply
[Sx,Sy] = gradient(s, gg.dx, gg.dy);        
Sx = Sx(:); Sy = Sy(:); h = h(:);

u = vv.u;                                   %Velocities
v = vv.v;             

exx = gg.du_x*u;                            %Strain Rates
eyy = gg.dv_y*v;
exy = 0.5*(gg.du_y*u + gg.dv_x*v);

edeff = sqrt(exx.^2 + eyy.^2 + exx.*eyy + exy.^2 + pp.n_rp.^2);

n = pp.n_Glen;                              %Ice Flow 
B = pp.c4;
C = aa.C;
nEff = 0.5 * B * edeff.^(1-n)/n; nEff_diag = spdiags(nEff(:),0,gg.nIJ,gg.nIJ);                                  

%% Field equations for velocities

A1 = gg.S_u*gg.dh_x*2*nEff_diag*h_diag*gg.du_x*(gg.S_u');




TMP1 = 2*h.*vEff.*exx + h.*vEff.*eyy; TMP1 = gg.nddx*TMP1;
TMP2 = h.*vEff.*exy; TMP2 = gg.nddy*TMP2;
[TMP3, ~] = slidingLaw(u,v,C,Sx,Sy,oo);
TMP4 = -pp.c5*h.*Sx;
FE1 = TMP1+TMP2+TMP3-TMP4;


TMP5 = 2*h.*vEff.*eyy + h.*vEff.*exx; TMP5 = gg.nddy*TMP5;
TMP6 = h.*vEff.*eyx; TMP6 = gg.nddx*TMP6;
[~, TMP7] = slidingLaw(u,v,C,Sx,Sy,oo);
TMP8 = -pp.c2*h.*Sy;
FE2 = TMP5+TMP6+TMP7-TMP8;

vv2=vv;
end


function [tbx, tby] = slidingLaw()

switch oo.sL
    case 'ismip'
        tbx = C.*u;
        tby = C.*v;
        
    case 'weertman'
        %To be implemented...
        %w = delta*(u.*Sx + v.*Sy);
        %Ub = sqrt(u.^2 + v.^2 + w.^2);
end


end

