function [] = ism_adjAD_generate( vv,aa, pp, gg, oo )
%Automatically differentiate necessary functions using Adigator, specifying
%input sizes from current simulation

%% Clear working directory of previously generated AD function
% fn = {'ism_slidinglaw_ADu','ism_slidinglaw_ADa', 'ism_inv_cost_ADu', 'ism_inv_cost_ADc', 'ism_inv_cost_ADa','ism_visc_AD',...
%     'ism_visc_ADu','ism_visc_diAD','ism_visc_diSADu','ism_visc_diSADc','ism_dim_Ddiag_ADc',...
%     'ism_dim_Ddiag_ADnEff','ism_cslip_form_ADc', 'ism_visc_diSAD','adi_deism_ADu','adi_deism_ADa'};
% for j = 1:numel(fn)
%     tmpA = char(strcat(fn(j), '.m')); tmpB = char(strcat(fn(j), '.mat'));
%     if exist(tmpA, 'file'), delete(tmpA); end; 
%     if exist(tmpB, 'file'), delete(tmpB); end; 
% end

%% Apply Automatic Differentiation
%Diff w.r.t acoeff

% acoeff = adigatorCreateDerivInput([gg.nha,1],'acoeff');
% %uv = adigatorCreateDerivInput([gg.nua+gg.nva,1], 'uv');
% SS = adigatorCreateDerivInput([1,1], 'SS');
% uv = adigatorCreateAuxInput([gg.nua+gg.nva,1]);
% Cb = adigatorCreateAuxInput([gg.nha,1]);
% F2 = adigatorCreateAuxInput([gg.nha,1]);
% oo.slidinglaw = 1;
% oo.inv_cst=1;
% oo.hybrid=1;
% adigator('adi_deism', {acoeff,uv,SS,Cb,F2,vv,aa,pp,gg,oo },'adi_deism_ADa')
% clear acoeff uv Cb F2;
%     
options = adigatorOptions('OVERWRITE',1);


%Sliding Law w.r.t velocity
uv = adigatorCreateDerivInput([gg.nua+gg.nva,1], 'U');
alpha = adigatorCreateAuxInput([gg.nha,1]);
Cb = adigatorCreateAuxInput([gg.nha,1]);
F2 = adigatorCreateAuxInput([gg.nha,1]);
adigator('ism_slidinglaw', {alpha,uv,Cb,F2,vv,aa,pp,gg,oo},'ism_slidinglaw_ADu',options)
clear U Cb alpha F2;

%Sliding Law w.r.t alpha
alpha = adigatorCreateDerivInput([gg.nha,1], 'alpha');
uv = adigatorCreateAuxInput([gg.nua+gg.nva,1]);
Cb = adigatorCreateAuxInput([gg.nha,1]);
F2 = adigatorCreateAuxInput([gg.nha,1]);
adigator('ism_slidinglaw', {alpha,uv,Cb,F2,vv,aa,pp,gg,oo},'ism_slidinglaw_ADa',options)
clear U Cb alpha F2;

%Cost Function w.r.t velocity
U = adigatorCreateDerivInput([gg.nua+gg.nva,1], 'U');
C = adigatorCreateAuxInput([gg.nha,1]);
alpha = adigatorCreateAuxInput([gg.nha,1]);
F1 = adigatorCreateAuxInput([gg.nha,1]);
F2 = adigatorCreateAuxInput([gg.nha,1]);
adigator('ism_inv_cost', {U,C,alpha,F1,F2,vv,aa,pp,gg,oo},'ism_inv_cost_ADu',options)
clear U C alpha F1 F2;

%Cost Function w.r.t Basal Slipperiness
C = adigatorCreateDerivInput([gg.nha,1],'C');
U = adigatorCreateAuxInput([gg.nua+gg.nva,1]);
alpha = adigatorCreateAuxInput([gg.nha,1]);
F1 = adigatorCreateAuxInput([gg.nha,1]);
F2 = adigatorCreateAuxInput([gg.nha,1]);
adigator('ism_inv_cost', {U,C,alpha,F1,F2,vv,aa,pp,gg,oo},'ism_inv_cost_ADc',options)
clear U C alpha F1 F2;

%Cost Function w.r.t Alpha
alpha = adigatorCreateDerivInput([gg.nha,1],'alpha');
U = adigatorCreateAuxInput([gg.nua+gg.nva,1]);
C = adigatorCreateAuxInput([gg.nha,1]);
F1 = adigatorCreateAuxInput([gg.nha,1]);
F2 = adigatorCreateAuxInput([gg.nha,1]);
adigator('ism_inv_cost', {U,C,alpha,F1,F2,vv,aa,pp,gg,oo},'ism_inv_cost_ADa',options)
clear U C alpha F1 F2;


%Viscosity w.r.t velocity
U = adigatorCreateDerivInput([gg.nua+gg.nva,1], 'U');
adigator('ism_visc', {U,vv,aa,pp,gg,oo},'ism_visc_ADu',options)
clear U;

if oo.hybrid
%Depth Integrated Viscosity w.r.t velocity
U = adigatorCreateDerivInput([gg.nua+gg.nva,1], 'U');
nEff_lyrs = adigatorCreateAuxInput([gg.nha,oo.nl+1]);
C = adigatorCreateAuxInput([gg.nha,1]);
adigator('ism_visc_diS', {U,nEff_lyrs,C,aa,pp,gg,oo},'ism_visc_diSADu',options)
clear U nEff C;

%Depth Integrated Viscosity w.r.t C
C = adigatorCreateDerivInput([gg.nha,1], 'C');
U = adigatorCreateAuxInput([gg.nua+gg.nva,1]);
nEff_lyrs = adigatorCreateAuxInput([gg.nha,oo.nl+1]);
adigator('ism_visc_diS', {U,nEff_lyrs,C,aa,pp,gg,oo},'ism_visc_diSADc',options)
clear U nEff C;
end

%Diagonal of A matrix w.r.t Basal slipperiness
C = adigatorCreateDerivInput([gg.nha,1], 'C');
nEff = adigatorCreateAuxInput([gg.nha,1]);
adigator('ism_dim_Ddiag', {C,nEff,aa,pp,gg,oo},'ism_dim_Ddiag_ADc',options)
clear C nEff;

%Diagonal of A matrix w.r.t depth integrated viscosity
nEff = adigatorCreateDerivInput([gg.nha,1], 'nEff');
C = adigatorCreateAuxInput([gg.nha,1]);
adigator('ism_dim_Ddiag', {C,nEff,aa,pp,gg,oo},'ism_dim_Ddiag_ADnEff',options)
clear C nEff;

%Basal Slipperiness form w.r.t Basal Slipperiness
C = adigatorCreateDerivInput([gg.nha,1],'C');
flag = adigatorCreateAuxInput([1,1]);
F2 = adigatorCreateAuxInput([gg.nha,1]);
adigator('ism_cslip_form', {flag,F2,C,aa,pp,gg,oo },'ism_cslip_form_ADc',options)
clear U C F1 F2;




end

