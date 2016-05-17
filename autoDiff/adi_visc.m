

z = adigatorCreateDerivInput([18327,1], 'z');
U = adigatorCreateAuxInput([36984,1]);
nEff = adigatorCreateAuxInput([18327,1]);
C = adigatorCreateAuxInput([18327,1]);


adigator('ism_visc', {z,U,nEff,C,aa,pp,gg,oo},'myderiv')





%Testing Routine
% U = [u;v];
% z = gg.S_h*aa.s(:);
% z = gg.S_h*aa.b(:);
% C = gg.S_h*aa.C(:);
% ism_visc(z,U,nEff,C,aa,pp,gg,oo);