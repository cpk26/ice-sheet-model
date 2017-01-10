function [ nEff ] = adi_deism_nEff(alpha,uv_init,C_init,F2_init,vv,aa,pp,gg,oo )


Cb = adi_slidinglaw(alpha,uv_init,C_init,F2_init,vv,aa,pp,gg,oo);

evFac = 1;
if oo.hybrid
evFac = (1 + (pp.c13*Cb).*F2_init);
end

C = Cb./evFac;

[nEff] = adi_visc_di(uv_init,C,aa,pp,gg,oo); %Updated Viscosity


end

