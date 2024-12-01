function y=substrate_kinetics(s,par)

if(strcmp(par.kinetics,'monod1'))
    %monod kinetics
    y=s./(par.km+s);
elseif(strcmp(par.kinetics,'monod2'))
    %difussion limited monod kinetics
    y=s./((par.km+par.vmax*par.km1)+s);
elseif(strcmp(par.kinetics,'haldane'))
    %haldane kinetics
    y=s./(par.km+s*(1+s./par.ki));
elseif(strcmp(par.kinetics,'haldanedifl'))
    %haldane kinetics + diffusion limiting
    y=s./((par.km+par.vmax*par.ri)+s*(1+s./par.ki));
end
