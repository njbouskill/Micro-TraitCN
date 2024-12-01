function vmax=vmax_tscaled(vmax0, temp, tmin, tmax)
%modeling temperature dependence of vmax, parabolic function of temperature
%dependence, this assumes there is an optimal temperature at which the
%activation and deactivation of the transporters are equal.

vmax=vmax0*(temp-tmin)*(temp-tmax);