function km=mm_km(km0, vmax,dd)
%compute substrate affinity,scaled with
% maximum surface sites and the diffusion rate
km=km0*(1+dd*vmax/km0);