function y=growth_droop(SS,SS_min)
%compute the growth factor based on droop kinetics
%SS is the substrate vector

%the Liebig's law of minimum is applied
y=min(1.0-SS_min./SS);

y=max([y,0]);%make sure it is no smaller than zero

