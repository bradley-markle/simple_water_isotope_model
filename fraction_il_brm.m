function [fraction_i, fraction_l] = fraction_il_brm(T,T_hi,T_lo)
% This is a function modified form the Kavanaugh isotope model, BRM 2016 
%input T should be in degC

% FRACTION__IL calculates the ice fraction of precipitation 
% for a given temperature T (expressed in Kelvin).  The liquid 
% fraction is thus (1 - f_i).

% The following has f_i = 1 at T <= -40 C and f_i = 0 at T> = -5C 
% and varying linearly between these two temperatures. 

% T_C = T - 273.15;
% T_lo = -30.0;
% T_hi = -10.0;

%if T_C <= T_lo
%  value = 1.0;
%elseif T_C >= T_hi
%  value = 0.0;
%else
%  value = 1-(T_C-T_lo)/(T_hi-T_lo);
%end

dT    = T_hi-T_lo;
T_0   = (T_lo+T_hi)/2;
arg   = (T-T_0)/(dT/4);

fraction_i = 0.5*(1-erf(arg));
fraction_l= 1-fraction_i;


% 
% %below is example of how to use it to calculate fractionation factors
% 
% T_C  = T-273.15;
% T_lo = -15.0;
% T_hi = -5.0;
% %F_lo = 1.06;  % i.e. value at low temperatures
% %F_hi = 1.02;  % i.e. value at high temperatures
% F_lo = 1.00;  % i.e. value at low temperatures
% F_hi = 1.00;  % i.e. value at high temperatures
% 
% if T_C <= T_lo
%   factor = F_lo;
% elseif T_C >= T_hi
%   factor = F_hi;
% else
%   factor = F_lo - (F_lo-F_hi)*(T_C-T_lo)/(T_hi-T_lo);
% end
% 
% factor = 1.0;
% 
% value1 = fraction_il(T)*alpha_D_iv(T) + (1-fraction_il(T))*alpha_D_lv(T);
% value2 = value1-1.0;
% value  = 1.0+factor*value2;