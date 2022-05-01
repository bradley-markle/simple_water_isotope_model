function [ss] = mixed_phased_supersaturation(T,P0,fraction_i, fraction_l)

% Example routine for calculating saturated pseudo adiabatic change
% follows iterative routine from Bakhshaii and Stull 2013
%
% Thus, the iterative solution proceeds as follows. For
% any starting [P1, T(P1)] such as at an initial point (Pi, Ti) =
% (100 kPa, uw) on the saturated adiabat, solve Eq. (7) for
% es, then solve Eq. (6) for rs, then solve Eq. (2) for b (using
% rs instead of r) or use the approximation b ? 1, then
% solve Eq. (8) for dT/dP. Next, assume that the saturated
% adiabat is approximately linear over a very small increment
% of pressure ?P = P2 - P1, giving T(P2) ~= T(P1) +
% (P2 - P1)*dT/dP. This new temperature and pressure can
% be used as input to the next small increment, with the
% process repeated until the final destination pressure Pf is
% reached.

%Eq. (7): e_s = e_0 exp( (17.67 *(T - 273.15)) / (T - 29.65)); T is in K

% except now we are going to use a different set of equations for the
% saturation vapor pressure that take into acount saturation over ice for
% T<0. these are from: %Murphy and Koop 2005 saturated vapor pressures for
% ice, eq (7) and (10), note that these give es in Pa not kPa!! add factor
% 1000 to below eqs
% es_MK(T>0) = exp(54.842763 - 6763.22 ./ TK(T>0) - 4.21 .* log(TK(T>0)) + 0.000367 .* TK(T>0) +...
%     tanh(0.0415 .* (TK(T>0) - 218.8)) .*  (53.878 - 1331.22 ./ TK(T>0) - 9.44523 .* log(TK(T>0)) + 0.014025 .* TK(T>0))) ;%with T in [K] and ew in [Pa]
% es_MK(T<=0) = exp(9.550426 - (5723.265./TK(T<=0)) + 3.53068.*log(TK(T<=0)) - 0.00728332 .*TK(T<=0)) ;%T in [K] and ei in [Pa]


%Eq. (6): r_s = eta * e_s / (P - e_s);

%Eq. (2): b = [1 + (r/eta)]/[1 + (r/c)] ~= (1 - 0.24 *r); however we
%replace r with r_s; or b ~=1

%Eq. (8): dT/dP = (b/P) * (R_d*T + L_v*R_s) / (C_pd + ((L_v^2)*r_s*eta*b)/(R_d*(T^2)));

%Ok the above calculates temp for a given pressure gradient. We want to
%calculate pressure for a given temp gradient. So we just need to do the
%inverted calculations of the last step.
% ?T = T2 - T1, P(T2) ~= P(T1) + (T2 - T1)*dP/dT
%Eq. (8) becomes: dP/dT = (P/b) *  (C_pd + ((L_v^2)*r_s*eta*b)/(R_d*(T^2))) / (R_d*T + L_v*r_s);

% Temperature profile.
% T = 10:-0.1:-40;%deg C
TK = T + 273.15;%Temp in K 
% Tcrit = 0;
% Pressure profile
P = nan(size(T));
% P0 = 100;%kPa, initial pressure
% P0 = 101.325;%kPa, initial pressure

%other profiles
e_s = nan(size(T));
r_s = nan(size(T));
b = nan(size(T));
dPdT=nan(1,length(T)-1);

%set up parameters;
e_0 = 0.6112; %kPa, at 0deg C
L_v = 2.501 * 10^6; %J kg^-1, at 0degC, approximate
R_d = 287.053; %J K^-1 kg^-1
C_pd = 1004; %J K^-1 kg^-1, at 0 degC, approximate
eta = 0.622; %kg_watervapor / kg_dryair = R_d/R_v; ratio of ideal gas constants
c = 0.5427; %kg_watervapor kg_dryair = Cpd/Cpv; ratio of specific heats at constant pressure at 0 degC, approximate 
C_pv = C_pd/c; %J K^-1 kg^-1, at 0 degC, approximate; checked this against Fig. 1 in MK05, looks good.

%Calculate Latent heat of ice, from Murphy and Koop 2005, eg (5), T > 30K, with L_i in J mol?1.
L_i =46782.5 + 35.8925.*TK - 0.07414.*TK.^2 + 541.5.*exp(-(TK ./123.75).^2);
L_i = L_i.*55.5084350618;%to convert J mol^?1 to J kg^-1; 1kg water = 55.508 mol
%calculate specific heat of ice
C_pi = -2.0572 + 0.14644.*TK + 0.06163.*TK .*exp(-(TK./125.1).^2);% T> 20 K., MK05 eg (4)
c_i=C_pd./C_pi;%ratio of specifc heats for ice and dry air

%calculate fraction of ice and supercooled water
% T_hi=0;
% T_lo=-30;
% [fraction_i, fraction_l] = fraction_il_brm(T,T_hi,T_lo);
L_eff=fraction_l.*L_v.*ones(size(T)) + fraction_i.*L_i;
c_eff=C_pd./(fraction_l.*C_pv.*ones(size(T)) + fraction_i.*C_pi);

%initialize variables
P(1) = P0;
% e_s(1) = e_0 * exp( (17.67 *(T(1) - 273.15)) / (T(1) - 29.65));
% r_s(1) = eta * e_s(1) / (P(1) - e_s(1));
% b(1) = (1 + (r_s(1)/eta))/(1 + (r_s(1)/c));
% dPdT(1) = (P(1)/b(1)) *  (C_pd + ((L_v^2)*r_s(1)*eta*b(1))/(R_d*(T(1)^2))) / (R_d*T(1) + L_v*r_s(1));

e_s_l = (1000^-1).*exp(54.842763 - 6763.22 ./ TK - 4.21 .* log(TK) + 0.000367 .* TK +...
    tanh(0.0415 .* (TK - 218.8)) .*  (53.878 - 1331.22 ./ TK - 9.44523 .* log(TK) + 0.014025 .* TK)) ;%with T in [K] and ew in [kPa]
e_s_i = (1000^-1).*exp(9.550426 - (5723.265./TK) + 3.53068.*log(TK) - 0.00728332 .*TK) ;%T in [K] and ei in [kPa] 

%mix ice and water
e_s=fraction_l.*e_s_l + fraction_i.*e_s_i;

for i = 1:length(T)-1
% e_s(i) = e_0 * exp( (17.67 *(TK(i) - 273.15)) / (TK(i) - 29.65));
%Murphy and Koop 2005 saturated vapor pressures for ice, eq (7) and (10)

r_s_l(i) = eta * e_s_l(i) / (P(i) - e_s_l(i));
r_s_i(i) = eta * e_s_i(i) / (P(i) - e_s_i(i));
r_s(i) = eta * e_s(i) / (P(i) - e_s(i));
b(i) = (1 + (r_s(i)/eta))/(1 + (r_s(i)/c_eff(i)));
dPdT(i) = (P(i)/b(i)) *  (C_pd + ((L_eff(i).^2)*r_s(i)*eta*b(i))/(R_d*(TK(i).^2))) / (R_d*TK(i) + L_eff(i)*r_s(i));

P(i+1)= P(i) + dPdT(i) *(TK(i+1)-TK(i));
end
% e_s(end) = (1000^-1).*exp(9.550426 - (5723.265./TK(i)) + 3.53068.*log(TK(i)) - 0.00728332 .*TK(i)) ;%T in [K] and ei in [kPa] factor of 1000 needed to turn Pa into kPa
r_s(length(T)) = eta .* e_s(end) ./ (P(end) - e_s(end));
r_s_i(length(T)) = eta .* e_s_i(end) ./ (P(end) - e_s_i(end));
r_s_l(length(T)) = eta .* e_s_l(end) ./ (P(end) - e_s_l(end));
b(end) = (1 + (r_s(end)/eta))/(1 + (r_s(end)/c_eff(end)));
ss=r_s./r_s_i;

% ss(T>0)=1;
f=r_s/r_s(1);
