function [f, P, e_s, r_s, rh] = isobaric_func_rh0(T,P0,rh0)

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
% Pressure profile
P = nan(size(T));
% P0 = 101.325;%kPa, initial pressure
% P0 = 100;%kPa, initial pressure

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

%initialize variables
% P(1) = P0;
P(:) = P0;

% e_s(1) = e_0 * exp( (17.67 *(T(1) - 273.15)) / (T(1) - 29.65));
% r_s(1) = eta * e_s(1) / (P(1) - e_s(1));
% b(1) = (1 + (r_s(1)/eta))/(1 + (r_s(1)/c));
% dPdT(1) = (P(1)/b(1)) *  (C_pd + ((L_v^2)*r_s(1)*eta*b(1))/(R_d*(T(1)^2))) / (R_d*T(1) + L_v*r_s(1));

%OK make this isobaric instead of pseudo adiabatic
for i = 1:length(T)-1
e_s(i) = e_0 * exp( (17.67 *(TK(i) - 273.15)) / (TK(i) - 29.65));
r_s(i) = eta * e_s(i) / (P(1) - e_s(i));%P always = P0
b(i) = (1 + (r_s(i)/eta))/(1 + (r_s(i)/c));

end
e_s(end) = e_0 * exp( (17.67 *(TK(end) - 273.15)) / (TK(end) - 29.65));
r_s(end) = eta * e_s(end) / (P(1) - e_s(end));
b(end) = (1 + (r_s(end)/eta))/(1 + (r_s(end)/c));

% f=r_s/r_s(1);
f=r_s./(r_s(1)*rh0);
f(f>1)=1;
rh=(r_s(1)*rh0)./r_s;
rh(f<1)=1;
