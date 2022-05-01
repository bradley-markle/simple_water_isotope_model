function [d18O_p, dD_p, d17O_p, d18O_p_ln, dD_p_ln, d17O_p_ln, dxs, dxs_ln, d17Oxs, T, S, supersat, r_s, e_s, f, P, rh, RD_c, RD_p, R18O_c, R18O_p, R17O_c, R17O_p] = distillation_2020(T0, T_end, dT, dD_v0, d18O_v0, d17Oxs_v0, RHn0, rh0, sst0, a, b, c, pathway)
% BRM March 2017, edited September 2017
% [d18O_p, dD_p, d17O_p, d18O_p_ln, dD_p_ln, d17O_p_ln, dxs, dxs_ln, d17Oxs, T, S, supersat, r_s, e_s, f, P,rh, RD_c, RD_p, R18O_c, R18O_p, R17O_c, R17O_p] = distillation_2018(T0, T_end, dT, dD_v0, d18O_v0, d17Oxs_v0, RHn0, rh0, sst0, a, b, c, pathway)

% The user may specify the type of pathway, either pathway = 'adiabatic' or
% 'isobaric'. The default is 'adiabatic', which is a pseudo adiabatic
% pathway, which accounts for mixed ice and liquid water cloud, and is
% consistent with the user defined supersaturation function. It's pretty
% slick.

%% Housekeeping
if ~exist('pathway','var')
    pathway = 'adiabatic';%default is pseudo adiabatic pathway
end

%% set up temperature grid
T = T0:-dT:T_end;
T0K = T0+273.15;
TK = T+273.15;% converts to deg K

%% determine fractions of ice and liquid water in cloud
%This is an old method for determining the mixing of ice and liquid
% T_hi=-10;
% T_hi=-10;
% T_lo=-50;
% T_lo=-35;

%good for b=0.0053,pure molecular diffusion
% T_hi=-20;
% T_lo=-35;
% [fraction_i, fraction_l] = fraction_il_brm(T,T_hi,T_lo);%this uses a KC03-like error function.

% OR ~~~ This is now my preferred method! BRM March 24 2017
% method = 'mid';
% method = 'iir';
method = 'adj';
[fraction_i, fraction_l] = fraction_il_brm_H10(T,method);%this uses polynomial fits based on sattelite data.
clear method
%% Determine f, fraction of water remaining in cloud
P0 = 101.325;%kPa, initial pressure. This is a reasonable surface pressure. 
% P0= 98;

if strcmp(pathway,'adiabatic')
%mixed ice super cooled water cloud
[supersat] = mixed_phased_supersaturation(T,P0,fraction_i, fraction_l);
%pseudo adiabatic pathway with consistent supersaturation

%this version assumes the air parcel reaches saturations immediately.
%This is also equivalent to assuming the relative humidity is constant at rh0 and
%enough moisture is removed at each stage to return (below saturation) to
%this RH.(initSat)
[f, P, e_s, r_s, dPdT,ss,ss_l] = pseudo_adiabat_function(T,P0,fraction_i, fraction_l,supersat,a,b,c);
rh=1.*ones(size(r_s));

% this version cools the parcel from the original rh0 to saturation (initRH) 
% [f, P, e_s, r_s, dPdT,ss,rh] = pseudo_adiabat_function_rh0(T,P0,rh0,fraction_i, fraction_l,supersat,a,b,c);


% this version cools the parcel from the original rh0 to saturation, but then overrains to less than saturation, defined by RH_in (initRH_overrain) 
% RH_in=0.80;
% [f, P, e_s, r_s, dPdT,ss,rh] = pseudo_adiabat_function_rh0_alt(T,P0,rh0,fraction_i, fraction_l,supersat,a,b,c,RH_in);
% clear RH_in

elseif strcmp(pathway,'isobaric')
%isobaric pathway
% [f, P, e_s, r_s] = isobaric_func(T,P0);
[f, P, e_s, r_s, rh] = isobaric_func_rh0(T,P0,rh0);

ss =a-b.*T-c.*(T.^2);
ss(ss<1)=1;
supersat=ss;
end

%% Equilibrium fractionation factors
%there are lots of options in here. It is currently set up for my prefered
%factors.

% a=exp((C1+(10^3)*C2/T +(10^6)*C3/(T^2))/1000);
%Table 3.2, Criss, Chapter 3, calculation of Fractionation Factor, a, for
%various circumstances, D=deuterium, O=oxygen, l=liquid, v=vapor, i=ice
%Type, phase A, phase B, C1, C2, C3, T min (K), T max (K)
% 'D' 'l' 'v' 52.612 -76.248 24.844 273 373;
% 'D' 'l' 'v' -100.0 0 15.013 258 273;
% 'D' 'i' 'v' -94.5 0 16.289 233 273;
% 'O' 'l' 'v' -2.0667 -0.4156 1.137 273 373;
% 'O' 'i' 'v' -28.224 11.839 0 240 273;
%Note: this table uses 273 K (not 273.15 K); we'll assume these are the
%same!

%X is Table 3.2 without "Type, Phase A, and Phase B"
X=[52.612 -76.248 24.844 273 373; -100.0 0 15.013 258 273; -94.5 0 16.289 233 273;...
    -2.0667 -0.4156 1.137 273 373; -28.224 11.839 0 240 273];

%Z is new equilibrium fractionation factor for -40 to 0C from Mads Ellehoj.
%I think these are wrong. Use newer facotrs from Lamb et al.
Z=[0.2133 -203.10 48888 233 273; 0.0831 -49.192 8312.5 233 273];


%OK now we have to account for the fraction of ice and water
    C1D_l=X(1,1); %For D above T>0, from Criss
    C2D_l=X(1,2); 
    C3D_l=X(1,3);
    C1D_l0=X(2,1); %For D below 0, from Criss
    C2D_l0=X(2,2);
    C3D_l0=X(2,3);
    C1O_l=X(4,1); %For O above T>0, from Criss
    C2O_l=X(4,2); 
    C3O_l=X(4,3);
    
    C1D_i=X(3,1); %For D and O below T<0, from Criss
    C2D_i=X(3,2); 
    C3D_i=X(3,3);
    C1O_i=X(5,1);
    C2O_i=X(5,2); 
    C3O_i=X(5,3);
     
%     C1D_i=Z(1,1); %For D and O below T<=0 using Mads Ellehoj Z matrix above
%     C2D_i=Z(1,2); 
%     C3D_i=Z(1,3);
%     C1O_i=Z(2,1);
%     C2O_i=Z(2,2); 
%     C3O_i=Z(2,3);  
%     


     aDe_l=exp((C1D_l+(10^3).*C2D_l./TK + (10^6).*C3D_l./(TK.^2))./1000);
%      aDe_i=exp((C1D_i+(10^3).*C2D_i./TK + (10^6).*C3D_i./(TK.^2))./1000);%for Criss values, ice
%      aDe_i=exp(C1D_i+C2D_i./TK + C3D_i./(TK.^2)); %new values from Mads Ellehoj,in KELVIN, vapor-solid
    aDe_i=exp(-0.0559 + ( 13525./(TK.^2)));%new value from Lamb et al, not yet published!!
    
% %BRM edit jan 24th 2022
% aDe_l=nan(size(TK));
% aDe_l(TK>=273)=exp((C1D_l+(10^3).*C2D_l./TK(TK>=273) + (10^6).*C3D_l./(TK(TK>=273).^2))./1000);%above zero
% aDe_l(TK<273)=exp((C1D_l0+(10^3).*C2D_l0./TK(TK<273) + (10^6).*C3D_l0./(TK(TK<273).^2))./1000);%below zero

      a18Oe_l=exp((C1O_l+(10^3).*C2O_l./TK + (10^6).*C3O_l./(TK.^2))./1000);
     a18Oe_i=exp((C1O_i+(10^3).*C2O_i./TK + (10^6).*C3O_i./(TK.^2))./1000);%for Criss values, ice
%      a18Oe_i=exp(C1O_i+C2O_i./TK + C3O_i./(TK.^2)); %new values from Mads Ellehoj,in KELVIN, vapor-solid

%      
     a17Oe_l=a18Oe_l.^0.529;
%      a17Oe_i=a18Oe_i.^0.529;
     a17Oe_i=a18Oe_i.^0.531;%maybe this is right?? Miller 2017, Precipitation regime influence on oxygen triple-isotope distributions in Antarctic precipitation and ice cores


%% Determine kinetic fractionation factor for condensation
% Jouzel Merlivat 78 gives D'/D (where ' denotes heavy) = 0.9723 for d18O and 0.9755 for dD vapor-liquid diffusion
%    DDdO = 1/0.9723; % = 1.0285 (or 28.5 permil more for H218O than H216O) at 20C for pure molecular diffusion, and 1.000 for pure turbulence
%    DDdD = 1/0.9755; % = 1.0251 at 20C  ***Note calculated on equal collision diameters is 1.0165

%!!! a_diff = (D/D')^n; they are not the same be, careful of everything, n=1 for no turbulence which is what Luz and Barkacn studied. You can use ph_diff for n~=1, since the n's don't matter too much. 
%below!!

%use luz and barkan data for diffusivity of vapor with no turbulence
%     DDdO = 1.0096; %From Luz & Barkan, 2010 pg 6282,
%     DDdO = 1.008; %from Uemura, this is wrt sea water, not appropriate the in air diffusivity
%     DDdD = 1.0164; %From Cappa 2003, dD @25C

% % JM 84
% DDdO = 1.0285;%From JM 84, ignoring ventilation effect. I *think* these are the values the used.
% DDdD = 1.0251;%From JM 84, ignoring ventilation effect. I *think* these


% are the values the used. Phi_diff = 0.88

%     DDdO = 1.028; %From Luz & Barkan, 2009 at 20C,  Fractionation of oxygen and hydrogen isotopes in evaporating water (pure molecular diffusion)
%     DDdD = 1.023; %From Luz & Barkan, 2009 at 20C Fractionation of oxygen and hydrogen isotopes in evaporating water (pure molecular diffusion)

    %      D_adiff = 1.023;   %D_adiff is same as DDdD ***no this is not true!***, notation from Luz & Barkan, 2009 Fractionation of oxygen and hydrogen isotopes in evaporating water
%      phi_diff = 0.88; %(D_adiff-1)/(DDdO-1); % Merlivat, 1978 had phi_diff of 0.88
%    DDdD = phi_diff *(DDdO -1)+1;%only ture if n=1!
  
%phi diff = (aDdiff -1)/(a18Odiff -1) = (DDdD^n -1)/(DDdO^n -1)
% phi_diffT = (1.25-0.02.*T); %Temperature dependence of phi_diff
% DDdD = phi_diffT *(DDdO -1)+1;%only ture if n=1!

% DDd17 = DDdO.^0.518;%does 0.518 hold under turbelent regimes?


% %~~~~~~comment/uncomment here~~~~~~
% %lets try a different approach
% % also, this notation sucks. I'm going to use D_O, D_17, and D_D, for
% % alpha_diff during transport. a_diff= DDd*^n, in the notation above
% %D_O can range from 1.0 (pure turbulense) to 1.0028 (pure molecular diffusion) (uemurea 2010 )
% 
% %pure molecular diffusion
% D_O = 1.0285;%From JM 84, ignoring ventilation effect. Pure moelcular diffusion
% % D_O = 1.018;
% % D_O = 1.017;
% % D_O = 1.007;
% % D_O = 1.025;
% % D_O=1.0;%pure turbulence
% 
% % D_D = 1.0251;%From JM 84, ignoring ventilation effect. This is the same as below.
% phi_diff = 0.88;%phi diff = (aDdiff -1)/(a18Odiff -1) = (DDdD^n -1)/(DDdO^n -1)
% D_D = phi_diff*(D_O -1)+1;
% % 
% D_17 = D_O^0.518; % from Barkan et al 2007
% 
% D_17 = D_O^0.514; % from  Robert Hellmann, Allan H. Harvey,
% % First?Principles Diffusivity Ratios for Kinetic Isotope Fractionation of
% % Water in Air 2020?
% 
% 
% %Jouzel and Merlivat gives ak = S/(asveq*D/D'*(S-1)+1) where ase = equilibrium fractionation coeeff with respect to the solid phase
% %Luz and Barkan, 2010 give adiff = (Dlight/Dheavy)^n (where n denotes turbulence, n=0 fully turbulent, n=1 no turbulence) 
% %~~~~~~comment/uncomment here~~~~~~


% BRM edit 4-19-2021
% Hellerman and Harvey 2020
% First‐Principles Diffusivity Ratios for Kinetic Isotope Fractionation of Water in Air
% https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020GL089999

T_star = TK./100;

D_r_HDO = 0.98258 - (0.02546./T_star) + (0.02421./(T_star.^(5/2)));

D_r_17 = 0.98284 + (0.003517./(T_star.^(1/2))) - (0.001996./(T_star.^(5/2)));

D_r_18 = 0.96671 + (0.007406./(T_star.^(1/2))) - (0.004861./(T_star.^(3)));


phi_diff_HH = (1 - D_r_HDO)./(1-D_r_18);

theta_diff_HH = log(D_r_17)./log(D_r_18);



D_O=1./D_r_18;
D_D=1./D_r_HDO;
D_17=1./D_r_17;
phi_diff = phi_diff_HH ;

%% initializaiton

% Initial Ocean Conditions
%smow values (actually vsmow) from wikipedia https://en.wikipedia.org/wiki/Vienna_Standard_Mean_Ocean_Water
%confirmed by Criss pg 21
R18O_smow=0.00200520;
RD_smow=0.00015576;
R17O_smow=0.0003799;


RD_v0 = (1+dD_v0./1000).*RD_smow;
R18O_v0 = (1+d18O_v0./1000).*R18O_smow;
d17O_v0=(exp(d17Oxs_v0+0.528.*log(1+d18O_v0./1000))-1).*1000;
R17O_v0 = ((d17O_v0./1000)+1).*R17O_smow;

%% DO I NEED TO INITIALIZE THESE????? TEST IT! brm march 2017
% alpha_tot = alpha_l*frac_l + alpha_i*frac_i
% alpha_i= alpha_i_eq*alpha_i_k; alpha_l= alpha_l_eq*alpha_l_k;

% S(1)=1;

%Set up initial kinetic frac factors for ice(_i) and liquid(_l).
    %The kinetic fractionation factors for liquid will always be 1 here. This
    %COULD be different if, for example, we considered re-evaporation of
    %precip. We will ignore this for now.

    aDk_l(1)=1; %Assume that condensation happens at 100% RH wrt liquid
    a18Ok_l(1)=1;
%     a17Ok_l(1)=a18Ok_l(1)^0.518;%this should be a18Ok_l^0.518, 
    a17Ok_l(1)=a18Ok_l(1)^0.514;%% from  Robert Hellmann, Allan H. Harvey,

    aDk_i(1)=1;%=S/(aDe_i*Diff_D*(S-1) +1); for efficiency start at =1. If we want to consider sources a T<T_crit, we'll have to change this.
%     aDk_i(1)=ss(1)./(aDe_i(1).*D_D(1).*(ss(1)-1)+1); %for efficiency start at =1. If we want to consider sources a T<T_crit, we'll have to change this.
    a18Ok_i(1)=1;
%     a17Ok_i(1)=a18Ok_i(1)^0.518;%this should be a18Ok_i^0.518, 
    a17Ok_i(1)=a18Ok_i(1)^0.514; % from  Robert Hellmann, Allan H. Harvey,
%     a17Ok_i(1)=a18Ok_i(1)^0.531;%this should be a18Ok_i^0.518, %maybe this is right?? Miller 2017, Precipitation regime influence on oxygen triple-isotope distributions in Antarctic precipitation and ice cores

aD(1)=(aDe_i(1)*aDk_i(1))*fraction_i(1) + (aDe_l(1)*aDk_l(1))*fraction_l(1);%for most starting cases this should = aDe_l(1)
a18O(1)=(a18Oe_i(1)*a18Ok_i(1))*fraction_i(1) + (a18Oe_l(1)*a18Ok_l(1))*fraction_l(1);
a17O(1)=(a17Oe_i(1)*a17Ok_i(1))*fraction_i(1) + (a17Oe_l(1)*a17Ok_l(1))*fraction_l(1);

%set up the cloud initial conidtions. Its just the initial vapor
R18O_c(1)=R18O_v0;        
R17O_c(1)=R17O_v0;
RD_c(1)=RD_v0;
%set up intial precip. This is in equilibrium with the cloud
R18O_p(1) = R18O_c(1) .*a18O(1);%RO precip = RO of cloud *alpha
d18O_p(1)= ((R18O_p(1)/R18O_smow) -1)*1000; %?18O of precip =(RO precip/R SMOW)-1)*1000
d18O_p_ln(1) = log(R18O_p(1)/R18O_smow)*1000; 

R17O_p(1) = R17O_c(1).*a17O(1);
d17O_p(1)= ((R17O_p(1)/R17O_smow) -1)*1000;
d17O_p_ln(1) = log(R17O_p(1)/R17O_smow)*1000;

RD_p(1) = RD_c(1) .*aD(1);%RD precip = RD of cloud *alpha
dD_p(1)= ((RD_p(1)/RD_smow) -1)*1000; %?D of precip =(RD precip/R SMOW)-1)*1000
dD_p_ln(1) = log(RD_p(1)/RD_smow)*1000;

dlnR18O(1)=0;
dlnRD(1)=0;
dlnR17O(1)=0;


%% run the numerical loop for Rayleigh Distilation with non-constant alphas and predetermined super saturation

%If we assume S is a function of temp. we can calculate all the fractionation factors now.
    % Let supersaturation be a simple function of T

S=ss;
% S=supersat;
%Alternate super saturation functions
% S_JM1 = 0.99 - 0.006.*(T);%EQ 15 JM84
% S_JM2 = 0.05 + 0.906 .*exp(-0.008.*T);  %EQ 16 JM84  
% S=S_JM1;

    %First D
    aDk_i = S./(aDe_i.*D_D.*(S-1)+1); %kinetic fractionation factor for ice-vapor
    aDk_l = (ss_l) ./ (aDe_i.*D_D.*((ss_l)-1)+1); %kinetic fractionation
%     factor for liquid-vapor. This is a little odd. We assume saturation
%     over liquid. But our PA path and ice/liquid frac imply sub
%     saaturation when there is very little liquid. Turns out it doesnt
%     matter at all. 
    aD = (aDe_i.*aDk_i).*fraction_i + (aDe_l.*aDk_l(1)).*fraction_l;%for aDk_l = 1, turn off above line
    %now 18O
    a18Ok_i = S./(a18Oe_i.*D_O.*(S-1)+1); %kinetic fractionation factor for ice-vapor
    a18Ok_l = (ss_l)./(a18Oe_i.*D_O.*((ss_l)-1)+1); %kinetic fractionation factor for ice-vapor
    a18O = (a18Oe_i.*a18Ok_i).*fraction_i + (a18Oe_l.*a18Ok_l(1)).*fraction_l;%for a18Ok_l = 1, turn off above line
    %now 17O
    a17Ok_i = S./(a17Oe_i.*D_17.*(S-1)+1); %kinetic fractionation factor for ice-vapor
    a17Ok_l = (ss_l)./(a17Oe_i.*D_17.*((ss_l)-1)+1); %kinetic fractionation factor for ice-vapor
    a17O = (a17Oe_i.*a17Ok_i).*fraction_i + (a17Oe_l.*a17Ok_l(1)).*fraction_l;%a17Ok_l = 1 , turn off above line.


% This for loop is for predetermined alphas, i.e. predetermined S
for tstep = 2:length(T)

%d ln(R) = (alpha - 1) d ln(f); this is the essential differential EQ of
%Rayleigh distillation for the vapor
dlnR18O(tstep)=(a18O(tstep)-1)*(log(f(tstep))-log(f(tstep-1)));%this is the change in vapor 
R18O_c(tstep)=exp(log(R18O_c(tstep-1))+dlnR18O(tstep));
if dlnR18O(tstep)~=0
R18O_p(tstep) = R18O_c(tstep) .*a18O(tstep);%RO precip = RO of cloud *alpha
else
    R18O_p(tstep)=NaN;%if there is no change in f then there is no precip
end
d18O_p(tstep)= ((R18O_p(tstep)/R18O_smow) -1)*1000; %delta 18O of precip =(RO precip/R SMOW)-1)*1000
d18O_p_ln(tstep) = log(R18O_p(tstep)/R18O_smow)*1000; 
    
dlnR17O(tstep)=(a17O(tstep)-1)*(log(f(tstep))-log(f(tstep-1)));%this is the change in vapor
R17O_c(tstep)=exp(log(R17O_c(tstep-1))+dlnR17O(tstep));
if dlnR17O(tstep)~=0
R17O_p(tstep) = R17O_c(tstep).*a17O(tstep);
else
    R17O_p(tstep)=NaN;%if there is no change in f then there is no precip
end
d17O_p(tstep)= ((R17O_p(tstep)/R17O_smow)-1)*1000;
d17O_p_ln(tstep) = log(R17O_p(tstep)/R17O_smow)*1000;

    
dlnRD(tstep)=(aD(tstep)-1)*(log(f(tstep))-log(f(tstep-1)));%this is the change in vapor
RD_c(tstep)=exp(log(RD_c(tstep-1))+dlnRD(tstep));
if dlnRD(tstep)~=0
RD_p(tstep) = RD_c(tstep) .*aD(tstep);%RD precip = RD of cloud *alpha
else
    RD_p(tstep)=NaN;%if there is no change in f then there is no precip
end
dD_p(tstep)= ((RD_p(tstep)/RD_smow)-1)*1000; %delta D of precip =(RD precip/R SMOW)-1)*1000
dD_p_ln(tstep) = log(RD_p(tstep)/RD_smow)*1000;

end
dxs= dD_p -(8.*d18O_p);
%log definition of dxs from Uemura Clim Past 2012
%ln(1 + dD) = -2.85 * 10^-2 * (ln(1 + d18O))^2 +8.47 * ln(1 + d18O) +13.3;%offset isn't needed except for modeling plots
dxs_ln = (dD_p_ln-(-2.85.*10^-2.*(d18O_p_ln).^2+8.47.*(d18O_p_ln)));%+13.3;
d17Oxs = (log(1+d17O_p./1000)-0.528*log(1+d18O_p./1000))*1e6;%units are weird     