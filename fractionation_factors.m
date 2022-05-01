function [aD, a18O, a17O, aDe_l, aDe_i, a18Oe_l, a18Oe_i, a17Oe_l,  a17Oe_i, aDk_i, a18Ok_i, a17Ok_i] = fractionation_factors(T,ss,fraction_i,fraction_l)
%cacluates fractionation factors for major isotopolgues based on
%temperature
%[aD, a18O, a17O, aDe_l, aDe_i, a18Oe_l, a18Oe_i, a17Oe_l,  a17Oe_i, aDk_i, a18Ok_i, a17Ok_i] = fractionation_factors(T,ss,fraction_i,fraction_l)

TK = T+273.15;% converts to deg K


%% Equilibrium fractionation factors
% a=exp((C1+(10^3)*C2/T +(10^6)*C3/(T^2))/1000);
%Table 3.2, Criss, Chapter 3, calculation of Fractionation Factor, a, for
%various circumstances, D=deuterium, O=oxygen, l=liquid, v=vapor, i=ice
%Type, phase A, phase B, C1, C2, C3, T min (K), T max (K)
% 'D' 'l' 'v' 52.612 -76.248 24.844 273 373;
% 'D' 'l' 'v' -1000.0 0 15.013 258 273;
% 'D' 'i' 'v' -94.5 0 16.289 233 273;
% 'O' 'l' 'v' -2.0667 -0.4156 1.137 273 373;
% 'O' 'i' 'v' -28.224 11.839 0 240 273;
%Note: this table uses 273 K (not 273.15 K); we'll assume these are the
%same!

%X is Table 3.2 without "Type, Phase A, and Phase B"
X=[52.612 -76.248 24.844 273 373; -100.0 0 15.013 258 273; -94.5 0 16.289 233 273;...
    -2.0667 -0.4156 1.137 273 373; -28.224 11.839 0 240 273];

%Z is new equilibrium fractionation factor for -40 to 0C from Mads Ellehoj
Z=[0.2133 -203.10 48888 233 273; 0.0831 -49.192 8312.5 233 273];


%OK now we have to account for the fraction of ice and water
    C1D_l=X(1,1); %For D and O above T>0, from Criss
    C2D_l=X(1,2); 
    C3D_l=X(1,3);
    C1O_l=X(4,1);
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
    
     a18Oe_l=exp((C1O_l+(10^3).*C2O_l./TK + (10^6).*C3O_l./(TK.^2))./1000);
     a18Oe_i=exp((C1O_i+(10^3).*C2O_i./TK + (10^6).*C3O_i./(TK.^2))./1000);%for Criss values, ice
%      a18Oe_i=exp(C1O_i+C2O_i./TK + C3O_i./(TK.^2)); %new values from Mads Ellehoj,in KELVIN, vapor-solid

%      
     a17Oe_l=a18Oe_l.^0.529;
     a17Oe_i=a18Oe_i.^0.529;

%% Determine kinetic fractionation factor for condensation
% Jouzel Merlivat 78 gives D'/D (where ' denotes heavy) = 0.9723 for d18O and 0.9755 for dD vapor-liquid diffusion
%    DDdO = 1/0.9723; % = 1.0285 (or 28.5 permil more for H218O than H216O) at 20C for pure molecular diffusion, and 1.000 for pure turbulence
%    DDdD = 1/0.9755; % = 1.0251 at 20C  ***Note calculated on equal collision diameters is 1.0165

%!!! a_diff = (D/D')^n; they are not the same be careful of everything, n=1 for no turbulence which is what Luz and Barkacn studied. You can use ph_diff for n~=1, since the n's don't matter too much. 
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



%lets try a different approach
% also, this notation sucks. I'm going to use D_O, D_17, and D_D, for
% alpha_diff during transport. a_diff= DDd*^n, in the notation above
%D_O can range from 1.0 (pure turbulense) to 1.0028 (pure molecular diffusion) (uemurea 2010 )

%pure molecular diffusion
D_O = 1.0285;%From JM 84, ignoring ventilation effect. Pure moelcular diffusion
% D_O = 1.017;
% D_O = 1.007;
% D_O = 1.01;
% D_O=1.0;%pure turbulence

% D_D = 1.0251;%From JM 84, ignoring ventilation effect. This is the same as below.
phi_diff = 0.88;%phi diff = (aDdiff -1)/(a18Odiff -1) = (DDdD^n -1)/(DDdO^n -1)
D_D = phi_diff*(D_O -1)+1;

D_17 = D_O^0.518; % from Barkan et al 2007


%%
%If we assume S is a function of temp. we can calculate all the fractionation factors now.
    % Let supersaturation be a simple function of T

S=ss;
% S=supersat;
%Alternate super saturation functions
% S_JM1 = 0.99 - 0.006.*(T);%EQ 15 JM84
% S_JM2 = 0.05 + 0.906 .*exp(-0.008.*T);  %EQ 16 JM84  
% S=S_JM1;

% for now ak_l=1;
aDk_l = 1;
a18Ok_l = 1;
a17Ok_l = 1;


    %First D
    aDk_i = S./(aDe_i.*D_D.*(S-1)+1); %kinetic fractionation factor for ice-vapor
    aD = (aDe_i.*aDk_i).*fraction_i + (aDe_l.*aDk_l(1)).*fraction_l;%aDk_l = 1 for now.
    %now 18O
    a18Ok_i = S./(a18Oe_i.*D_O.*(S-1)+1); %kinetic fractionation factor for ice-vapor
    a18O = (a18Oe_i.*a18Ok_i).*fraction_i + (a18Oe_l.*a18Ok_l(1)).*fraction_l;%a18Ok_l = 1 for now.
    %now 17O
    a17Ok_i = S./(a17Oe_i.*D_17.*(S-1)+1); %kinetic fractionation factor for ice-vapor
    a17O = (a17Oe_i.*a17Ok_i).*fraction_i + (a17Oe_l.*a17Ok_l(1)).*fraction_l;%a18Ok_l = 1 for now.








%Jouzel and Merlivat gives ak = S/(asveq*D/D'*(S-1)+1) where ase = equilibrium fractionation coeeff with respect to the solid phase
%Luz and Barkan, 2010 give adiff = (Dlight/Dheavy)^n (where n denotes turbulence, n=0 fully turbulent, n=1 no turbulence) 

