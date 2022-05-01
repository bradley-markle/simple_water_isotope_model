function [dD_v0, d18O_v0, d17Oxs_v0,RHn0,rh0,sst0] =evaporation_2022(T0,SST0,RH0,closure,reanalysis, SH, season, sw, alpha, unc)
% BRM March 2017: 
% [dD_v0, d18O_v0, d17Oxs_v0,RHn0,rh0,sst0] = evaporation_model(T0,SST0,RH0,closure,reanalysis, SH);
% evaporation_model is a function that solves for the delta values of inital vapor
% evaporated from the ocean. Specify the initial surface air temperaure, T0, in degrees C.
% The user can specify the inital sea surface temperaure (in degrees C), SST0, and initial
% Relative Humidity, RH0, as a number from 0 to 1.
% The user may specify which closure assumption you wish to use. Enter
% either 'local' or 'global'. The default is 'local'.
% To run model in climatiology mode, leave SST0 =[], and RH0 =[], or just don't enter them.
% In this set up the model will calculate SST0 and RH0 using best fits from
% reanalysis climatology given the surface temperature, T0, and the
% hemisphere of choice (discussed below). 
% The user may specifiy which reanalysis climatology to use, either 'ncep'
% or 'era'. The default is 'ncep'. The climatologies are very similar, but
% ERA has a slightly different normalized relative humidity since the skin
% temperature is available as a different field than the surface air
% temperature.
% To specify the Southern Hemisphere, enter SH = 1. To sepcify
% the Northern Hemisphere enter SH = 0; The default is SH = 1.
% There are many, many different assumptions and parameters that can be
% tested and changed within the code.
%% House Keeping
%set up defaults
% if nargin == 3
%     closure = 'local';
%     reanalysis = 'ncep';
%     SH = 1;
% end
% if nargin == 4
%     reanalysis = 'ncep';
%     SH = 1;
% end
% if nargin < 6
%     SH = 1;
% end
if ~exist('SST0','var')
    SST0 = [];%default is to use correlations
end
if ~exist('RH0','var')
    RH0 = [];%default is to use correlations
end
if ~exist('closure','var')
    closure = 'local';%default closure assumption is local
end

if ~exist('reanalysis','var')
    reanalysis = 'ncep';%default reanalysis for correlation is NCEP
end

if ~exist('SH','var')%default hemisphere is Southern
    SH = 1;
end

if ~exist('season','var')
    season = 'annual';%default seasonality is annual average
end

if strcmp(reanalysis,'ncep')
    season = 'annual'; %the NCEP correlations only support annual averages
end
%% Set up initial conditions
T0K=T0+273.15;% converts to deg K

P0 = 101.325;%kPa, initial pressure
% P0=98;
%% calculate SST, RH, RHn. If not specified use reanalysis correlations

if isempty(SST0) && isempty(RH0)
%     if strcmp(reanalysis,'ncep')
%     [rh0, delta_rh0, sst0, delta_sst0, RHn0, deltarhn] = T_RH_RHn_NCEP(T0, SH);%use NCEP data for correlations
% %     [rh0, delta_rh0, sst0, delta_sst0, RHn0, deltarhn] = T_RH_RHn_NCEP_alt(T0, SH,0);%use NCEP data for correlations
% %     offset=5;%offset the realtive humidity
% %     [rh0, delta_rh0, sst0, delta_sst0, RHn0, deltarhn] = T_RH_RHn_NCEP_offset(T0, SH, offset);
% 
%     elseif strcmp(reanalysis,'era')
%     [rh0, delta_rh0, sst0, delta_sst0, RHn0, deltarhn0, T_skin0, delta_T_skin0] = T_RH_RHn_ERA(T0, SH,season);%use ERA data for correlations
% %         [rh0, delta_rh0, sst0, delta_sst0, RHn0, deltarhn0] = T_RH_RHn_ERA_alt(T0, SH,season,0);%use ERA data for correlations
% %         offset=5;
% %         [rh0, delta_rh0, sst0, delta_sst0, RHn0, deltarhn0] = T_RH_RHn_ERA_offset(T0, SH,season,offset);%use ERA data for correlations
% 
%     end
[rh0, delta_rh0, sst0, delta_sst0, RHn0, deltarhn0] = T_RH_RHn_2022(T0, SH, reanalysis, 'spline', season, 0);

    rh0=rh0/100;
end

if ~isempty(SST0) && isempty(RH0)
    sst0 = SST0;
    sst0K=sst0+273.15;% converts to deg K   
%calculate RH
%     if strcmp(reanalysis,'ncep')
%     [rh0, delta_rh0, ~, ~, RHn0, deltarhn] = T_RH_RHn_NCEP(T0, SH);%use NCEP data for correlations
%     elseif strcmp(reanalysis,'era')
%     [rh0, delta_rh0, ~, ~, RHn0, deltarhn0, T_skin0, delta_T_skin0] = T_RH_RHn_ERA(T0, SH,season);%use ERA data for correlations
%     end
    [rh0, delta_rh0, ~, ~, RHn0, deltarhn0] = T_RH_RHn_2022(T0, SH, reanalysis, 'spline', season, 0);

    rh0=rh0/100;
end

if isempty(SST0) && ~isempty(RH0)
%calculate sst
%     if strcmp(reanalysis,'ncep')
%     [~, ~, sst0, delta_sst0, ~, ~] = T_RH_RHn_NCEP(T0, SH);%use NCEP data for correlations
%     elseif strcmp(reanalysis,'era')
%     [~, ~, sst0, delta_sst0, ~, ~, T_skin0, delta_T_skin0] = T_RH_RHn_ERA(T0, SH,season);%use ERA data for correlations
%     end
    [~, ~, sst0, delta_sst0,~, ~] = T_RH_RHn_2022(T0, SH, reanalysis, 'spline', season, 0);

    
    sst0K=sst0+273.15;% converts to K   
%calculate RHn
rh0=RH0;
    e_s_skin0 = (1000^-1).*exp(54.842763 - 6763.22 ./ T0K - 4.21 .* log(T0K) + 0.000367 .* T0K +...
    tanh(0.0415 .* (T0K - 218.8)) .*  (53.878 - 1331.22 ./T0K - 9.44523 .* log(T0K) + 0.014025 .* T0K)) ;%with T in [K] and ew in [kPa]
    e_s_SST0 = (1000^-1).*exp(54.842763 - 6763.22 ./ sst0K - 4.21 .* log(sst0K) + 0.000367 .* sst0K +...
    tanh(0.0415 .* (sst0K - 218.8)) .*  (53.878 - 1331.22 ./ sst0K - 9.44523 .* log(sst0K) + 0.014025 .* sst0K)) ;%with T in [K] and ew in [kPa
RHn0=real((rh0.*e_s_skin0)./e_s_SST0);%normalized realitive humidity, using T0 as skin temp
end

if ~isempty(SST0) && ~isempty(RH0)
    sst0 = SST0;
    sst0K=sst0+273.15;% converts to deg K
%calculate RHn
rh0=RH0;
    e_s_skin0 = (1000^-1).*exp(54.842763 - 6763.22 ./ T0K - 4.21 .* log(T0K) + 0.000367 .* T0K +...
    tanh(0.0415 .* (T0K - 218.8)) .*  (53.878 - 1331.22 ./T0K - 9.44523 .* log(T0K) + 0.014025 .* T0K)) ;%with T in [K] and ew in [kPa]
    e_s_SST0 = (1000^-1).*exp(54.842763 - 6763.22 ./ sst0K - 4.21 .* log(sst0K) + 0.000367 .* sst0K +...
    tanh(0.0415 .* (sst0K - 218.8)) .*  (53.878 - 1331.22 ./ sst0K - 9.44523 .* log(sst0K) + 0.014025 .* sst0K)) ;%with T in [K] and ew in [kPa
RHn0=real((rh0.*e_s_skin0)./e_s_SST0);%normalized realitive humidity, using T0 as skin temp
end

%this is for testing the sensitivity to the reanalysis correlations
%     sst0=sst0-delta_sst0;%for testing uncertainty
%     RHn=RHn+deltaRHn;%for testing uncertainty
%     rh0=rh0-delta_rh0;%for testing uncertainty

if ~exist('unc','var')
    unc= 0;%
end
if unc ==1
    sst0=sst0+delta_sst0;
    RHn0=RHn0;
    rh0=rh0; 
elseif unc ==2
    sst0=sst0-delta_sst0;
    RHn0=RHn0;
    rh0=rh0; 
elseif unc ==3
    sst0=sst0;
    RHn0=RHn0+deltarhn0;
    rh0=rh0+delta_rh0;
elseif unc ==4
    sst0=sst0;
    RHn0=RHn0-deltarhn0;
    rh0=rh0-delta_rh0;
end


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
    
%Z is new equilibrium fractionation factor for -40 to 0C from Mads Ellehoj,
%I thnk it is wrong. Use newer factors from Lamb et al.
% Z=[0.2133 -203.10 48888 233 273; 0.0831 -49.192 8312.5 233 273];
%     C1D_i=Z(1,1); %For D and O below T<=0 using Mads Ellehoj Z matrix above
%     C2D_i=Z(1,2); 
%     C3D_i=Z(1,3);
%     C1O_i=Z(2,1);
%     C2O_i=Z(2,2); 
%     C3O_i=Z(2,3);  

% Calculate fractionaation factors for D, 18O, 17O for ice and liquid

     aDe_l=exp((C1D_l+(10^3).*C2D_l./T0K + (10^6).*C3D_l./(T0K.^2))./1000);% Liquid fractionation factor 
%      aDe_i=exp((C1D_i+(10^3).*C2D_i./T0K + (10^6).*C3D_i./(T0K.^2))./1000);%for Criss values, ice
%      aDe_i=exp(C1D_i+C2D_i./T0K + C3D_i./(T0K.^2)); %new values from Mads Ellehoj,in KELVIN, vapor-solid
    aDe_i=exp(-0.0559 + ( 13525./(T0K.^2)));%new value from Lamb et al, not yet published!!
    
     a18Oe_l=exp((C1O_l+(10^3).*C2O_l./T0K + (10^6).*C3O_l./(T0K.^2))./1000);% Liquid fractionation factor
     a18Oe_i=exp((C1O_i+(10^3).*C2O_i./T0K + (10^6).*C3O_i./(T0K.^2))./1000);%for Criss values, ice
%      a18Oe_i=exp(C1O_i+C2O_i./T0K + C3O_i./(T0K.^2)); %new values from Mads Ellehoj,in KELVIN, vapor-solid
      
     a17Oe_l=a18Oe_l.^0.529;
     a17Oe_i=a18Oe_i.^0.529;
%% Initial evaporation

%For evaporation we're gonna use diffusivities
 
%We're going use a_evap valued from Luz et al rather than Criss.
%However, we're not gonna use their a_diff values, because they were for
%non turbulent flow. We'll use Uemura's values from over the ocean for
%that.

% Initial Ocean Conditions
%smow values (actually vsmow) from wikipedia https://en.wikipedia.org/wiki/Vienna_Standard_Mean_Ocean_Water
%confirmed by Criss pg 21
R18O_smow=0.00200520;
RD_smow=0.00015576;
R17O_smow=0.0003799;

% however Southern ocean ocean water is not the same as mean ocean water,
% it may have a delta value, or it may change in time.

if ~exist('sw','var')||isempty(sw)
    sw = -0.3;%[per mil]
end
d18O_ocean=sw;
% d18O_ocean=0;
% d18O_ocean=0.7;%LGM like conditions
% d18O_ocean=-0.3; %delta 18O value of Southern Ocean, from http://data.giss.nasa.gov/o18data/
% d18O_ocean=0.3; %delta 18O value more like low latitudes, from http://data.giss.nasa.gov/o18data/

% make variable d18O_ocean based on NASA GISS data in which T0<10
% d18O_sw~=-0.2 and T0>10 d18O_sw~=+0.2 
% d18O_ocean=zeros(size(T0));
% d18O_ocean(T0<10)=-0.3;
% d18O_ocean(T0>10)=0.3;


% d18O_ocean=-1;
Rwi_18O=(1+d18O_ocean./1000).*R18O_smow;%initial R value of ocean water

% dD_ocean=0; %delta D value of Ocean
% dD_ocean=8*d18O_ocean-12; %delta D value of Ocean
% dD_ocean=-1;%just for sensitivity testing 
% dD_ocean=8*d18O_ocean; %assumption following Uemura 2012
%   %alternate assumption. Use Uemura 2012 relationship for d18Oln of SW to
%   %dDln of SW. convert from standard d18O. Doesn't make too big of a
%   %difference.
% dD_ln_ocean=-2.85 * 10^-2 * (log(1 + d18O_ocean./1000).*1000).^2+ 8.47 * log(1 + d18O_ocean./1000).*1000;
% dD_ocean=(exp(dD_ln_ocean/1000)-1)*1000;

dD_ocean=d18Osw_to_dDsw(d18O_ocean);%uses slope of available data

Rwi_D=(1+dD_ocean./1000).*RD_smow;

%initial d17O of ocean 
% d17O_ocean=0;
d17O_ocean= ((exp(0+0.529*(log(d18O_ocean/1000+1))))-1)*1000;
Rwi_17O=(1+d17O_ocean/1000)*R17O_smow;


%% Initial fractionation

% Calculate initial vapor from ocean, using relative humidity
% Mostly Following Criss

% ~~~~~~~~~~~~~~This is all from Criss
%this is table 4.2 from Criss, for a_evap0, calculated from eq 4.67, for
%n=0.5, 
%temps 0, 10, 20, 100 degC
crissTable4_2=[1.122 1.0280; 1.107 1.0270; 1.094 1.0261; 1.0355 1.0212];
%!! note that Luz and Barkan 2009 have different alpha evap values. maybe
%check out.

%% First find a_evap
%~~~~~~~~~~~~~~~~~~~~~~~~
%   Option 1) follow Criss
% a_evap=Rw/Re=(a_evap0 (1-h)) /( 1- a_eq h Rv/Rw)

%%---turn on/off here
a18O_evap=1.0045;%from Criss pg 162-163
aD_evap=1.0267;%from Criss pg 162-163
% a17O_evap=a18O_evap.^0.518;
%%---turn on/off here

%~~~~~~~~~~~~~~~~~~~~~~~~
% %   Option 2) use values from Luz, Barkan et al 2009

% % DO NOT Use these. These are under VERY VERY low humidity. 

% %I made the following peicewise linear fits to the a_evap data in their
% %Table 2
% 
% %%---turn on/off here
% if sst0<20.1
% aD_evap = (1.1478-0.0019.*sst0);
% elseif sst0>=20.1
% aD_evap = (1.1331-0.0012.*sst0);
% end
% 
% if sst0<20.1
% a18O_evap = (1.0382-(6.5347e-05)*sst0);
% elseif sst0>=20.1
% a18O_evap = (1.0387-(8.6294-05)*sst0);
% end
% %%---turn on/off here

%% Next find a_evap0 or a_eq*a_diff

%~~~~~~~~~~~~~~~~~~~~~~~~
%   Option 1) follow Criss
%Criss Values!
%values from Criss

% %%---turn on/off here
% aD_evap0=interpPH([-30 0 10 20 100]',[crissTable4_2(1,1); crissTable4_2(:,1)],T0); %from Criss Table 4.2
% a18O_evap0=interpPH([-30 0 10 20 100]',[crissTable4_2(1,2); crissTable4_2(:,2)],T0);%from Criss Table 4.2
% 
% %this uses above values and Uemura's formula for 17O
% a18O_diff=a18O_evap0./a18Oe_l;
% a17O_diff=a18O_diff.^0.518;% Uemura 2010
% a17O_evap0=a17O_diff.*a17Oe_l;
% %%---turn on/off here

%~~~~~~~~~~~~~~~~~~~~~~~~
%   Option 2) follow Uemura values for adiff, calculate a_evap0, use Luz
%   and Barkan phi_diff

%alternate method for a_evap0, using a_diff from Barkan and Luz, and Uemura,  % a_evap0 = a_diff*a_eq in Barkan and Luz 2007
%try values form Uemura
% I think adiff values from uemura are a_evap0? No! a_evap0 = a_diff*a_eq in Barkan and Luz 2007
% % a18O_diff=1.008;%Uemura 2010; this is Uemura's mean value. However the range seems to be from 1.005 to 1.011; stated uncertainty is 1.0083 ± 0.0018
% a18O_diff=1.007; %also from Uemura 2010 but optimized for dxs instead of o17xs.
% a17O_diff=a18Odiff^0.518;


if ~exist('alpha','var')||isempty(alpha)
    alpha = 1.009;
end
a18O_diff=alpha;


% a18O_diff=1.006;
% a18O_diff=1.007;%from Uemura 2010 optimized for dxs
% a18O_diff=1.008;%from Uemura 2010 optimized for d17Oxs
% a18O_diff=1.0076;%Pfahl and Wernli [2009] 
% a18O_diff=1.009;
% a18O_diff=1.010;
% a18O_diff=1.011;


%new tests June 2019. Make a_diff T dependant
% a18O_diff=nan(size(T0));
% a18O_diff(T0>15)=1.005;
% a18O_diff(T0<=15)=1.009;
% a18O_diff(T0<=0)=1.020;

%make it linearly dependant on temp.
% a18O_diff=nan(size(T0));
% y1=1.011;
% y2=1.006;
% x1=0;
% x2=30;
% m=(y2-y1)/(x2-x1);
% b=y1-(m*x1);
% a18O_diff=m*T0 + b;

%make a function that transitions betwee 1.006 and 1.009
% a18O_diff=nan(size(T0));
% a18O_diff=1.006 + (0.0015.*(1-erf((T0-15)/(30/4))));
% a18O_diff=1.007 + (0.0010.*(1-erf((T0-7)/(30/4))));
% a18O_diff=1.006 + (0.002.*(1-erf((T0-7)/(30/4))));
% a18O_diff=1.007 + (0.001.*(1-erf((T0-10)/(30/4))));
% a18O_diff=1.006 + (0.0014.*(1-erf((T0-25)/(30/4))));
% a18O_diff=1.006 + (0.0015.*(1-erf((T0-20)/(30/4))));
% a18O_diff=1.006 + (0.0015.*(1-erf((T0-17)/(30/4))));
% a18O_diff=1.0055 + (0.00175.*(1-erf((T0-17)/(30/4))));%pretty good

% a18O_diff=1.006 + (0.0015.*(1-erf((T0-15)/(30/4))))...
%     +(-0.0015.*(1-erf((T0-20)/(30/4))));
% a18O_diff=1.007 + (0.0015.*(1-erf((T0-7)/(30/4)))) + (-0.0015.*(1-erf((T0-15)/(30/4))));

% a18O_diff=1.009 + (0.002.*(1-erf((T0-20)/(30/4)))) + (-0.002.*(1-erf((T0-30)/(30/4))));

%next calculate aD_diff, use phi_diff from Luz and Barkan, et al 2009

% phi_diff
% (Merlivat, 1978) suggest phi_diff = 0.88, for evaporating conditions 
% phi_diff_sst = (1.25-0.02.*sst0);% from Luz, Barkan, et al 2009, this is what they suggest for phi_diff at cold temps!
% phi_diff_sst = (1.25-0.02.*T);% from Luz, Barkan, et al 2009, this is what they suggest for phi_diff at cold temps!

% %%%~~~!!! turn on/off here!!!~~~~~
% %I made my own ph_diff as a temperature dependant quantity from the Luz,
% %Barkan et al 2009 data for their obs at 10, 20.1, and 39.8 deg C
% phi_diff_T0=nan(size(T0));
% phi_diff_sst=nan(size(T0));
% 
% phi_diff_T0(T0<20.1) = (1.2778-0.0218.*T0(T0<20.1));
% phi_diff_T0(T0>=20.1) = (0.9318-0.0046.*T0(T0>=20.1));
% 
% phi_diff_sst(sst0<20.1) = (1.2778-0.0218.*sst0(sst0<20.1));
% phi_diff_sst(sst0>=20.1) = (0.9318-0.0046.*sst0(sst0>=20.1));
% 
% aD_diff = phi_diff_sst .*(a18O_diff -1)+1;%from Luz, Barkan et al 2009
% % aD_diff = phi_diff_T0 .*(a18O_diff -1)+1;%from Luz, Barkan et al 2009
% % aD_diff = 0.88.*(a18O_diff -1)+1;%Merlivat, 1978
% 
% 
% a18O_evap0 = a18O_diff.*a18Oe_l;
% aD_evap0 = aD_diff.*aDe_l;
% 
% 
% a17O_diff=a18O_diff.^0.518;% Uemura 2010
% a17O_evap0=a17O_diff.*a17Oe_l;
% %%%~~~!!! turn on/off here!!!~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~

%   Option 3) a_diff use values from Luz, and Barkan, 2009

% These may not be great since they are under labratory conditions. But
% whatever.
% 
% %%---turn on/off here
% if sst0<20.1
% aD_diff = (1.0344-(5.6832e-04)*sst0);
% elseif sst0>=20.1
% aD_diff = (1.0255-(1.2284e-04)*sst0);
% end
% 
% if sst0<20.1
% a18O_diff = (1.0267+(4.0594e-05)*sst0);
% elseif sst0>=20.1
% a18O_diff = (1.0276-(1.5228e-06)*sst0);
% end
% 
% 
% 
% a18O_evap0 = a18O_diff*a18Oe_l;
% aD_evap0 = aD_diff*aDe_l;
% 
% % a18O_evap0=1.007;
% % aD_evap0= phi_diff_sst *(a18O_evap0 -1)+1;
% % a17O_evap0=a18O_evap0.^0.518;% this is wrong!
% a17O_diff=a18O_diff.^0.518;% Uemura 2010
% a17O_evap0=a17O_diff*a17Oe_l;
% %%---turn on/off here


%   Option 4) results from Hellerman and Harvey 2020
% First‐Principles Diffusivity Ratios for Kinetic Isotope Fractionation of Water in Air
% https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020GL089999


T_star = T0K./100;

D_r_HDO = 0.98258 - (0.02546./T_star) + (0.02421./(T_star.^(5/2)));

D_r_17 = 0.98284 + (0.003517./(T_star.^(1/2))) - (0.001996./(T_star.^(5/2)));

D_r_18 = 0.96671 + (0.007406./(T_star.^(1/2))) - (0.004861./(T_star.^(3)));


phi_diff_HH = (1 - D_r_HDO)./(1-D_r_18);

theta_diff_HH = log(D_r_17)./log(D_r_18);

% % Option 4a) use a18_diff from above, calcualte the rest.
% aD_diff = phi_diff_HH  .*(a18O_diff -1)+1;
% a18O_evap0 = a18O_diff.*a18Oe_l;
% aD_evap0 = aD_diff.*aDe_l;
% 
% 
% a17O_diff=a18O_diff.^0.518;% Uemura 2010
% a17O_evap0=a17O_diff.*a17Oe_l;

% %Option 4b) calculate alll following HH
N=0.27;%equivalent-ish to a18O_diff=1.009
% N=0.304;
% N=0.302;
% N=0.32;

a18O_diff=(1./D_r_18).^N;
% aD_diff = phi_diff_HH  .*(a18O_diff -1)+1;% ah, this does this slightly
% wrong
aD_diff=(1./D_r_HDO).^N;

a18O_evap0 = a18O_diff.*a18Oe_l;
aD_evap0 = aD_diff.*aDe_l;

% a17O_diff=a18O_diff.^0.518;% Uemura 2010
a17O_diff=a18O_diff.^(log(1./D_r_17)./log(1./D_r_18));%using HH2020, correctly.... I think!
a17O_evap0=a17O_diff.*a17Oe_l;



%% Ok now we do the evaporation
if strcmp(closure,'global')
%now we will use a manipulation of Eg 4.31, Criss, pg 153
%Rv=Rw * (1- a_evap0*(1-h)/a_evap)/(a_eq*h)
% %~~~~~~~~~~~~~~~~
% 
xxx=(aD_evap0.*(1-RHn0))./aD_evap;
yyy=(aDe_l.*RHn0);
% yyy=(0.9989.*aDe_l*RHn);%factor to take into accoun molaity of sea water eq 4.40c Criss pg 160, M=0.550 mol kg^-1
zzz=(1-xxx)./yyy;
RD_v0=Rwi_D.*zzz;
dD_v0= ((RD_v0./RD_smow) -1).*1000;
clear xxx yyy zzz

xxx=(a18O_evap0.*(1-RHn0))./a18O_evap;
yyy=(a18Oe_l.*RHn0);
% yyy=(0.9989.*a18Oe_l*RHn);%factor to take into accoun molaity of sea water eq 4.40c Criss pg 160, M=0.550 mol kg^-1, oh shit this seems to make a difference
zzz=(1-xxx)./yyy;
R18O_v0=Rwi_18O.*zzz;
d18O_v0= ((R18O_v0./R18O_smow) -1).*1000;
clear xxx yyy zzz

% xxx=(a17O_evap0*(1-RHn))/a17O_evap;
% yyy=(a17Oe_l*RHn);
% % yyy=(0.9989.*a17Oe_l*RHn);%factor to take into accoun molaity of sea water eq 4.40c Criss pg 160, M=0.550 mol kg^-1,
% zzz=(1-xxx)/yyy;
% R17O_v0=Rwi_17O*zzz;
% d17O_v0= ((R17O_v0/R17O_smow) -1)*1000;
% d17Oxs_v0 = (log(1+d17O_v0./1000)-0.528*log(1+d18O_v0./1000));%*1e6;%units are weird
% clear xxx yyy zzz
end



%~~~~~~~~~~~~~~~~
if strcmp(closure,'local')
%alright lets up date this. We'll Equation 1) from Risi et al 2010
%Supplement, which is from Craig and Gordon 1965
% using the closure assumption, leads to EQ 3:
%R_bl = R_oce / (a_eq * (a_k +RHn * (1-a_k)))
%in my vars:
% R_v0 = Rwi / (ae_l * (a_diff + RHN * (1-a_diff)))

% %~~~~~~

%for d18O
R18O_v0 = Rwi_18O ./ (a18Oe_l .* (a18O_diff + RHn0 .* (1-a18O_diff)));%ok good news. this is the same as from criss
d18O_v0= ((R18O_v0./R18O_smow) -1).*1000;

%for dD
RD_v0 = Rwi_D ./ (aDe_l .* (aD_diff + RHn0 .* (1-aD_diff)));%ok this one is different than criss...
dD_v0= ((RD_v0./RD_smow) -1).*1000;%delta value ends up being identical however.
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  

%the above for d17O xs does not agree with Eq 4 in Landais et al 2008: Record of d18Oand17O-excess in ice from Vostok Antarctica duringthe last 150,000 years
% they use: -log((a18Oe_l^(0.529))*((a18O_diff^(0.518))*(1-RHn)+RHn)) +0.528*log((a18Oe_l)*((a18O_diff)*(1-RHn)+RHn))
% the problem may be a17O_evap. I don't think I did that right.

% %check this!! I dont think this uses R17O_ocean... BRM 11/20/16
% d17Oxs_v0=-log((a18Oe_l^0.529).*((a18O_diff.^0.518)*(1-RHn)+RHn))+0.528.*log((a18Oe_l).*((a18O_diff)*(1-RHn)+RHn));
% d17O_v0= ((exp(d17Oxs_v0+0.528*(log(d18O_v0/1000+1))))-1)*1000;
% R17O_v0=1+d17O_v0/1000;

%use eq 4 from landais et al 2008/ eq 2 from uemura et al 2010
d17Oxs_v0=-log((a18Oe_l.^0.529).*((a18O_diff.^0.518).*(1-RHn0)+RHn0))+0.528.*log((a18Oe_l).*((a18O_diff).*(1-RHn0)+RHn0));
d17O_v0=(exp(d17Oxs_v0+0.528.*log(1+d18O_v0./1000))-1).*1000;
R17O_v0=((d17O_v0./1000)+1).*R17O_smow;

end
