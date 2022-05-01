function [T_site, T_source, RH_source, d18O_site, dD_site, d18Oln_site, dDln_site, dxs_site, d17O_xs_site, dlnU_site, r_s_site, P_site] = simple_water_isotope_model_2019(Tsite, Tsource, RHsource, a, b, c, closure, reanalysis, SH, season)
% BRM March 2017
% [T_site, T_source, RH_source, d18O_site, dD_site, d18Oln_site, dDln_site, dxs_site, d17O_xs_site, dlnU_site] = simple_water_istope_model(Tsite, Tsource, RHsource, a, b, c, closure, reanalysis, SH, season);
% 
% Enter T_site as a three element vector of the lowest site (condensation) temperature, the highest site temperature, and the step, in degrees C: Tsite = [Tsite_lo Tsite_hi dTsite];
% Enter T_source as a three element vector of the lowest source (surface air) temperature, the highest source temperature, and the step, in degrees C: Tsource = [Tcource_lo Tsource_hi dTcource];
% Enter RH at the source as a three element vector of the lowest source RH (from 0 to 1), the highest source RHsource, and the step: RHsource = [RH_lo RH_hi dRH]. 
% OR to use the climatology (reanalysis) RH set RHsource = [];
%
%
% The folowing is an example for a simple run of the model:
% a = 1; b = 0.00525; c = 0; %tuned super saturation parameters
% Tsite = [-40 -20 1];
% Tsource = [10 15 5];
% [T_site, T_source, RH_source, d18O_site, dD_site, d18Oln_site, dDln_site, dxs_site, d17O_xs_site, dlnU_site] = simple_water_istope_model_2018(Tsite, Tsource, [], a, b, c);

%% Housekeeping

if ~exist('closure','var')
    closure = 'local';%default closure assumption is local
end

if ~exist('reanalysis','var')
    reanalysis = 'ncep';%default reanalysis for correlations is NCEP
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
%define pathway. I could make this an input, but I don't know why you'd
%want to cahnge it often
pathway = 'adiabatic';
% pathway = 'isobaric';

%% Set up model parameter gridds

Tsite_lo = Tsite(1);
Tsite_hi = Tsite(2);
dTsite = Tsite(3);

Tsource_lo = Tsource(1);
Tsource_hi = Tsource(2);
dTsource = Tsource(3);

T_site=[Tsite_lo:dTsite:Tsite_hi]; %This creates a vector of temperatures for the end of the temp gradient
T_source=[Tsource_lo:dTsource:Tsource_hi]; %This creates a vector of temperatures for the start of the temp gradient

if isempty(RHsource)
    RH0=[];
    RH_source=[];
else
RHsource_lo = RHsource(1);
RHsource_hi = RHsource(2);
dRHsource = RHsource(3);
RH_source=[RHsource_lo:dRHsource:RHsource_hi];
end

%% Set up misc. other initial parameters

%Uemura Fit
%ln(1 + dD) = -2.85 * 10^-2 * (ln(1 + d18O))^2 +8.47 * ln(1 + d18O) +13.3;%offset isn't needed except for modeling plots
dDln_d18Oln_fit_U= [-2.85*10^-2 8.47 13.3];




dT = 0.1;%This is the temperature step for the temperature gradient
SST0 = [];
%These are a set of variables that we define elsewhere
% RH0 = [];
% closure = 'local';
% reanalysis = 'ncep';
% SH =1;
% a = 1;
% b = 0.00525;
% c = 0;

%% run loops for climatology RH and SST
if isempty(RHsource) && isempty(SST0)
for i=1:length(T_source)
    for j=1:length(T_site)
        T0=T_source(i);
        T_end=T_site(j);
if T0 < T_end
    T0=NaN;
    T_end=NaN;
end

        [dD_v0, d18O_v0, d17Oxs_v0,RHn0,rh0,sst0] = evaporation(T0,SST0,RH0,closure,reanalysis, SH,season);
%         [d18O_p, dD_p, d17O_p, d18O_p_ln, dD_p_ln, d17O_p_ln, dxs, dxs_ln, d17Oxs, T, S, supersat, r_s] = distillation(T0, T_end, dT, dD_v0, d18O_v0, d17Oxs_v0, RHn0, rh0, sst0, a, b, c,pathway);
        [d18O_p, dD_p, d17O_p, d18O_p_ln, dD_p_ln, d17O_p_ln, dxs, dxs_ln, d17Oxs, T, S, supersat, r_s, e_s, f, P] = distillation_2018(T0, T_end, dT, dD_v0, d18O_v0, d17Oxs_v0, RHn0, rh0, sst0, a, b, c,pathway);
%         [d18O_p, dD_p, d17O_p, d18O_p_ln, dD_p_ln, d17O_p_ln, dxs, dxs_ln, d17Oxs, T, S, supersat, r_s, e_s, f, P] = distillation_2019(T0, T_end, dT, dD_v0, d18O_v0, d17Oxs_v0, RHn0, rh0, sst0, P_s, pathway);

model_dln_U= dD_p_ln - polyval(dDln_d18Oln_fit_U,d18O_p_ln) +dDln_d18Oln_fit_U(end);%remove the intercept
%these are various dxs definitions:
% model_dln_2= dD_p_ln - polyval(dDln_d18Oln_fit2,d18O_p_ln) +dDln_d18Oln_fit2(end);
% model_dln_3= dD_p_ln - polyval(dDln_d18Oln_fit3,d18O_p_ln) +dDln_d18Oln_fit3(end);
% model_dln_2plus= dD_p_ln - polyval(dDln_d18Oln_fit_plus2,d18O_p_ln) +dDln_d18Oln_fit_plus2(end);
% model_dln_3plus= dD_p_ln - polyval(dDln_d18Oln_fit_plus3,d18O_p_ln) +dDln_d18Oln_fit_plus3(end);

d17O_xs = ((log(d17O_p./1000+1))-0.528.*((log(d18O_p./1000+1))))*1e6; 

d18O_site(i,j)=d18O_p(end);
dD_site(i,j)=dD_p(end);
d18Oln_site(i,j)=d18O_p_ln(end);
dDln_site(i,j)=dD_p_ln(end);
dxs_site(i,j)=dxs(end);
d17O_xs_site(i,j)=d17O_xs(end);
dlnU_site(i,j)=model_dln_U(end);
% dln2_site(i,j)=model_dln_2(end);
% dln3_site(i,j)=model_dln_3(end);
% dln2plus_site(i,j)=model_dln_2plus(end);
% dln3plus_site(i,j)=model_dln_3plus(end);
r_s_site(i,j)=r_s(end);
P_site(i,j)=P(end);


% e_s_site(i,j)=e_s(end);

% supersat_test(i,j)=supersat(end);
% T_test(i,j)=T_end;

clear d18O_p dD_p d17O_p d18O_p_ln dD_p_ln d17O_p_ln dxs dxs_ln T S ss d17O_xs model_dln_U r_s P e_s%model_dln_2 model_dln_3 model_dln_2plus model_dln_3plus 
    end
end
end

%% run loops for climatology SST and specified RH
if ~isempty(RHsource) && isempty(SST0)
for k=1:length(RH_source)
    for i=1:length(T_source)
        for j=1:length(T_site)
        T0=T_source(i);
        T_end=T_site(j);
        RH0=RH_source(k);
        if T0 < T_end
        T0=NaN;
        T_end=NaN;
        RH0=NaN;
        end
        [dD_v0, d18O_v0, d17Oxs_v0,RHn0,rh0,sst0] = evaporation(T0,SST0,RH0,closure,reanalysis, SH,season);
        [d18O_p, dD_p, d17O_p, d18O_p_ln, dD_p_ln, d17O_p_ln, dxs, dxs_ln, d17Oxs, T, S, supersat, r_s, e_s, f, P] = distillation_2018(T0, T_end, dT, dD_v0, d18O_v0, d17Oxs_v0, RHn0, rh0, sst0, a, b, c,pathway);
%         [d18O_p, dD_p, d17O_p, d18O_p_ln, dD_p_ln, d17O_p_ln, dxs, dxs_ln, d17Oxs, T, S, supersat, r_s, e_s, f, P] = distillation_2019(T0, T_end, dT, dD_v0, d18O_v0, d17Oxs_v0, RHn0, rh0, sst0, P_s, pathway);


model_dln_U= dD_p_ln - polyval(dDln_d18Oln_fit_U,d18O_p_ln) +dDln_d18Oln_fit_U(end);%remove the intercept
% model_dln_2= dD_p_ln - polyval(dDln_d18Oln_fit2,d18O_p_ln) +dDln_d18Oln_fit2(end);
% model_dln_3= dD_p_ln - polyval(dDln_d18Oln_fit3,d18O_p_ln) +dDln_d18Oln_fit3(end);
% model_dln_2plus= dD_p_ln - polyval(dDln_d18Oln_fit_plus2,d18O_p_ln) +dDln_d18Oln_fit_plus2(end);
% model_dln_3plus= dD_p_ln - polyval(dDln_d18Oln_fit_plus3,d18O_p_ln) +dDln_d18Oln_fit_plus3(end);

d17O_xs = ((log(d17O_p./1000+1))-0.528.*((log(d18O_p./1000+1))))*1e6; 

d18O_site(i,j,k)=d18O_p(end);
dD_site(i,j,k)=dD_p(end);
d18Oln_site(i,j,k)=d18O_p_ln(end);
dDln_site(i,j,k)=dD_p_ln(end);
dxs_site(i,j,k)=dxs(end);
d17O_xs_site(i,j,k)=d17O_xs(end);
dlnU_site(i,j,k)=model_dln_U(end);
% dln2_site(i,j)=model_dln_2(end);
% dln3_site(i,j)=model_dln_3(end);
% dln2plus_site(i,j)=model_dln_2plus(end);
% dln3plus_site(i,j)=model_dln_3plus(end);
r_s_site(i,j,k)=r_s(end);
P_site(i,j,k)=P(end);

clear d18O_p dD_p d17O_p d18O_p_ln dD_p_ln d17O_p_ln dxs dxs_ln T S ss d17O_xs model_dln_U r_s P %model_dln_2 model_dln_3 model_dln_2plus model_dln_3plus 
        end
    end
end
end
