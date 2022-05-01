%Range of temperatures for final depositions site
Tsite_hi=0; 
Tsite_lo=-25; 

dTsite=1; 

Tsite=[Tsite_lo Tsite_hi dTsite];

%Range of temperatures for initial evaporation
Tsource_hi=10; 
Tsource_lo=0; 

dTsource=1;

Tsource=[Tsource_lo Tsource_hi dTsource];

%Initial evaporation relative humidity (leave empty to use climatology)

% RHsource=[.8 .8 .1];
RHsource=[];

% SH=0;%Noerthern hemisphere
SH=1;%Southern hemisphere

% closure = 'global';
closure = 'local';

reanalysis = 'ncep';
% reanalysis = 'era';

season = 'annual';
% season = 'DJF';
% season = 'JJA';
% season = 'MAM';
% season = 'SON';

% pathway='isobaric';
pathway = 'adiabatic';

%tune the super saturation parameters
a=1;
b=0.00525;
c=0.00000;

[T_site, T_source, RH_source, d18O_site, dD_site, d18Oln_site, dDln_site, dxs_site, d17O_xs_site, dlnU_site, r_s_site, P_site] = simple_water_isotope_model_2020(Tsite, Tsource, RHsource, a, b, c, closure, reanalysis, SH, season);

figure
hold on
plot(T_site,d18O_site,'k')

figure
hold on
plot(T_site,d17O_xs_site,'k')