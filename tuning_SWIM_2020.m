% cd('/Volumes/Macintosh HD')%edit 09 2018
% cd('/Users/Bradley/Documents/work/wais_dxs_paper/WDC_dxs_2016/redefining_dxs_in_antarctica/')
% cd('Users/Bradley/Documents/work/wais_dxs_paper/WDC_dxs_2016/redefining_dxs_in_antarctica/')
% cd('/Volumes/Macintosh HD/Users/Bradley/Documents/work/wais_dxs_paper/WDC_dxs_2016/redefining_dxs_in_antarctica')
% cd('/Volumes/Macintosh HD/Users/Bradley/Documents/work/wais_dxs_paper/WDC_dxs_2016/redefining_dxs_in_antarctica')
% compile_surface_data_4
% compile_surface_data_5
 cd '/Users/bradley/work/Data/iso_data'
compile_surface_data_2020

compile_deep_d17O_data_2020
% close all

Colors=[0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

% cd /Users/Bradley/Documents/MATLAB/simple_water_isotope_model
% cd Users/Bradley/Documents/MATLAB/simple_water_isotope_model
% cd '/Users/bradleymarkle/Documents/MATLAB/simple_water_isotope_model'
cd '/Users/bradley/Documents/MATLAB/simple_water_isotope_model'

%%

fig(1,'units','inches','width',10,'height',8,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(landais_d18O_cor,landais_d17Oxs,'ko')
plot(wais_snowpit_d18O,wais_snowpit_d17Oxs,'ko')
% plot(modern_comp_d18O,modern_comp_O17xs_PD,'k.')
% plot(modern_comp_d18O,modern_comp_O17xs_EH,'kx')
plot(DomeA_d18O,DomeA_17O,'k.')
xlabel(['\delta^{18}O (',char(8240),')'])
ylabel(['\delta^{17}O_{excess} (per meg)'])
plot(wais_d18O,wais_d17Oxs,'k.')
plot(siple_d18O,siple_d17Oxs,'k.')
plot(taylor_d18O,taylor_d17Oxs,'k.')
plot(vostok_d18O,vostok_d17Oxs,'k.')
plot(vostok_d18O,vostok_d17Oxs_org,'ko')
% plot(vostok_d18O,vostok_d17Oxs_offset,'kx')
plot(talos_d18O,talos_d17Oxs,'k.')
plot(edc_d18O,edc_d17Oxs,'k.')


fig(2,'units','inches','width',10,'height',8,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(GNIP_d18O_ln,GNIP_dlnU,'k.')
plot(MD08_d18O_ln,MD08_dlnU,'k.')
plot(Dome_A_d18O_ln,Dome_A_dlnU,'ko')
plot(DomeA_d18O_ln,DomeA_dlnU,'k.')
plot(Dahe_d18O_ln,Dahe_dlnU,'k.')
% plot(GMWL_d18O_ln,GMWL_dlnU,'k')
plot(wais_snowpit_d18O_ln,wais_snowpit_dlnU,'k.')
plot(modern_comp_d18O_ln,modern_comp_dlnU,'k.')
plot(landais_d18O_ln,landais_dlnU,'ko')
plot(JIRP_d18O_ln,JIRP_dlnU,'k.')
plot(wais_shallow_cores_d18O_ln,wais_shallow_cores_dlnU,'k.')
xlabel(['\delta^{18}O_{ln} (',char(8240),')'])
ylabel(['d_{ln}^{U} (',char(8240),')'])

fig(3,'units','inches','width',10,'height',8,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(GNIP_d18O_ln,GNIP_dD_ln,'k.')
plot(MD08_d18O_ln,MD08_dD_ln,'k.')
plot(Dome_A_d18O_ln,Dome_A_dD_ln,'ko')
plot(DomeA_d18O_ln,DomeA_dD_ln,'k.')
plot(Dahe_d18O_ln,Dahe_dD_ln,'k.')
% plot(GMWL_d18O_ln,GMWL_dlnU,'k')
plot(wais_snowpit_d18O_ln,wais_snowpit_dD_ln,'k.')
plot(modern_comp_d18O_ln,modern_comp_dD_ln,'k.')
plot(landais_d18O_ln,landais_dD_ln,'ko')
plot(JIRP_d18O_ln,JIRP_dD_ln,'c.')
plot(wais_shallow_cores_d18O_ln,wais_shallow_cores_dD_ln,'k.')
xlabel(['\delta^{18}O_{ln} (',char(8240),')'])
ylabel(['\deltaD_{ln} (',char(8240),')'])

fig(4,'units','inches','width',10,'height',8,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(GNIP_d18O,GNIP_dD,'.')
plot(MD08_d18O,MD08_dD,'.')
plot(Dome_A_d18O,Dome_A_dD,'o')
plot(DomeA_d18O,DomeA_dD,'.')
plot(Dahe_d18O,Dahe_dD,'.')
plot(wais_snowpit_d18O,wais_snowpit_dD,'.')
plot(modern_comp_d18O,modern_comp_dD,'.')
% plot(landais_d18O,landais_dD,'o')
plot(JIRP_d18O,JIRP_dD,'.')
xlabel(['\delta^{18}O (',char(8240),')'])
ylabel(['\deltaD (',char(8240),')'])
plot([-60 0],8.*[-60 0]+10,'k')


fig(400,'units','inches','width',10,'height',8,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(GNIP_dD_ln,GNIP_dlnU,'k.')
plot(MD08_dD_ln,MD08_dlnU,'k.')
plot(Dome_A_dD_ln,Dome_A_dlnU,'ko')
plot(DomeA_dD_ln,DomeA_dlnU,'k.')
plot(Dahe_dD_ln,Dahe_dlnU,'k.')
% plot(GMWL_dD_ln,GMWL_dlnU,'k')
plot(wais_snowpit_dD_ln,wais_snowpit_dlnU,'k.')
plot(modern_comp_dD_ln,modern_comp_dlnU,'k.')
plot(landais_dD_ln,landais_dlnU,'ko')
plot(JIRP_dD_ln,JIRP_dlnU,'c.')
plot(wais_shallow_cores_dD_ln,wais_shallow_cores_dlnU,'k.')
xlabel(['\delta D_{ln} (',char(8240),')'])
ylabel(['d_{ln}^{U} (',char(8240),')'])


%%


load('/Users/bradley/Documents/MATLAB/simple_water_isotope_model/data/GNIP/GNIP_Cryo_short.txt');
GNIP_Cryo_d18O=GNIP_Cryo_short(:,5);
GNIP_Cryo_dD=GNIP_Cryo_short(:,6);
GNIP_Cryo_T=GNIP_Cryo_short(:,7);

load('/Users/bradley/Documents/MATLAB/simple_water_isotope_model/data/GNIP/GNIP_monthly_SH_short.txt');
GNIP_SH_d18O=GNIP_monthly_SH_short(:,8);
GNIP_SH_dD=GNIP_monthly_SH_short(:,9);
GNIP_SH_T=GNIP_monthly_SH_short(:,11);

load('/Users/bradley/Documents/MATLAB/simple_water_isotope_model/data/GNIP/GNIP_monthly_NH_combined_short.txt');
GNIP_NH_d18O=GNIP_monthly_NH_combined_short(:,8);
GNIP_NH_dD=GNIP_monthly_NH_combined_short(:,9);
GNIP_NH_T=GNIP_monthly_NH_combined_short(:,11);

% figure;
% hold on
% plot(GNIP_Cryo_T,GNIP_Cryo_d18O,'.')
% plot(GNIP_NH_T,GNIP_NH_d18O,'.')
% plot(GNIP_SH_T,GNIP_SH_d18O,'.')

% % % combine all into one
GNIP_d18O_monthly=[GNIP_Cryo_d18O; GNIP_SH_d18O; GNIP_NH_d18O];
GNIP_dD_monthly=[GNIP_Cryo_dD; GNIP_SH_dD; GNIP_NH_dD];
GNIP_T_monthly=[GNIP_Cryo_T; GNIP_SH_T; GNIP_NH_T];

% combine SH and cryo into one
% GNIP_d18O=[GNIP_Cryo_d18O; GNIP_SH_d18O];
% GNIP_dD=[GNIP_Cryo_dD; GNIP_SH_dD];
% GNIP_T=[GNIP_Cryo_T; GNIP_SH_T];

GNIP_d18O_ln_monthly=log(GNIP_d18O_monthly./1000+1)*1000;
GNIP_dD_ln_monthly=log(GNIP_dD_monthly./1000+1)*1000;

% [PU,So] = polyfit([GNIP_d18O_ln; MD08_d18O_ln],[GNIP_dD_ln; MD08_dD_ln],2);%this give Uemuras fit. !!NOT EXACTly
PU = [-2.85*10^-2 8.47 13.3];
GNIP_dlnU_monthly = GNIP_dD_ln_monthly-(polyval(PU,GNIP_d18O_ln_monthly)-PU(end));

GNIP_dln_monthly=GNIP_dD_ln_monthly-(-2.85*10^(-2) .* GNIP_d18O_ln_monthly.^2 + 8.47 .* GNIP_d18O_ln_monthly);
GNIP_dxs_monthly=GNIP_dD_monthly-8*GNIP_d18O_monthly;


figure(2)
hold on
% plot(GNIP_d18O_ln_monthly,GNIP_dlnU_monthly,'b.')
plot(GNIP_d18O_ln_monthly,GNIP_dln_monthly,'k.')


fig(100,'units','inches','width',10,'height',8,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(GNIP_d18O_monthly,GNIP_dln_monthly,'k.')
% plot(GNIP_d18O,GNIP_dlnU,'.','Color',[0.8 0.8 0.8])
plot(GNIP_d18O,GNIP_dlnU,'.','Color','r')
% plot(MD08_d18O,MD08_dlnU,'.','Color',[0.6 0.6 0.6])%repeat
% plot(Dome_A_dD_ln,Dome_A_dlnU,'ko')
plot(DomeA_d18O,DomeA_dlnU,'.','Color',[0.0 0.0 0.1])
% plot(Dahe_d18O,Dahe_dlnU,'.','Color',[1 0.8 0.8])%repeat
% plot(GMWL_dD_ln,GMWL_dlnU,'k')
plot(wais_snowpit_d18O,wais_snowpit_dlnU,'k.')
plot(modern_comp_d18O,modern_comp_dlnU,'k.')
plot(landais_d18O_cor,landais_dlnU,'k.')
plot(JIRP_d18O,JIRP_dlnU,'k.')
plot(wais_shallow_cores_d18O,wais_shallow_cores_dlnU,'k.')
xlabel(['\delta^{18}O (',char(8240),')'])
ylabel(['d_{ln}^{U} (',char(8240),')'])


fig(1000,'units','inches','width',13,'height',11,'font','Helvetica','fontsize',16,'border','on');
subplot(211)
hold on
plot(GNIP_d18O_monthly,GNIP_dln_monthly,'k.')
% plot(GNIP_d18O,GNIP_dlnU,'.','Color',[0.8 0.8 0.8])
plot(GNIP_d18O,GNIP_dlnU,'.','Color','r')
% plot(MD08_d18O,MD08_dlnU,'.','Color',[0.6 0.6 0.6])%repeat
% plot(Dome_A_dD_ln,Dome_A_dlnU,'ko')
plot(DomeA_d18O,DomeA_dlnU,'.','Color',[0.0 0.0 0.1])
% plot(Dahe_d18O,Dahe_dlnU,'.','Color',[1 0.8 0.8])%repeat
% plot(GMWL_dD_ln,GMWL_dlnU,'k')
plot(wais_snowpit_d18O,wais_snowpit_dlnU,'k.')
plot(modern_comp_d18O,modern_comp_dlnU,'k.')
plot(landais_d18O_cor,landais_dlnU,'k.')
plot(JIRP_d18O,JIRP_dlnU,'c.')
plot(wais_shallow_cores_d18O,wais_shallow_cores_dlnU,'k.')
xlabel(['\delta^{18}O (',char(8240),')'])
ylabel(['d_{ln}^{U} (',char(8240),')'])

subplot(212)
hold on
plot(GNIP_d18O_monthly,GNIP_dxs_monthly,'k.')
% plot(GNIP_d18O,GNIP_dxs,'.','Color','r')
plot(DomeA_d18O,DomeA_dxs,'.','Color',[0.0 0.0 0.1])
plot(wais_snowpit_d18O,wais_snowpit_dD-8.*wais_snowpit_d18O,'k.')
% plot(modern_comp_d18O,modern_comp_dxs,'k.')
plot(JIRP_d18O,JIRP_dD-8.*JIRP_d18O,'c.')
plot(wais_shallow_cores_d18O,wais_shallow_cores_dD-8.*wais_shallow_cores_d18O,'k.')
xlabel(['\delta^{18}O (',char(8240),')'])
ylabel(['d_{xs} (',char(8240),')'])

%%
Tsite_hi=24; 
Tsite_lo=-57; 
% dTsite=1; 
dTsite=1; 

Tsite=[Tsite_lo Tsite_hi dTsite];

% Tsource_hi=28; 
% Tsource_lo=0; 
Tsource_hi=28; 
Tsource_lo=2; 
dTsource=4;
% dTsource=2;

Tsource=[Tsource_lo Tsource_hi dTsource];

RHsource=[];
% RHsource=[.9 .9 .1];
% RHsource=[.8 .8 .1];
% RHsource=[.7 .7 .1];
% RHsource=[0.7 0.9 0.1];

% SH=0;

SH=1;
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

%tun the super saturation parameters
%for adiabatic
% a=1;
% % b=0.00525;%works well for ncep,local, BEST
% % b=0.005;%works well for ncep,local
% % b=0.021;%for very low diffusivity
% % b=0.025;%
% b=0.00525;
% % b= 0.0055;
% c=0.00001;
% % c=-0.0000;
% % c=0.000005;
% % c=0;


%for isobaric
a=1;
% a=1.001;
% b=0.00525;
% b=0.003;
% b=0.0030;
b=0.00507;

% c=0.00001;
c=0.00000;

% [T_site, T_source, RH_source, d18O_site, dD_site, d18Oln_site, dDln_site, dxs_site, d17O_xs_site, dlnU_site, r_s_site] = simple_water_istope_model(Tsite, Tsource, RHsource, a, b, c, closure, reanalysis, SH, season);
% [T_site, T_source, RH_source, d18O_site, dD_site, d18Oln_site, dDln_site, dxs_site, d17O_xs_site, dlnU_site, r_s_site] = simple_water_istope_model_testing(Tsite, Tsource, RHsource, a, b, c, closure, reanalysis, SH, season);
% [T_site, T_source, RH_source, d18O_site, dD_site, d18Oln_site, dDln_site, dxs_site, d17O_xs_site, dlnU_site, r_s_site] = simple_water_istope_model_2018(Tsite, Tsource, RHsource, a, b, c, closure, reanalysis, SH, season);
% [T_site, T_source, RH_source, d18O_site, dD_site, d18Oln_site, dDln_site, dxs_site, d17O_xs_site, dlnU_site, r_s_site] = simple_water_istope_model_2019(Tsite, Tsource, RHsource, a, b, c, closure, reanalysis, SH, season);
% [T_site, T_source, RH_source, d18O_site, dD_site, d18Oln_site, dDln_site, dxs_site, d17O_xs_site, dlnU_site, r_s_site] = simple_water_isotope_model_2019(Tsite, Tsource, RHsource, a, b, c, closure, reanalysis, SH, season);
[T_site, T_source, RH_source, d18O_site, dD_site, d18Oln_site, dDln_site, dxs_site, d17O_xs_site, dlnU_site, r_s_site] = simple_water_isotope_model_2020(Tsite, Tsource, RHsource, a, b, c, closure, reanalysis, SH, season);

[m, n, q] =size(d18O_site);

Color=[1 .5 .5];

if q > 1
    for i = 1:q
    figure(1)
    plot(squeeze(d18O_site(:,:,i)),squeeze(d17O_xs_site(:,:,i)),'o','Color',Color)
    figure(2)
    plot(squeeze(d18Oln_site(:,:,i)),squeeze(dlnU_site(:,:,i)),'o','Color',Color)
    figure(3)
    plot(squeeze(d18Oln_site(:,:,i)),squeeze(dDln_site(:,:,i)),'o','Color',Color)
    figure(4)
    plot(squeeze(dDln_site(:,:,i)),squeeze(dlnU_site(:,:,i)),'o','Color',Color)
    end
else
    figure(1)
    plot(d18O_site,d17O_xs_site,'o','Color',Color)
    figure(2)
    plot(d18Oln_site,dlnU_site,'o','Color',Color)
    figure(3)
    plot(d18Oln_site,dDln_site,'o','Color',Color)
     figure(4)
    plot(dDln_site,dlnU_site,'o','Color',Color)
    
    
    figure(101)
    hold on
   plot(d18Oln_site,dlnU_site,'o','Color',Color)

        figure(201)
    hold on
   plot(d18O_site,dxs_site,'o','Color',Color)


   figure(1000)
   subplot(211)
       scatter(d18O_site,dlnU_site,10,'filled','MarkerFaceColor',Color)
       plot(d18O_site',dlnU_site','Color',Color)

      subplot(212)
    scatter(d18O_site,dxs_site,10,'filled','MarkerFaceColor',Color)
       plot(d18O_site',dxs_site','Color',Color)

    
 figure(5)
 hold on
 contour(d18O_site,dlnU_site,repmat(T_source,n,1)',T_source,'ShowText','on','Linewidth',2)
% colormap(flipud(newcmap_RdBu))
    
% T_source_grid=repmat(T_source,n,1);T_source_grid=T_source_grid';
% T_site_grid=repmat(T_site,m,1);
%     figure(1)
%     scatter(d18O_site(:),d17O_xs_site(:),1e5.*r_s_site(:),'ko')
%     figure(2)
%     scatter(d18Oln_site(:),dlnU_site(:),1e5.*r_s_site(:),'ko')
%     figure(3)
%     scatter(d18Oln_site(:),dDln_site(:),1e5.*r_s_site(:),'ko')
end


[newcmap_RdBu]=(cbrewer_BRM('div', 'RdBu',50));



fig(100,'units','inches','width',10,'height',8,'font','Helvetica','fontsize',16,'border','on');
hold on
 contour(d18O_site,dlnU_site,repmat(T_source,n,1)',T_source,'ShowText','on','Linewidth',2,'Linewidth',3)
%  contour(d18O_site,dlnU_site,repmat(T_source,n,1)'-mean(T_source),T_source-mean(T_source),'ShowText','on','Linewidth',3)
caxis([-28 28])
colormap(flipud(newcmap_RdBu))
xlabel(['\delta^{18}O (',char(8240),')'])
ylabel(['d_{ln} (',char(8240),')'])
colorbar

fig(44,'units','inches','width',10,'height',8,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(T_site,d18O_site,'.','Color',Color)
% plot(T_site./.69,d18O_site,'k.')
% plot(T_site,d18O_site(T_source==12,:),'Color',Color,'linewidth',2)
ylabel(['\delta^{18}O_{ln} (',char(8240),')'])
xlabel('Condensation Temperature (\circC)')

 figure(8)
 hold on
 contour(d18O_site,dlnU_site,repmat(T_site,m,1),T_site,'ShowText','on','Linewidth',2)
 
%  
%  figure;
%  hold on
% xxx= repmat(T_site,m,1);
% %  plot(xxx,r_s_site,'.')
%   plot(xxx(:,1:end-1),diff(r_s_site,[],2),'.')
 %%
 
 
yyy=repmat(T_source,n,1)';  
xxx= repmat(T_site,m,1);


fig(200,'units','inches','width',10,'height',8,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(GNIP_d18O_monthly,GNIP_dln_monthly,'k.')
plot(GNIP_d18O,GNIP_dlnU,'.','Color','k')
plot(DomeA_d18O,DomeA_dlnU,'.','Color','k')
plot(wais_snowpit_d18O,wais_snowpit_dlnU,'k.')
plot(modern_comp_d18O,modern_comp_dlnU,'k.')
plot(landais_d18O_cor,landais_dlnU,'k.')
plot(JIRP_d18O,JIRP_dlnU,'k.')
plot(wais_shallow_cores_d18O,wais_shallow_cores_dlnU,'k.')
xlabel(['\delta^{18}O (',char(8240),')'])
ylabel(['d_{ln}^{U} (',char(8240),')'])

% scatter(d18O_site(:),dlnU_site(:),r_s_site(:).*10e3,yyy(:))
 scatter(d18O_site(:),dlnU_site(:),20,yyy(:))

figure(300)
subplot(211)
hold on
plot(GNIP_d18O_monthly,GNIP_dln_monthly,'k.')
plot(GNIP_d18O,GNIP_dlnU,'.','Color','k')
plot(DomeA_d18O,DomeA_dlnU,'.','Color','k')
plot(wais_snowpit_d18O,wais_snowpit_dlnU,'k.')
plot(modern_comp_d18O,modern_comp_dlnU,'k.')
plot(landais_d18O_cor,landais_dlnU,'k.')
plot(JIRP_d18O,JIRP_dlnU,'k.')
plot(wais_shallow_cores_d18O,wais_shallow_cores_dlnU,'k.')
xlabel(['\delta^{18}O (',char(8240),')'])
ylabel(['d_{ln}^{U} (',char(8240),')'])
scatter(d18O_site(:),dlnU_site(:),20,yyy(:))
subplot(212)
hold on
plot(GNIP_d18O_monthly,GNIP_dln_monthly,'k.')
plot(GNIP_d18O,GNIP_dlnU,'.','Color','k')
plot(DomeA_d18O,DomeA_dlnU,'.','Color','k')
plot(wais_snowpit_d18O,wais_snowpit_dlnU,'k.')
plot(modern_comp_d18O,modern_comp_dlnU,'k.')
plot(landais_d18O_cor,landais_dlnU,'k.')
plot(JIRP_d18O,JIRP_dlnU,'k.')
plot(wais_shallow_cores_d18O,wais_shallow_cores_dlnU,'k.')
xlabel(['\delta^{18}O (',char(8240),')'])
ylabel(['d_{ln}^{U} (',char(8240),')'])
scatter(d18O_site(:),dlnU_site(:),20,xxx(:))

%%

[cmap_Tsource]=(cbrewer_BRM('seq', 'YlOrRd',100));
cmap_Tsource_cliped=cmap_Tsource(20:100,:);
[cmap_Tsite]=flipud(cbrewer_BRM('seq', 'Blues',100));



SS=20;
fig(4,'units','inches','width',10,'height',8,'font','Helvetica','fontsize',16,'border','on');
hold on
scatter(landais_d18O_cor,landais_d17Oxs,SS,'filled','MarkerFaceColor','k')
plot(wais_snowpit_d18O,wais_snowpit_d17Oxs,'*b')
scatter(modern_comp_d18O,modern_comp_O17xs_PD,SS,'filled','MarkerFaceColor','k')
scatter(modern_comp_d18O,modern_comp_O17xs_EH,SS,'filled','MarkerFaceColor','k')
scatter(DomeA_d18O,DomeA_17O,SS,'filled','MarkerFaceColor','k')
scatter(vostok_d18O,vostok_d17Oxs_offset,SS,'filled','MarkerFaceColor','c')
scatter(edc_d18O,edc_d17Oxs,SS,'filled','MarkerFaceColor','c')


% plot(landais_d18O_cor,landais_d17Oxs,'ko')
% plot(wais_snowpit_d18O,wais_snowpit_d17Oxs,'ko')
% plot(modern_comp_d18O,modern_comp_O17xs_PD,'k.')
% plot(modern_comp_d18O,modern_comp_O17xs_EH,'ko')
% plot(DomeA_d18O,DomeA_17O,'k.')
xlabel(['\delta^{18}O (',char(8240),')'])
ylabel(['\delta^{17}O_{excess} (per meg)'])
% plot(wais_d18O,wais_d17Oxs,'.')
% plot(siple_d18O,siple_d17Oxs,'.')
% plot(taylor_d18O,taylor_d17Oxs,'.')
% plot(vostok_d18O,vostok_d17Oxs,'k.')
% plot(vostok_d18O,vostok_d17Oxs_org,'ko')
% plot(vostok_d18O,vostok_d17Oxs_offset,'k.')
% plot(talos_d18O,talos_d17Oxs,'.')
% plot(edc_d18O,edc_d17Oxs,'k.')


fig(5,'units','inches','width',10,'height',8,'font','Helvetica','fontsize',16,'border','on');
hold on
scatter(GNIP_d18O_ln,GNIP_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(MD08_d18O_ln,MD08_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(DomeA_d18O_ln,DomeA_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(Dahe_d18O_ln,Dahe_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(wais_snowpit_d18O_ln,wais_snowpit_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(modern_comp_d18O_ln,modern_comp_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(JIRP_d18O_ln,JIRP_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(wais_shallow_cores_d18O_ln,wais_shallow_cores_dlnU,SS,'filled','MarkerFaceColor','k')
% plot(GNIP_d18O_ln,GNIP_dlnU,'k.')
% plot(MD08_d18O_ln,MD08_dlnU,'k.')
% plot(Dome_A_d18O_ln,Dome_A_dlnU,'o')
% plot(DomeA_d18O_ln,DomeA_dlnU,'k.')
% plot(Dahe_d18O_ln,Dahe_dlnU,'k.')
% % plot(GMWL_d18O_ln,GMWL_dlnU,'k')
% plot(wais_snowpit_d18O_ln,wais_snowpit_dlnU,'ko')
% plot(modern_comp_d18O_ln,modern_comp_dlnU,'k.')
% plot(landais_d18O_ln,landais_dlnU,'o')
% plot(JIRP_d18O_ln,JIRP_dlnU,'k.')
% plot(wais_shallow_cores_d18O_ln,wais_shallow_cores_dlnU,'k.')
xlabel(['\delta''','^{18}O (',char(8240),')'])
ylabel(['d_{ln} (',char(8240),')'])

fig(6,'units','inches','width',10,'height',8,'font','Helvetica','fontsize',16,'border','on');
hold on
scatter(GNIP_d18O_ln,GNIP_dD_ln,SS,'filled','MarkerFaceColor','k')
scatter(MD08_d18O_ln,MD08_dD_ln,SS,'filled','MarkerFaceColor','k')
scatter(DomeA_d18O_ln,DomeA_dD_ln,SS,'filled','MarkerFaceColor','k')
scatter(Dahe_d18O_ln,Dahe_dD_ln,SS,'filled','MarkerFaceColor','k')
scatter(wais_snowpit_d18O_ln,wais_snowpit_dD_ln,SS,'filled','MarkerFaceColor','k')
scatter(modern_comp_d18O_ln,modern_comp_dD_ln,SS,'filled','MarkerFaceColor','k')
scatter(JIRP_d18O_ln,JIRP_dD_ln,SS,'filled','MarkerFaceColor','k')
scatter(wais_shallow_cores_d18O_ln,wais_shallow_cores_dD_ln,SS,'filled','MarkerFaceColor','k')
% plot(GNIP_d18O_ln,GNIP_dD_ln,'k.')
% plot(MD08_d18O_ln,MD08_dD_ln,'k.')
% % plot(Dome_A_d18O_ln,Dome_A_dD_ln,'o')
% plot(DomeA_d18O_ln,DomeA_dD_ln,'k.')
% plot(Dahe_d18O_ln,Dahe_dD_ln,'k.')
% % plot(GMWL_d18O_ln,GMWL_dlnU,'k')
% % plot(wais_snowpit_d18O_ln,wais_snowpit_dD_ln,'k.')
% plot(modern_comp_d18O_ln,modern_comp_dD_ln,'k.')
% % plot(landais_d18O_ln,landais_dD_ln,'o')
% plot(JIRP_d18O_ln,JIRP_dD_ln,'k.')
% plot(wais_shallow_cores_d18O_ln,wais_shallow_cores_dD_ln,'k.')
xlabel(['\delta''','^{18}O (',char(8240),')'])
ylabel(['\delta''','D(',char(8240),')'])

fig(7,'units','inches','width',10,'height',8,'font','Helvetica','fontsize',16,'border','on');
hold on
scatter(GNIP_d18O,GNIP_dD-8.*GNIP_d18O,SS,'filled','MarkerFaceColor','k')
scatter(MD08_d18O,MD08_dD-8.*MD08_d18O,SS,'filled','MarkerFaceColor','k')
scatter(DomeA_d18O,DomeA_dD-8.*DomeA_d18O,SS,'filled','MarkerFaceColor','k')
scatter(Dahe_d18O,Dahe_dD-8.*Dahe_d18O,SS,'filled','MarkerFaceColor','k')
scatter(wais_snowpit_d18O,wais_snowpit_dD-8.*wais_snowpit_d18O,SS,'filled','MarkerFaceColor','k')
scatter(modern_comp_d18O,modern_comp_dD-8.*modern_comp_d18O,SS,'filled','MarkerFaceColor','k')
scatter(JIRP_d18O,JIRP_dD-8.*JIRP_d18O,SS,'filled','MarkerFaceColor','k')
scatter(wais_shallow_cores_d18O,wais_shallow_cores_dD-8.*wais_shallow_cores_d18O,SS,'filled','MarkerFaceColor','k')
% plot(GNIP_d18O,GNIP_dD-8.*GNIP_d18O,'k.')
% plot(MD08_d18O,MD08_dD-8.*MD08_d18O,'k.')
% plot(DomeA_d18O,DomeA_dD-8.*DomeA_d18O,'k.')
% plot(Dahe_d18O,Dahe_dD-8.*Dahe_d18O,'k.')
% plot(wais_snowpit_d18O,wais_snowpit_dD-8.*wais_snowpit_d18O,'k.')
% plot(modern_comp_d18O,modern_comp_dD-8.*modern_comp_d18O,'k.')
% plot(JIRP_d18O,JIRP_dD-8.*JIRP_d18O,'k.')
% plot(wais_shallow_cores_d18O,wais_shallow_cores_dD-8.*wais_shallow_cores_d18O,'k.')
xlabel(['\delta''','^{18}O (',char(8240),')'])
ylabel(['d_{xs}(',char(8240),')'])

% fig(4,'units','inches','width',10,'height',8,'font','Helvetica','fontsize',16,'border','on');
% hold on
% plot(GNIP_d18O,GNIP_dD,'.')
% plot(MD08_d18O,MD08_dD,'.')
% plot(Dome_A_d18O,Dome_A_dD,'o')
% plot(DomeA_d18O,DomeA_dD,'.')
% plot(Dahe_d18O,Dahe_dD,'.')
% plot(wais_snowpit_d18O,wais_snowpit_dD,'.')
% plot(modern_comp_d18O,modern_comp_dD,'.')
% % plot(landais_d18O,landais_dD,'o')
% plot(JIRP_d18O,JIRP_dD,'.')
% xlabel(['\delta^{18}O_{ln} (',char(8240),')'])
% ylabel(['\deltaD_{ln} (',char(8240),')'])
% plot([-60 0],8.*[-60 0],'k')


fig(8,'units','inches','width',10,'height',8,'font','Helvetica','fontsize',16,'border','on');
hold on
scatter(GNIP_d18O_ln,GNIP_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(MD08_d18O_ln,MD08_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(DomeA_d18O_ln,DomeA_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(Dahe_d18O_ln,Dahe_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(wais_snowpit_d18O_ln,wais_snowpit_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(modern_comp_d18O_ln,modern_comp_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(JIRP_d18O_ln,JIRP_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(wais_shallow_cores_d18O_ln,wais_shallow_cores_dlnU,SS,'filled','MarkerFaceColor','k')
%%

[m, n, q] =size(d18O_site);

T_source_grid=repmat(T_source,n,1);T_source_grid=T_source_grid';
T_site_grid=repmat(T_site,m,1);
S=30;
LW=1.5;
    figure(4)
    scatter(d18O_site(:),d17O_xs_site(:),S,T_source_grid(:),'LineWidth',LW)
    colormap(cmap_Tsource_cliped)
    colorbar
    figure(5)
    scatter(d18Oln_site(:),dlnU_site(:),S,T_source_grid(:),'LineWidth',LW)
    colormap(cmap_Tsource_cliped)
    colorbar
    figure(6)
    scatter(d18Oln_site(:),dDln_site(:),S,T_source_grid(:),'LineWidth',LW)
    colormap(cmap_Tsource_cliped)
    colorbar

    figure(7)
    scatter(d18O_site(:),dxs_site(:),S,T_source_grid(:),'LineWidth',LW)
    colormap(cmap_Tsource_cliped)
    colorbar
    
    
    
    figure
        scatter(d18O_site(:),dlnU_site(:),S,T_source_grid(:)-nanmean(T_source>5),'LineWidth',LW)
%     colormap(cmap_Tsource_cliped)
    colorbar
    
    %%
    
    
%     load('/Volumes/Macintosh HD/Users/Bradley/Documents/work/Data/Masson_Delmotte_surface_short_2.txt')
    load('/Users/bradley/work/Data/iso_data/Masson_Delmotte_surface_short_2.txt')
MD08_d18O=Masson_Delmotte_surface_short_2(:,9);
MD08_dD=Masson_Delmotte_surface_short_2(:,8);


MD08_lat=Masson_Delmotte_surface_short_2(:,1);
MD08_elev=Masson_Delmotte_surface_short_2(:,3);
MD08_temp=Masson_Delmotte_surface_short_2(:,6);
MD08_acc=Masson_Delmotte_surface_short_2(:,7);

% figure;plot(MD08_lat,MD08_d18O,'.')
% 
% figure;plot(MD08_elev,MD08_d18O,'.')
% figure;plot(MD08_temp,MD08_d18O,'.')

% MD_T_dD_d18O=load('/Users/Bradley/Documents/work/wais_dxs_paper/WDC_dxs_2016/redefining_dxs_in_antarctica/papers_and_data/MD08_T_dD_d18O_data.txt');
% T_MD=MD_T_dD_d18O(:,1);
% dD_MD=MD_T_dD_d18O(:,2);
% d18O_MD=MD_T_dD_d18O(:,3);

T_MD=MD08_temp;
dD_MD=MD08_dD;
d18O_MD=MD08_d18O;
[T_c] = Ts_to_Tc(T_MD);


[newcmap]=flipud(cbrewer_BRM('seq', 'YlOrRd',length(T_source)+1));
newcmap=flipud(newcmap);
fig('units','inches','width',13,'height',6,'font','Helvetica','fontsize',16,'border','on');
subplot(121)
hold on
h5=plot(T_c,d18O_MD,'o');
% h1=scatter(T_c,d18O_MD,'filled');

ylabel(['\delta^{18}O (',char(8240),')'])
xlabel('Condensation Temperature (\circC)')

subplot(122)
hold on
plot(T_c,dD_MD,'o')
% scatter(T_c,dD_MD,'filled')

ylabel(['\deltaD (',char(8240),')'])
xlabel('Condensation Temperature (\circC)')
for i = 1:length(T_source)
subplot(121)
hold on
plot(T_site,d18O_site(i,:),'Color',newcmap(i+1,:),'linewidth',2.5)
subplot(122)
hold on
plot(T_site,dD_site(i,:),'Color',newcmap(i+1,:),'linewidth',2.5)
end
subplot(121)
% legend([h5],{'Observations (MD08)'},'Location','SouthEast')
text(-58,0,'\bfa)')
% xlim([-75 10])
% ylim([-65 0])

subplot(122)
h1=plot(0,0,'Color',newcmap(2,:),'linewidth',2.5);
h2=plot(0,0,'Color',newcmap(3,:),'linewidth',2.5);
h3=plot(0,0,'Color',newcmap(4,:),'linewidth',2.5);
h4=plot(0,0,'Color',newcmap(5,:),'linewidth',2.5);
% legend([h1 h2 h3 h4],...
%     {['T_{0}=',num2str(T_source(1)),'\circ'];['T_{0}=',num2str(T_source(2)),'\circ'];['T_{0}=',num2str(T_source(3)),'\circ'];['T_{0}=',num2str(T_source(4)),'\circ']}...
%     ,'Location','SouthEast');
legend([h5 h1 h2 h3 h4],...
    {['Observations (MD08)']; ['T_{0}=',num2str(T_source(1)),'\circ'];['T_{0}=',num2str(T_source(2)),'\circ'];['T_{0}=',num2str(T_source(3)),'\circ'];['T_{0}=',num2str(T_source(4)),'\circ']}...
    ,'Location','SouthEast');

text(-58,0,'\bfb)')
% xlim([-75 10])
% ylim([-500 -50])




[newcmap]=flipud(cbrewer_BRM('seq', 'YlOrRd',length(T_source)+1));
newcmap=flipud(newcmap);
fig('units','inches','width',13,'height',6,'font','Helvetica','fontsize',16,'border','on');
subplot(121)
hold on
% h5=plot(T_c,d18O_MD,'o');
% h1=scatter(T_c,d18O_MD,'filled');

ylabel(['\delta^{18}O (',char(8240),')'])
xlabel('Condensation Temperature (\circC)')

subplot(122)
hold on
% plot(T_c,dD_MD,'o')
% scatter(T_c,dD_MD,'filled')

ylabel(['\deltaD (',char(8240),')'])
xlabel('Condensation Temperature (\circC)')
for i = 1:length(T_source)
subplot(121)
hold on
plot(T_site,d18O_site(i,:),'Color',newcmap(i+1,:),'linewidth',2.5)
subplot(122)
hold on
plot(T_site,dD_site(i,:),'Color',newcmap(i+1,:),'linewidth',2.5)
end
subplot(121)
% legend([h5],{'Observations (MD08)'},'Location','SouthEast')
text(-58,0,'\bfa)')
% xlim([-75 10])
% ylim([-65 0])

subplot(122)
h1=plot(0,0,'Color',newcmap(2,:),'linewidth',2.5);
h2=plot(0,0,'Color',newcmap(3,:),'linewidth',2.5);
h3=plot(0,0,'Color',newcmap(4,:),'linewidth',2.5);
h4=plot(0,0,'Color',newcmap(5,:),'linewidth',2.5);
legend([h1 h2 h3 h4],...
    {['T_{0}=',num2str(T_source(1)),'\circ'];['T_{0}=',num2str(T_source(2)),'\circ'];['T_{0}=',num2str(T_source(3)),'\circ'];['T_{0}=',num2str(T_source(4)),'\circ']}...
    ,'Location','SouthEast');
% legend([h5 h1 h2 h3 h4],...
%     {['Observations (MD08)']; ['T_{0}=',num2str(T_source(1)),'\circ'];['T_{0}=',num2str(T_source(2)),'\circ'];['T_{0}=',num2str(T_source(3)),'\circ'];['T_{0}=',num2str(T_source(4)),'\circ']}...
%     ,'Location','SouthEast');

text(-58,0,'\bfb)')
% xlim([-75 10])
% ylim([-500 -50])

% [newcmap]=flipud(cbrewer('seq', 'YlOrRd',length(T_source)));
% newcmap=flipud(newcmap);
% [T_surf] = Ts_to_Tc(T_site,1);
% T_inversion = (T_site + 2)./0.67;%inversion temp Jouzel and Merlivat 84 (intercept estmated from figure 8, MD 08)
% T_cond=T_site;
% 
% fig('units','inches','width',13,'height',6,'font','Helvetica','fontsize',16,'border','on');
% subplot(121)
% hold on
% plot(T_MD,d18O_MD,'o')
% ylabel('\delta^{18}O')
% xlabel('Surface Temperature (\circC)')
% 
% plot(T_surf,d18O_site(3,:),'Color',newcmap(3,:),'linewidth',2.5)
% plot(T_inversion,d18O_site(3,:),'Color',newcmap(4,:),'linewidth',2.5)
% plot(T_cond,d18O_site(3,:),'Color','k','linewidth',2.5)
% 
% subplot(122)
% hold on
% plot(T_MD,dD_MD,'o')
% ylabel('\deltaD')
% xlabel('Surface Temperature (\circC)')
% 
% plot(T_surf,dD_site(3,:),'Color',newcmap(3,:),'linewidth',2.5)
% plot(T_inversion,dD_site(3,:),'Color',newcmap(4,:),'linewidth',2.5)
% plot(T_cond,dD_site(3,:),'Color','k','linewidth',2.5)

[T_surf] = Ts_to_Tc(T_site,1);
T_inversion = (T_site + 2)./0.67;%inversion temp Jouzel and Merlivat 84 (intercept estmated from figure 8, MD 08)
T_cond=T_site;
T_test=(T_site + 5)./0.75;

fig('units','inches','width',13,'height',6,'font','Helvetica','fontsize',16,'border','on');
subplot(121)
hold on
ylabel(['\delta^{18}O (',char(8240),')'])
xlabel('Surface Temperature (\circC)')
% fill([T_inversion  T_inversion(end) fliplr(T_cond) T_cond(1)], [d18O_site(3,:) d18O_site(3,end)  fliplr(d18O_site(3,:)) d18O_site(3,1)],newcmap(2,:))
% fill([T_inversion  T_inversion(end) fliplr(T_cond) T_cond(1)], [d18O_site(3,:) d18O_site(3,end)  fliplr(d18O_site(3,:)) d18O_site(3,1)],[0.8 0.8 0.8])
h2=plot(T_cond,d18O_site(2,:),'Color','k','linewidth',2.5);
h1=plot(T_MD,d18O_MD,'o','Color',[0    0.4470    0.7410]);
% scatter(T_MD,d18O_MD,'filled','MarkerFaceColor',[0    0.4470    0.7410])

% plot(T_inversion,d18O_site(2,:),'Color',newcmap(4,:),'linewidth',2.5)
h3=plot(T_inversion,d18O_site(2,:),'--k','linewidth',2.5);
h4=plot(T_surf,d18O_site(2,:),'Color','c','linewidth',2.5);
xlim([-75 10])
ylim([-65 0])
% plot(T_test,d18O_site(2,:),'Color','c','linewidth',2.5)
l1=legend([h1 h2 h3 h4],{'observations (MD08)','T_{c}=T{s}','T_{c}=T{inv}','T_{c}=model'});
set(l1,'Position',[    0.3478    0.1921    0.1330    0.2164])
text(-70,0,'\bfa)')


subplot(122)
hold on
ylabel(['\deltaD (',char(8240),')'])
xlabel('Surface Temperature (\circC)')
% fill([T_inversion  T_inversion(end) fliplr(T_cond) T_cond(1)], [dD_site(3,:) dD_site(3,end)  fliplr(dD_site(3,:)) dD_site(3,1)],newcmap(2,:))
% fill([T_inversion  T_inversion(end) fliplr(T_cond) T_cond(1)], [dD_site(3,:) dD_site(3,end)  fliplr(dD_site(3,:)) dD_site(3,1)],[0.8 0.8 0.8])
plot(T_cond,dD_site(2,:),'Color','k','linewidth',2.5)
plot(T_MD,dD_MD,'o','Color',[0    0.4470    0.7410])
% plot(T_inversion,dD_site(2,:),'Color',newcmap(4,:),'linewidth',2.5)
plot(T_inversion,dD_site(2,:),'--k','linewidth',2.5)
plot(T_surf,dD_site(2,:),'Color','c','linewidth',2.5)
xlim([-75 10])
ylim([-500 -50])
% plot(T_test,dD_site(2,:),'Color','c','linewidth',2.5)
text(-70,-50,'\bfb)')



[newcmap]=flipud(cbrewer_BRM('seq', 'YlOrRd',length(T_source)+1));
newcmap=flipud(newcmap);
fig('units','inches','width',13,'height',6,'font','Helvetica','fontsize',16,'border','on');
subplot(121)
hold on
% h1=scatter(T_c,d18O_MD,'filled');

ylabel(['\delta^{18}O (',char(8240),')'])
xlabel('Surface Temperature (\circC)')

subplot(122)
hold on
% scatter(T_c,dD_MD,'filled')

ylabel(['\deltaD (',char(8240),')'])
xlabel('Surface Temperature (\circC)')
for i = 1:length(T_source)
subplot(121)
hold on
plot(T_surf,d18O_site(i,:),'Color',newcmap(i+1,:),'linewidth',2.5)
subplot(122)
hold on
plot(T_surf,dD_site(i,:),'Color',newcmap(i+1,:),'linewidth',2.5)
end
subplot(121)
plot(T_MD,d18O_MD,'.k');

subplot(122)
h5=plot(T_MD,dD_MD,'.k');

subplot(121)
% legend([h5],{'Observations (MD08)'},'Location','SouthEast')
xlim([-75 10])
ylim([-65 0])
text(-70,0,'\bfa)')

subplot(122)
h1=plot(0,0,'Color',newcmap(2,:),'linewidth',2.5);
h2=plot(0,0,'Color',newcmap(3,:),'linewidth',2.5);
h3=plot(0,0,'Color',newcmap(4,:),'linewidth',2.5);
h4=plot(0,0,'Color',newcmap(5,:),'linewidth',2.5);
% legend([h1 h2 h3 h4],...
%     {['T_{0}=',num2str(T_source(1)),'\circ'];['T_{0}=',num2str(T_source(2)),'\circ'];['T_{0}=',num2str(T_source(3)),'\circ'];['T_{0}=',num2str(T_source(4)),'\circ']}...
%     ,'Location','SouthEast');
legend([h5 h1 h2 h3 h4],...
    {['Observations (MD08)']; ['T_{0}=',num2str(T_source(1)),'\circ'];['T_{0}=',num2str(T_source(2)),'\circ'];['T_{0}=',num2str(T_source(3)),'\circ'];['T_{0}=',num2str(T_source(4)),'\circ']}...
    ,'Location','SouthEast');

xlim([-75 10])
ylim([-500 -50])
text(-70,-50,'\bfb)')
%%

%%

fig('units','inches','width',10,'height',10,'font','Helvetica','fontsize',16,'border','on');
subplot('Position',[.12 .55 .85 .4])
hold on
% plot(GNIP_d18O_ln_monthly,GNIP_dln_monthly,'.','color',[0.8 0.8 0.8])


    scatter(d18Oln_site(:),dlnU_site(:),S,T_source_grid(:),'LineWidth',LW)
%     contour(d18Oln_site,dlnU_site,T_source_grid)
    colormap(cmap_Tsource_cliped)
    colorbar

scatter(GNIP_d18O_ln,GNIP_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(MD08_d18O_ln,MD08_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(DomeA_d18O_ln,DomeA_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(Dahe_d18O_ln,Dahe_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(wais_snowpit_d18O_ln,wais_snowpit_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(modern_comp_d18O_ln,modern_comp_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(JIRP_d18O_ln,JIRP_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(wais_shallow_cores_d18O_ln,wais_shallow_cores_dlnU,SS,'filled','MarkerFaceColor','k')



xlabel(['\delta''','^{18}O (',char(8240),')'])
ylabel(['d_{ln} (',char(8240),')'])
xlim([-70 0])
ylim([-15 50])
text(-65,50,'\bfa)')

subplot('Position',[.12 .12 .35 .3])
hold on
 scatter(d18Oln_site(:),dDln_site(:),S,T_source_grid(:),'LineWidth',LW)
    colormap(cmap_Tsource_cliped)
%     colorbar
% plot(GNIP_d18O_ln_monthly,GNIP_dD_ln_monthly,'.','color',[0.8 0.8 0.8])

scatter(GNIP_d18O_ln,GNIP_dD_ln,SS,'filled','MarkerFaceColor','k')
scatter(MD08_d18O_ln,MD08_dD_ln,SS,'filled','MarkerFaceColor','k')
scatter(DomeA_d18O_ln,DomeA_dD_ln,SS,'filled','MarkerFaceColor','k')
scatter(Dahe_d18O_ln,Dahe_dD_ln,SS,'filled','MarkerFaceColor','k')
scatter(wais_snowpit_d18O_ln,wais_snowpit_dD_ln,SS,'filled','MarkerFaceColor','k')
scatter(modern_comp_d18O_ln,modern_comp_dD_ln,SS,'filled','MarkerFaceColor','k')
scatter(JIRP_d18O_ln,JIRP_dD_ln,SS,'filled','MarkerFaceColor','k')
scatter(wais_shallow_cores_d18O_ln,wais_shallow_cores_dD_ln,SS,'filled','MarkerFaceColor','k')
xlabel(['\delta''','^{18}O (',char(8240),')'])
ylabel(['\delta''','D(',char(8240),')'])
xlim([-70 0])
text(-67,200,'\bfb)')

subplot('Position',[.6 .12 .35 .3])
hold on
% plot(GNIP_d18O_ln_monthly,GNIP_dxs_monthly,'.','color',[0.8 0.8 0.8])

    scatter(d18O_site(:),dxs_site(:),S,T_source_grid(:),'LineWidth',LW)
    colormap(cmap_Tsource_cliped)
%     colorbar

scatter(GNIP_d18O,GNIP_dD-8.*GNIP_d18O,SS,'filled','MarkerFaceColor','k')
scatter(MD08_d18O,MD08_dD-8.*MD08_d18O,SS,'filled','MarkerFaceColor','k')
scatter(DomeA_d18O,DomeA_dD-8.*DomeA_d18O,SS,'filled','MarkerFaceColor','k')
scatter(Dahe_d18O,Dahe_dD-8.*Dahe_d18O,SS,'filled','MarkerFaceColor','k')
scatter(wais_snowpit_d18O,wais_snowpit_dD-8.*wais_snowpit_d18O,SS,'filled','MarkerFaceColor','k')
scatter(modern_comp_d18O,modern_comp_dD-8.*modern_comp_d18O,SS,'filled','MarkerFaceColor','k')
scatter(JIRP_d18O,JIRP_dD-8.*JIRP_d18O,SS,'filled','MarkerFaceColor','k')
scatter(wais_shallow_cores_d18O,wais_shallow_cores_dD-8.*wais_shallow_cores_d18O,SS,'filled','MarkerFaceColor','k')

xlabel(['\delta''','^{18}O (',char(8240),')'])
ylabel(['d_{xs}(',char(8240),')'])
xlim([-70 0])
ylim([-20 40])
text(-65,40,'\bfc)')


% fig('units','inches','width',7,'height',5,'font','Helvetica','fontsize',16,'border','on');
% hold on
% scatter(GNIP_d18O,GNIP_dD-8.*GNIP_d18O,SS,'filled','MarkerFaceColor','k')
% scatter(MD08_d18O,MD08_dD-8.*MD08_d18O,SS,'filled','MarkerFaceColor','k')
% scatter(DomeA_d18O,DomeA_dD-8.*DomeA_d18O,SS,'filled','MarkerFaceColor','k')
% scatter(Dahe_d18O,Dahe_dD-8.*Dahe_d18O,SS,'filled','MarkerFaceColor','k')
% scatter(wais_snowpit_d18O,wais_snowpit_dD-8.*wais_snowpit_d18O,SS,'filled','MarkerFaceColor','k')
% scatter(modern_comp_d18O,modern_comp_dD-8.*modern_comp_d18O,SS,'filled','MarkerFaceColor','k')
% scatter(JIRP_d18O,JIRP_dD-8.*JIRP_d18O,SS,'filled','MarkerFaceColor','k')
% scatter(wais_shallow_cores_d18O,wais_shallow_cores_dD-8.*wais_shallow_cores_d18O,SS,'filled','MarkerFaceColor','k')
% xlabel(['\delta^{18}O (',char(8240),')'])
% ylabel(['d_{xs}(',char(8240),')'])
% xlim([-70 0])
% ylim([-20 40])
% 
% fig('units','inches','width',7,'height',5,'font','Helvetica','fontsize',16,'border','on');
% hold on
% scatter(GNIP_d18O,GNIP_dD,SS,'filled','MarkerFaceColor','k')
% scatter(MD08_d18O,MD08_dD,SS,'filled','MarkerFaceColor','k')
% scatter(DomeA_d18O,DomeA_dD,SS,'filled','MarkerFaceColor','k')
% scatter(Dahe_d18O,Dahe_dD,SS,'filled','MarkerFaceColor','k')
% scatter(wais_snowpit_d18O,wais_snowpit_dD,SS,'filled','MarkerFaceColor','k')
% scatter(modern_comp_d18O,modern_comp_dD,SS,'filled','MarkerFaceColor','k')
% scatter(JIRP_d18O,JIRP_dD,SS,'filled','MarkerFaceColor','k')
% scatter(wais_shallow_cores_d18O,wais_shallow_cores_dD,SS,'filled','MarkerFaceColor','k')
% xlabel(['\delta^{18}O (',char(8240),')'])
% ylabel(['\delta D(',char(8240),')'])
% xlim([-70 0])

%%


fig(111,'units','inches','width',7,'height',5,'font','Helvetica','fontsize',16,'border','on');
hold on
scatter(GNIP_d18O_ln,GNIP_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(MD08_d18O_ln,MD08_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(DomeA_d18O_ln,DomeA_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(Dahe_d18O_ln,Dahe_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(wais_snowpit_d18O_ln,wais_snowpit_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(modern_comp_d18O_ln,modern_comp_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(JIRP_d18O_ln,JIRP_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(wais_shallow_cores_d18O_ln,wais_shallow_cores_dlnU,SS,'filled','MarkerFaceColor','k')

%     scatter(d18Oln_site(:),dlnU_site(:),S,[.2 0 1],'LineWidth',LW)
%     scatter(d18Oln_site(:),dlnU_site(:),S,[1 0 .2],'LineWidth',LW)
%     scatter(d18Oln_site(:),dlnU_site(:),S,[.2 .8 .2],'LineWidth',LW)
%     contour(d18Oln_site,dlnU_site,T_source_grid,'LineWidth',2)
%     h1=contour(d18Oln_site,dlnU_site,T_source_grid,'LineWidth',2,'LineColor',Colors(1,:));
%     h2=contour(d18Oln_site,dlnU_site,T_source_grid,'LineWidth',2,'LineColor',Colors(2,:));
    h3=contour(d18Oln_site,dlnU_site,T_source_grid,'LineWidth',2,'LineColor',Colors(3,:));

%     colormap(cmap_Tsource_cliped)
%     colorbar
xlabel(['\delta''','^{18}O (',char(8240),')'])
ylabel(['d_{ln} (',char(8240),')'])
xlim([-65 0])
ylim([-15 50])

% legend([h3, h1, h2],{['b=0.003' 'b=0.006' 'b=0.00525']})

%%
b_003 = 0.003;
[T_site_003, T_source_003, RH_source_003, d18O_site_003, dD_site_003, d18Oln_site_003, dDln_site_003, dxs_site_003, d17O_xs_site_003, dlnU_site_003, r_s_site_003] = simple_water_isotope_model_2020(Tsite, Tsource, RHsource, a, b_003, c, closure, reanalysis, SH, season);
b_007 = 0.007;
[T_site_007, T_source_007, RH_source_007, d18O_site_007, dD_site_007, d18Oln_site_007, dDln_site_007, dxs_site_007, d17O_xs_site_007, dlnU_site_007, r_s_site_007] = simple_water_isotope_model_2020(Tsite, Tsource, RHsource, a, b_007, c, closure, reanalysis, SH, season);



fig('units','inches','width',10,'height',10,'font','Helvetica','fontsize',16,'border','on');
subplot(311)
hold on

    scatter(d18Oln_site(:),dlnU_site(:),S,T_source_grid(:),'LineWidth',LW)
%     contour(d18Oln_site,dlnU_site,T_source_grid)
    colormap(cmap_Tsource_cliped)
    colorbar

scatter(GNIP_d18O_ln,GNIP_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(MD08_d18O_ln,MD08_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(DomeA_d18O_ln,DomeA_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(Dahe_d18O_ln,Dahe_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(wais_snowpit_d18O_ln,wais_snowpit_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(modern_comp_d18O_ln,modern_comp_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(JIRP_d18O_ln,JIRP_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(wais_shallow_cores_d18O_ln,wais_shallow_cores_dlnU,SS,'filled','MarkerFaceColor','k')

xlabel(['\delta''','^{18}O (',char(8240),')'])
ylabel(['d_{ln} (',char(8240),')'])
xlim([-70 0])
ylim([-15 50])
text(-65,50,'\bfa)')

subplot(312)
hold on

    scatter(d18Oln_site_003(:),dlnU_site_003(:),S,T_source_grid(:),'LineWidth',LW)
%     contour(d18Oln_site,dlnU_site,T_source_grid)
    colormap(cmap_Tsource_cliped)
    colorbar

scatter(GNIP_d18O_ln,GNIP_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(MD08_d18O_ln,MD08_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(DomeA_d18O_ln,DomeA_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(Dahe_d18O_ln,Dahe_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(wais_snowpit_d18O_ln,wais_snowpit_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(modern_comp_d18O_ln,modern_comp_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(JIRP_d18O_ln,JIRP_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(wais_shallow_cores_d18O_ln,wais_shallow_cores_dlnU,SS,'filled','MarkerFaceColor','k')

xlabel(['\delta''','^{18}O (',char(8240),')'])
ylabel(['d_{ln} (',char(8240),')'])
xlim([-70 0])
ylim([-15 50])
text(-65,50,'\bfa)')


subplot(313)
hold on

    scatter(d18Oln_site_007(:),dlnU_site_007(:),S,T_source_grid(:),'LineWidth',LW)
%     contour(d18Oln_site,dlnU_site,T_source_grid)
    colormap(cmap_Tsource_cliped)
    colorbar

scatter(GNIP_d18O_ln,GNIP_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(MD08_d18O_ln,MD08_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(DomeA_d18O_ln,DomeA_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(Dahe_d18O_ln,Dahe_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(wais_snowpit_d18O_ln,wais_snowpit_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(modern_comp_d18O_ln,modern_comp_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(JIRP_d18O_ln,JIRP_dlnU,SS,'filled','MarkerFaceColor','k')
scatter(wais_shallow_cores_d18O_ln,wais_shallow_cores_dlnU,SS,'filled','MarkerFaceColor','k')

xlabel(['\delta''','^{18}O (',char(8240),')'])
ylabel(['d_{ln} (',char(8240),')'])
xlim([-70 0])
ylim([-15 50])
text(-65,50,'\bfa)')

%%

fig('units','inches','width',10,'height',10,'font','Helvetica','fontsize',16,'border','on');
subplot(311)
hold on

    scatter(d18Oln_site(:),d17O_xs_site(:),S,T_source_grid(:),'LineWidth',LW)
%     contour(d18Oln_site,dlnU_site,T_source_grid)
    colormap(cmap_Tsource_cliped)
    colorbar



plot(landais_d18O_cor,landais_d17Oxs,'ko')
plot(wais_snowpit_d18O,wais_snowpit_d17Oxs,'ko')
plot(modern_comp_d18O,modern_comp_O17xs_PD,'k.')
plot(modern_comp_d18O,modern_comp_O17xs_EH,'kx')
plot(DomeA_d18O,DomeA_17O,'k.')
xlabel(['\delta^{18}O (',char(8240),')'])
ylabel(['\delta^{17}O_{excess} (per meg)'])
plot(wais_d18O,wais_d17Oxs,'k.')
plot(siple_d18O,siple_d17Oxs,'k.')
plot(taylor_d18O,taylor_d17Oxs,'k.')
plot(vostok_d18O,vostok_d17Oxs,'k.')
plot(vostok_d18O,vostok_d17Oxs_org,'ko')
plot(vostok_d18O,vostok_d17Oxs_offset,'kx')
plot(talos_d18O,talos_d17Oxs,'k.')
plot(edc_d18O,edc_d17Oxs,'k.')


subplot(312)
hold on

    scatter(d18Oln_site_003(:),d17O_xs_site_003(:),S,T_source_grid(:),'LineWidth',LW)
%     contour(d18Oln_site,dlnU_site,T_source_grid)
    colormap(cmap_Tsource_cliped)
    colorbar



plot(landais_d18O_cor,landais_d17Oxs,'ko')
plot(wais_snowpit_d18O,wais_snowpit_d17Oxs,'ko')
plot(modern_comp_d18O,modern_comp_O17xs_PD,'k.')
plot(modern_comp_d18O,modern_comp_O17xs_EH,'kx')
plot(DomeA_d18O,DomeA_17O,'k.')
xlabel(['\delta^{18}O (',char(8240),')'])
ylabel(['\delta^{17}O_{excess} (per meg)'])
plot(wais_d18O,wais_d17Oxs,'k.')
plot(siple_d18O,siple_d17Oxs,'k.')
plot(taylor_d18O,taylor_d17Oxs,'k.')
plot(vostok_d18O,vostok_d17Oxs,'k.')
plot(vostok_d18O,vostok_d17Oxs_org,'ko')
plot(vostok_d18O,vostok_d17Oxs_offset,'kx')
plot(talos_d18O,talos_d17Oxs,'k.')
plot(edc_d18O,edc_d17Oxs,'k.')



subplot(313)
hold on

    scatter(d18Oln_site_007(:),d17O_xs_site_007(:),S,T_source_grid(:),'LineWidth',LW)
%     contour(d18Oln_site,dlnU_site,T_source_grid)
    colormap(cmap_Tsource_cliped)
    colorbar



plot(landais_d18O_cor,landais_d17Oxs,'ko')
plot(wais_snowpit_d18O,wais_snowpit_d17Oxs,'ko')
plot(modern_comp_d18O,modern_comp_O17xs_PD,'k.')
plot(modern_comp_d18O,modern_comp_O17xs_EH,'kx')
plot(DomeA_d18O,DomeA_17O,'k.')
xlabel(['\delta^{18}O (',char(8240),')'])
ylabel(['\delta^{17}O_{excess} (per meg)'])
plot(wais_d18O,wais_d17Oxs,'k.')
plot(siple_d18O,siple_d17Oxs,'k.')
plot(taylor_d18O,taylor_d17Oxs,'k.')
plot(vostok_d18O,vostok_d17Oxs,'k.')
plot(vostok_d18O,vostok_d17Oxs_org,'ko')
plot(vostok_d18O,vostok_d17Oxs_offset,'kx')
plot(talos_d18O,talos_d17Oxs,'k.')
plot(edc_d18O,edc_d17Oxs,'k.')

