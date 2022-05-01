function [dD_sw] = d18Osw_to_dDsw(d18O_sw)
%This function calculates the dD of seawater given a specified d18O of
%seawater, based on modern correlation between d18O_sw and dD_sw
%observations

% folder = '/Users/Bradley/Documents/work/wais_dxs_paper/WDC_dxs_2016/redefining_dxs_in_antarctica/papers_and_data/';
folder='./data/'; 
% load('/Users/Bradley/Documents/work/wais_dxs_paper/WDC_dxs_2016/redefining_dxs_in_antarctica/papers_and_data/NASA_GISS_sw_d18O_dD_short.txt')
load(fullfile(folder,'NASA_GISS_sw_d18O_dD_short.txt'))
seawater_lat=NASA_GISS_sw_d18O_dD_short(:,2);
seawater_d18O=NASA_GISS_sw_d18O_dD_short(:,5);
seawater_dD=NASA_GISS_sw_d18O_dD_short(:,6);

ind=find(seawater_dD~=-99.9 & seawater_d18O~=-99.9);
seawater_lat=seawater_lat(ind);
seawater_d18O=seawater_d18O(ind);
seawater_dD=seawater_dD(ind);


[P,S] = polyfit(seawater_d18O,seawater_dD,1);
model_d18O_sw=[-0.5:0.01:0.5];
model_dD_sw= polyval(P,model_d18O_sw);

% figure
% plot(seawater_lat,seawater_d18O,'.')
% 
% figure
% hold on
% plot(seawater_d18O,seawater_dD,'.')
% plot(model_d18O_sw,model_dD_sw)

dD_sw= polyval(P,d18O_sw);

end



%%Come back to this at a later date:
%connect the d18O_sw to latitude of moisture source:

%set up for Southern Hemisphere

% % [RHn, delta_rh0, sst0, delta_sst0, rh0] = T_RH_RHn_func(T0, SH)
% %load NCEP SST data (NOAA OI SST data set)
% SST=ncread('/Users/Bradley/Documents/work/wais_dxs_paper/WDC_dxs_2016/redefining_dxs_in_antarctica/papers_and_data/ncep_data/yycompos.wQ7zPQ7Pct.nc','sst');
% lat_SST=ncread('/Users/Bradley/Documents/work/wais_dxs_paper/WDC_dxs_2016/redefining_dxs_in_antarctica/papers_and_data/ncep_data/yycompos.wQ7zPQ7Pct.nc','lat');
% lon_SST=ncread('/Users/Bradley/Documents/work/wais_dxs_paper/WDC_dxs_2016/redefining_dxs_in_antarctica/papers_and_data/ncep_data/yycompos.wQ7zPQ7Pct.nc','lon');
% SST(SST<-2)=nan;
% [lats_SST,lons_SST] = meshgrid(lat_SST,lon_SST);
% SST_line=reshape(SST(lats_SST<0),1,[]);
% lats_line=reshape(lats_SST(lats_SST<0),1,[]);
% [P_sst_lat,S_sst_lat, mu_sst_lats] = polyfit(lats_line,SST_line,9);
% lats_model=[-90:0];
% [sst_model, delta_sst] = polyval(P_sst_lat,lats_model,S_sst_lat, mu_sst_lats);
% figure;
% hold on
% plot(lats_line,SST_line,'.')
% plot(lats_model,sst_model)

