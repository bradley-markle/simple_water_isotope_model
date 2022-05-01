function [d18O_ln_cor,dD_ln_cor,dln_cor, d18O_cor, dD_cor]= seawater_cor_ln(d18O,dD,iso_age,d18Osw_i)
%Dec 19th 2016, Uh oh! I may have been doing this wrong! I think I have been using the
%d18Osw from Bintanja (also got the name wrong), But per Uemura (and logic)
%I should be using the ice sheet component, only. That is once the effect
%of deep sea ocean temperature is removed!!!

% load ~/Documents/work/WAIS/Data/Bintaja.mat
load ./data/Bintaja.mat
% d18Osw=interp1(Bintaja.time.*1000,Bintaja.d18Osw,iso_age);%this is wrong!
d18Osw=interp1(Bintaja.time.*1000,Bintaja.iso_ice,iso_age);
% d18Osw=interp1([min(iso_age); Bintaja.time.*1000],[0; Bintaja.iso_ice],iso_age);%BRM edit June 2020, alow data to go past 1950

% d18Osw=d18Osw';

offset=Bintaja.iso_ice(1)- (d18Osw_i);

d18Osw=d18Osw-offset;
d18Oln_sw=(log(1 + d18Osw./1000).*1000);
d18Oln_sw_i=(log(1 + d18Osw_i./1000).*1000);
deltad18Oln_sw= d18Oln_sw - d18Oln_sw_i;


%This is the method from jouzel 2003. It is wrong.
% d18O_cor=d18O-d18Osw.*((1+d18O./1000)./(1+d18Osw./1000));
% dD_cor=dD-dDsw.*((1+dD./1000)./(1+dDsw./1000));
% dxs_cor=dD_cor-8*d18O_cor;

% from Uemura:ln(1 + ?18O)corr = ln(1 + ?18Oice) ? ln(1 + ?18OSW); ln(1 + ?D)corr = ln(1 + ?Dice) ? ln(1 + ?DSW).
%ln(?Dsw + 1) = ?2.85 × 10?2 × (ln(1 + ?18Osw))2+ 8.47 × ln(1 + ?18Osw).
dDln_sw = -2.85 * 10^-2 * d18Oln_sw.^2+ 8.47 * d18Oln_sw;
dDln_sw_i = -2.85 * 10^-2 * d18Oln_sw_i.^2+ 8.47 * d18Oln_sw_i;
deltadDln_sw= dDln_sw - dDln_sw_i;



d18O_ln_cor = (log(1 + (d18O./1000)).*1000) - deltad18Oln_sw;
dD_ln_cor = (log(1 + (dD./1000)).*1000) - deltadDln_sw;

dln_cor= dD_ln_cor- (-2.85 * 10^-2 * (d18O_ln_cor).^2 +8.47 * d18O_ln_cor);

d18O_cor=(exp(d18O_ln_cor./1000)-1).*1000;
dD_cor=(exp(dD_ln_cor./1000)-1).*1000;

% dln= (log(1 + (dD./1000)).*1000)- (-2.85 * 10^-2 * (log(1 + (d18O./1000)).*1000).^2 +8.47 * (log(1 + (d18O./1000)).*1000));
% this part uses the scaling form Bintanja put the temporal variability from
% the ice core!
%  [P,S] = polyfit(d18O(~isnan(d18Osw)),(d18Osw(~isnan(d18Osw))),1); 
%  SW_test=polyval(P,d18O);
%  dDsw=8* SW_test;
% 
%  d18O_cor=d18O- SW_test.*((1+d18O./1000)./(1+ SW_test./1000));
% dD_cor=dD-dDsw.*((1+dD./1000)./(1+dDsw./1000));
% dxs_cor=dD_cor-8*d18O_cor;
%  d18O_ln_cor = log(1 + (d18O./1000)).*1000 - log(1 +  SW_test./1000).*1000;
% dDsw_ln = -2.85 * 10^-2 * (log(1 +  SW_test./1000).*1000).^2+ 8.47 * log(1 +  SW_test./1000).*1000;
% dD_ln_cor = log(1 + (dD./1000)).*1000 - dDsw_ln;
% dln_cor= dD_ln_cor- (-2.85 * 10^-2 * (d18O_ln_cor).^2 +8.47 * d18O_ln_cor);

end

