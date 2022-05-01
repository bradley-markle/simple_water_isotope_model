function [T_cond_ice_core, T_source_ice_core, rs_ice_core] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file)


    
    %% set up ice core variables
ice_core_d18O_ln=(log(1 + ice_core_d18O./1000).*1000);
ice_core_dD_ln=(log(1 + ice_core_dD./1000).*1000);
ice_core_dln=ice_core_dD_ln-(-2.85*10^(-2) .* ice_core_d18O_ln.^2 + 8.47 .* ice_core_d18O_ln);
ice_core_dxs=ice_core_dD-8*ice_core_d18O;
    
    %% load existing run


    folder='./SWIM_results/';
    load(fullfile(folder,file))
%~~~~~~~~~~~~~~~

zz=repmat(T_source, length(T_site),1)';
% soursetest = griddata(d18Oln_site,dlnU_site,zz,-40,0)
ZZ=repmat(T_site, length(T_source),1);
% sitetest = griddata(d18Oln_site,dlnU_site,ZZ,-40,0)
%~~~~~~~~~~~~~~~

nan_flags =isnan(d18Oln_site);

d18Oln_site(nan_flags) = [];
dDln_site(nan_flags) = [];
dlnU_site(nan_flags) = [];
dxs_site(nan_flags) = [];
ZZ(nan_flags) = [];
zz(nan_flags) = [];

r_s_site(nan_flags) = [];

%~~~~~~~
%weird test: to account for mixing add offsets to recontruction maps.
% d18Oln_site=d18Oln_site-0.2;
% dlnU_site=dlnU_site+4;

%~~~~~~~~~~~~~~~
if method == 1
% T_cond_ice_core=griddata(d18Oln_site,dlnU_site,ZZ,ice_core_d18O_ln,ice_core_dln);
% T_source_ice_core=griddata(d18Oln_site,dlnU_site,zz,ice_core_d18O_ln,ice_core_dln);
T_cond_ice_core=griddata(d18Oln_site,dlnU_site,ZZ,ice_core_d18O_ln,ice_core_dln,'natural');
T_source_ice_core=griddata(d18Oln_site,dlnU_site,zz,ice_core_d18O_ln,ice_core_dln,'natural');
% [T_surf_ice_core] = Ts_to_Tc(T_cond_ice_core,1);
rs_ice_core=griddata(d18Oln_site,dlnU_site,r_s_site,ice_core_d18O_ln,ice_core_dln);

elseif method == 2
T_cond_ice_core=griddata(d18Oln_site,dDln_site,ZZ,ice_core_d18O_ln,ice_core_dD_ln);
T_source_ice_core=griddata(d18Oln_site,dDln_site,zz,ice_core_d18O_ln,ice_core_dD_ln);
% T_cond_ice_core=griddata(d18Oln_site,dDln_site,ZZ,ice_core_d18O_ln,ice_core_dD_ln,'natural');
% T_source_ice_core=griddata(d18Oln_site,dDln_site,zz,ice_core_d18O_ln,ice_core_dD_ln,'natural');
% [T_surf_ice_core] = Ts_to_Tc(T_cond_ice_core,1);

rs_ice_core=griddata(d18Oln_site,dDln_site,r_s_site,ice_core_d18O_ln,ice_core_dln);

elseif method == 3
T_cond_ice_core=griddata(d18Oln_site,dxs_site,ZZ,ice_core_d18O,ice_core_dxs);
T_source_ice_core=griddata(d18Oln_site,dxs_site,zz,ice_core_d18O,ice_core_dxs);
% T_cond_ice_core=griddata(d18Oln_site,dxs_site,ZZ,ice_core_d18O,ice_core_dxs,'natural');
% T_source_ice_core=griddata(d18Oln_site,dxs_site,zz,ice_core_d18O,ice_core_dxs,'natural');
% [T_surf_ice_core] = Ts_to_Tc(T_cond_ice_core,1);

rs_ice_core=griddata(d18Oln_site,dxs_site,r_s_site,ice_core_d18O_ln,ice_core_dln);

end  
end