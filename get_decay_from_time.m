function [BC]=get_decay_from_time(t, BC, PARAMS, MAT, GRID,rho,p) 
%% pressure for layers
if t < PARAMS.t_start
    BC.p_input=0;
    BC.Kx_bot_p=BC.Kx_bot_bkgd;
elseif t >= PARAMS.t_start
    p_lith=MAT(1).rho0*PARAMS.gy*BC.h_air + MAT(2).rho0*PARAMS.gy*(BC.h_ice)+MAT(4).rho0*PARAMS.gy*BC.h_rock;
    
    BC.bot_profile_p=BC.xx*0 + p_lith + BC.pressure_tanh.*BC.p.*exp(-(t-PARAMS.t_start)/BC.t_rel); %exponential pressure 
    BC.Kx_bot_p=BC.Kx_bot_transient;
    BC.p_input=BC.p*exp(-(t-PARAMS.t_start)/BC.t_rel); % exp decay for visualization purposes
end

