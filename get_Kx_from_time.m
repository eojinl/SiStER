function [BC]=get_Kx_from_time(t, PARAMS, BC, GRID) 

%% set up for blister
if t < PARAMS.t_start
   BC.Kx_bot = BC.Kx_bot_bkgd; % no slip everywhere
end

if t >= PARAMS.t_start
   BC.Kx_bot= zeros(1,floor(GRID.x/GRID.dx)+1); % free traction everywhere edit
end
