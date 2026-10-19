%% Plotting code 
%% full stress field at two timesteps
L = 5e3;
h_air = 0.5e3;
x0=BC.xsize/2;
[ice] = find((h_air-1 < y) & (y < h_air + h_ice));
patch_l = xsize/2 - L/2;
patch_r = xsize/2 + L/2;
[margin] = find((patch_l - 2*h_ice < x) & (x < patch_r + 2*h_ice));
n = 256;
cmap = interp1([0 0.5 1], [0 0 1; 1 1 1; 1 0 0], linspace(0,1,n));
fig = figure();
fig.Units = 'centimeters';
fig.Position = [0 0 40 20]; % taller now for 2 panels
fig.PaperSize = [40 30];
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
% ---------- Panel 1: time = 0 ----------
time1 = 5;
fname1 = strcat("t",num2str(time1),"_d",num2str((PARAMS.t_end-PARAMS.t_start)/(3600*24)),"_G",sprintf('%0.1g',MAT(2).G),"_ep",num2str(epsilon),"_mu",sprintf('%0.1g',eta_ice),"_etaw", sprintf('%0.1g',eta_water),"_etar", sprintf('%0.1g',eta_rock),"_p",sprintf('%0.1g', BC.p),'_ell',sprintf('%0.1g',BC.slip_patch_size/2),'_viscoelastic.mat');
data1 = load(fname1);
ax1 = nexttile(tl);
xdata1 = x(margin) - x0;
imagesc(ax1, xdata1, y(ice) - h_air, data1.sxxOLD(ice, margin));
title(ax1, strcat('(B) Horizontal normal deviatoric stress $\tau_{xx}$, t=', num2str(0), ' hrs'), Interpreter='latex', FontSize=38);
ylabel(ax1, 'depth (m)', Interpreter='latex', FontSize=40);
xline(ax1, [- L/2, L/2], ':', 'LineWidth', 2, 'HandleVisibility', 'off');
colormap(gca, flipud(brewermap(256, 'RdBu')))
clim(ax1, [-5e5 5e5]);
xlim(ax1, [xdata1(1), xdata1(end)]);
set(ax1, 'YTick', 0:200:h_ice, 'YDir', 'reverse');
set(ax1, 'XTick', (2.5e4:1e3:3.5e4)-x0);
set(ax1, 'XTickLabel', []); % <-- no x tick labels
ax1.FontSize = 35;
ax1.YLabel.FontSize = 40;
ax1.Title.FontSize = 40;
% no xlabel on this panel
% ---------- Panel 2: time = 5 ----------
time2 =20;
fname2 = strcat("t",num2str(time2),"_d",num2str((PARAMS.t_end-PARAMS.t_start)/(3600*24)),"_G",sprintf('%0.1g',MAT(2).G),"_ep",num2str(epsilon),"_mu",sprintf('%0.1g',eta_ice),"_etaw", sprintf('%0.1g',eta_water),"_etar", sprintf('%0.1g',eta_rock),"_p",sprintf('%0.1g', BC.p),'_ell',sprintf('%0.1g',BC.slip_patch_size/2),'_viscoelastic.mat');
data2 = load(fname2);
ax2 = nexttile(tl);
imagesc(ax2, x(margin)-x0, y(ice) - h_air, data2.sxxOLD(ice, margin));
title(ax2, strcat('(C) Horizontal normal deviatoric stress $\tau_{xx}$, t=', num2str(15.8), ' hrs'), Interpreter='latex', FontSize=38);
ylabel(ax2, 'depth (m)', Interpreter='latex', FontSize=40);
xlabel(ax2, 'x (m)', Interpreter='latex', FontSize=40); % xlabel only here
xline(ax2, [- L/2, L/2], ':', 'LineWidth', 2, 'HandleVisibility', 'off');
colormap(gca, flipud(brewermap(256, 'RdBu')))
clim(ax2, [-5e5 5e5]);
xlim(ax2, [xdata1(1), xdata1(end)]);
set(ax2, 'YTick', 0:200:h_ice, 'YDir', 'reverse');
set(ax2, 'XTick', (2.5e4:1e3:3.5e4)-x0); % keeps labels this time
ax2.FontSize = 35;
ax2.XLabel.FontSize = 40;
ax2.YLabel.FontSize = 40;
% ---------- Shared colorbar ----------
c = colorbar(ax2);
c.Layout.Tile = 'east';
c.Label.String = '$\tau_{xx}$ (Pa)';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 35;
c.FontSize = 35;
drawnow;

%% DAA tau_xx for layers
Nt=80;
t_cut = 3; %cutoff
skip  = 4; %take every nth file after cutoff

BC.xsize=6e4;
t_vals = [3:t_cut, t_cut+skip:skip:Nt];
n_lines = numel(t_vals);
cmap = flipud(parula(n_lines));
offset=BC.xsize/2/1e3;

[ice]=find((h_air<y) & (y<h_air+h_ice));
abs_time=zeros(1, Nt);
t_times = zeros(1, n_lines); 

fig = figure();
fig.Units    = 'centimeters';
fig.Position = [0 0 30 15];
ax = gca;
hold(ax, 'on');

for i=1:n_lines % time loop
    t=t_vals(i);
    fname=strcat("t",num2str(t),"_d",num2str((PARAMS.t_end-PARAMS.t_start)/(3600*24)),"_G",sprintf('%0.1g',MAT(2).G),"_ep",num2str(epsilon),"_mu",sprintf('%0.1g',eta_ice),"_etaw", sprintf('%0.1g',eta_water),"_etar", sprintf('%0.1g',eta_rock),"_p",sprintf('%0.1g', BC.p),'_ell',sprintf('%0.1g',BC.slip_patch_size/2),'_viscoelastic.mat');
    data = load(fname);
    time=data.time;

    sxxDepthAvg_ice = (data.sxxOLD(ice(1),:)); %surface
    plot(ax, data.X(1,:)./1e3-offset, sxxDepthAvg_ice, 'LineWidth',1.5, 'Color',cmap(i,:));
    abs_time(i)=data.time;
    t_times(i) = (data.time - PARAMS.t_start)/3600;  
end

BC.slip_patch_size=5e3;

patch_l = BC.xsize/2 - BC.slip_patch_size/2;
patch_r = BC.xsize/2 + BC.slip_patch_size/2;

s0 = 0.25*MAT(2).rho0*PARAMS.gy*tand(1)*(BC.slip_patch_size);

no_slip_x = [0:GRID.dx(1):xsize/2-BC.slip_patch_size/2];
no_slip_l = -s0*exp(-sqrt(2/3)/h_ice*(no_slip_x));
no_slip_r = s0*exp(-sqrt(2/3)/h_ice*no_slip_x);

free_slip_x = [0:GRID.dx(2):BC.slip_patch_size/2];
free_slip_l = -s0*2/BC.slip_patch_size*free_slip_x;
free_slip_r = -s0*2/BC.slip_patch_size*(-free_slip_x);

% Analytical surface stress profile along the top surface
p=1e5*exp(-(abs_time(3)-PARAMS.t_start)/BC.t_rel);
L=BC.slip_patch_size;
x_local = linspace(-L/2, L/2, 1000);
x_eval = x_local + L/2;
patch_center = (patch_l + patch_r) / 2;
tau_xx_top_analytical = (p/(3*h_ice^2))*(6*L*x_eval - 6*x_eval.^2 - L^2);

% --- colorbar in place of the per-line legend ---
colormap(ax, cmap);
caxis(ax, [min(t_times), max(t_times)]);
cb = colorbar(ax);
cb.TickLabelInterpreter = 'latex';
cb.Label.Interpreter = 'latex';
cb.Label.String = 'Time (hr)';
cb.Label.FontSize = 28;
cb.FontSize = 28;

% --- axes formatting: larger fonts throughout ---
ax.FontSize             = 16;
ax.TickLabelInterpreter = 'latex';
ax.Box                  = 'on';
ax.LineWidth            = 0.75;
ax.XMinorTick           = 'on';
ax.YMinorTick           = 'on';
ylim(ax, [-5e5, 5e5])
yline(ax, 0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.75, 'HandleVisibility', 'off')

xlabel(ax, 'Horizontal position (km)',       Interpreter='latex', FontSize=28)
ylabel(ax, 'Normal stress $\tau_{xx}$ (Pa)', Interpreter='latex', FontSize=28)
title(ax, 'Surface horizontal normal stress', Interpreter='latex', FontSize=28)

% shaded patch region
patch_l_km = patch_l / 1e3-offset;
patch_r_km = patch_r / 1e3-offset;
yl = ylim(ax);
fill(ax, [patch_l_km patch_r_km patch_r_km patch_l_km], ...
     [yl(1) yl(1) yl(2) yl(2)], ...
     [0.9 0.9 0.9], ...
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.3, ...
     'HandleVisibility', 'off')
uistack(findobj(ax,'Type','patch'), 'bottom')
hold(ax, 'off')
xlim(ax, [-15,15]);