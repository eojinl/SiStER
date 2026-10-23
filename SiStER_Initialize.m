% SiStER Initialize

PARAMS.Nphase = Nphase; % for convenience

% construct staggered grids
[X,Y,x,y,xc,yc,dx,dy,Nx,Ny] = SiStER_initialize_grid(xsize,ysize,GRID);

% initialize marker arrays and positions
[xm, ym] = SiStER_initialize_marker_positions(xsize,ysize,dx,dy,Mquad);

% locate markers with respect to grid
[qd,icn,jcn] = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy);

% assign marker phases
[im] = SiStER_initialize_marker_phases(Nphase,GEOM,xm,ym);
% initialize marker plastic strain (to zero) and strain rate (to one)
ep=zeros(size(xm));
epNH=ep;
epsIIm=ones(size(xm));

% initialize marker stresses
sxxm=zeros(size(xm));
sxym=sxxm;

% initialize marker index (a unique number to identify and track each marker)
idm=1:length(xm);

% initialize temperature structure on nodes
T=PARAMS.a0+PARAMS.a1*Y+PARAMS.a2*Y.^2+PARAMS.a3*Y.^3;
T=T+PARAMS.amp*sin(2*pi*X/PARAMS.lam);
if PARAMS.ynTreset==1 % reset T=T0 in top layer
    T(T<PARAMS.T0)=PARAMS.T0;
end
% pass initial nodal T to markers
[Tm]=SiStER_interp_shear_nodes_to_markers(T,x,y,xm,ym,icn,jcn);
Tm0=Tm;

% initialize nodal strain rate and other useful arrays
EXX=zeros(size(X));
EXY=zeros(size(X));
vx=zeros(size(X));
vy=zeros(size(X));
v=vx;
p=1e12*ones(size(EXX));  %initialize to be high so plasticity doesnt activate at t=1, pit=1;
etan_new=zeros(Ny,Nx);
%-------------------------------------------------------------------------
% initialize dt_m small to keep things elastic & no plasticity at t=1, G.Ito
%-------------------------------------------------------------------------
if (exist('dt_m','var')==0);
    dt_m=1e2;
end

% initialize marker chain to track base of layer 1 (sticky layer)
Ntopo=PARAMS.Ntopo_markers;
topo_x=linspace(0,xsize,Ntopo);
topo_y=GEOM(1).bot*ones(size(topo_x))+2*GRID.dy(1); % HERE for interface between ice/air
topo_bot_y=GEOM(2).bot*ones(size(topo_x))-GRID.dy(1); % HERE for bottom of ice
topo_bot_y0=topo_bot_y;
topo_bot_x=topo_x;
topo_bot_x0=topo_bot_x;
topo_marker_spacing=mean(diff(topo_x)); % initial mean spacing of topography markers
exx_topo=0*topo_x;
% figure;
% plot(topo_x,topo_y)
% title('Topography')

%% Debugging regridding code
% topoL=topo_y(1);
% topoR=topo_y(end);
% topo_x_regrid=linspace(0,xsize,Ntopo); 
% topo_y_regrid=interp1(topo_x,topo_y,topo_x_regrid(2:end-1));
% topo_y_regrid=[topoL topo_y_regrid topoR];
% plot(topo_x_regrid,topo_y_regrid)
% title('Regridded topography')
% xlim([-0.1e4,3e4])
% ylim([679.5,680.5])


