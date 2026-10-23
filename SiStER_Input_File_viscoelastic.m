% SiStER_Input_File


% DURATION OF SIMULATION AND FREQUENCY OF OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%
% Nt=20; % max number of time iterations
dt_out=1; % output files every "dt_out" iterations

% 
% DOMAIN SIZE AND GRIDDING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xsize=6e4; 
ysize=1.8e3; 
BC.xsize=xsize;

h_ice=1e3; %height of ice
h_air=0.5e3; %height of air
h_rock=0.3e3; %height of rock
epsilon=0; % change slope of bed; unit of degrees

% DEFINE SETUP
BC.slip_patch_size=5e3;

BC.h_air=h_air;
BC.h_ice=h_ice;
BC.h_rock=h_rock;
% material viscosities 
eta_air=1e5;
eta_ice=1e15;
BC.G=2e9; %1e18 for viscous; 2e9 for viscoelastic
G_ice=BC.G; 
eta_rock=1e16;
eta_water=1e8;
ice_density=920;
rock_density=1000;

% for Glen's flow
nonNewtonian=true; 
A=1.2e-25;
n_flow=3;

maxwell_t = eta_ice/G_ice;
% set pressure function 
BC.p=0; %initial value 
BC.t_rel=2e4; %relaxation time for exponential

% gridding- from 0 to GRID.x(1), grid size is GRID.dx(1)
% from GRID.x(1) to GRID.x(2), grid size is GRID.dx(1) etc...
% same for y

grid=45; % 45, 90, 180
GRID.dx(1)=1.8e3/grid*2.5;
GRID.x(1)=xsize/2-BC.slip_patch_size/2-5*h_ice;
GRID.dx(2)=GRID.dx(1)/1;
GRID.x(2)=xsize/2+BC.slip_patch_size/2+5*h_ice;
GRID.dx(3)=GRID.dx(1);
GRID.x(3)=xsize;
GRID.dy(1)=GRID.dx(1)/2;
GRID.y(1)=ysize;

% LAGRANGIAN MARKERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mquad=6; % number of markers in the smallest quadrant
Mquad_crit=3; % minimum number of markers allowed in smallest quadrant (for reseeding)

% GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nphase=4; % number of phases

% phase 1: air
GEOM(1).type=1; % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(1).top=0;
GEOM(1).bot=h_air;

% phase 2: ice
GEOM(2).type=1; % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(2).top=h_air; 
GEOM(2).bot=h_air+h_ice;

% phase 3: rock 
GEOM(3).type=1; % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(3).top=h_air+h_ice; 
GEOM(3).bot=h_air+h_ice+h_rock;

% phase 4: water 
GEOM(4).type=3; % 3 = rectangle
GEOM(4).top=h_air+h_ice;
GEOM(4).bot=h_air+h_ice+h_rock;
GEOM(4).left=xsize/2-BC.slip_patch_size/2;
GEOM(4).right=xsize/2+BC.slip_patch_size/2;


% MATERIAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% creep laws of the form: pre^(-1/n)*epsII^((1-n)/n)*exp(E/(nRT))
% harmonically averaging diffusion creep, dislocation creep 
% (and plastic creep to simulate brittle failure)

% phase 1: air
MAT(1).phase=1;
% density parameters
MAT(1).rho0=1.225;
MAT(1).alpha=0;
% thermal parameters
MAT(1).k=3;
MAT(1).cp=1000;
% elasticity 
MAT(1).G=1e18;
% diffusion creep parameters
MAT(1).pre_diff=0.5/eta_air;
MAT(1).Ediff=0;
MAT(1).ndiff=1;
% dislocation creep parameters
MAT(1).pre_disc=0.5/eta_air;
MAT(1).Edisc=0;
MAT(1).ndisc=1;
% plasticity
MAT(1).mu=0.6;
MAT(1).mumin=0.3;
MAT(1).Cmax=40e6;
MAT(1).Cmin=0.01e6;
MAT(1).ecrit=0.1;

% phase 2: ice
MAT(2).phase=2;
% density parameters
MAT(2).rho0=ice_density;
MAT(2).alpha=0;
% thermal parameters
MAT(2).k=3;
MAT(2).cp=1000;
% elasticity 
MAT(2).G=G_ice;
if nonNewtonian
    % diffusion creep parameters
    MAT(2).pre_diff=A;% Glen's flow-law parameter
    MAT(2).Ediff=0;
    MAT(2).ndiff=n_flow;
    % dislocation creep parameters
    MAT(2).pre_disc=A;% Glen's flow-law parameter
    MAT(2).Edisc=0;
    MAT(2).ndisc=n_flow;
else
    % diffusion creep parameters
    MAT(2).pre_diff=0.5/eta_ice; %newtonian case is harmonic avg of pre_diff and pre_disc
    MAT(2).Ediff=0;
    MAT(2).ndiff=1;
    % dislocation creep parameters
    MAT(2).pre_disc=0.5/eta_ice;
    MAT(2).Edisc=0;
    MAT(2).ndisc=1;
end
% plasticity
MAT(2).mu=0.6;
MAT(2).mumin=0.3;
MAT(2).Cmax=40e6;
MAT(2).Cmin=0.01e6;
MAT(2).ecrit=0.1;

% phase 3: rock
MAT(3).phase=3;
% density parameters
MAT(3).rho0=rock_density; %2700 for rock
MAT(3).alpha=0;
% thermal parameters
MAT(3).k=3;
MAT(3).cp=1000;
% elasticity 
MAT(3).G=1e18; 
% diffusion creep parameters
MAT(3).pre_diff=0.5/eta_rock; %1e16 for rock
MAT(3).Ediff=0;
MAT(3).ndiff=1;
% dislocation creep parameters
MAT(3).pre_disc=0.5/eta_rock;
MAT(3).Edisc=0;
MAT(3).ndisc=1;
% plasticity
MAT(3).mu=0.6;
MAT(3).mumin=0.3;
MAT(3).Cmax=40e6;
MAT(3).Cmin=0.01e6;
MAT(3).ecrit=0.1;

% phase 4: water
MAT(4).phase=4;
% density parameters
MAT(4).rho0=rock_density; %1000 for water
MAT(4).alpha=0;
% thermal parameters
MAT(4).k=3;
MAT(4).cp=1000;
% elasticity 
MAT(4).G=1e18; 
% diffusion creep parameters
MAT(4).pre_diff=0.5/eta_rock; %initialize as rock viscosity, update after spin-up time
MAT(4).Ediff=0;
MAT(4).ndiff=1;
% dislocation creep parameters
MAT(4).pre_disc=0.5/eta_rock; %initialize as rock viscosity, update after spin-up time
MAT(4).Edisc=0;
MAT(4).ndisc=1;
% plasticity
MAT(4).mu=0.6;
MAT(4).mumin=0.3;
MAT(4).Cmax=40e6;
MAT(4).Cmin=0.01e6;
MAT(4).ecrit=0.1;

% ADDITIONAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAMS.YNElast=1; % elasticity on (1) or off (0)
PARAMS.YNPlas=0; % plasticity on (1) or off (0)
PARAMS.tau_heal=1e12; % healing time for plasticity (s)
PARAMS.epsII_from_stress=1; % get strain rate from stresses (1, default) or from velocity field (0)
PARAMS.gx=-9.8*sind(epsilon); % gravity along x, slope for window; originally 0
PARAMS.gy=9.8*cosd(epsilon); % gravity along y, slope for window; orignally 9.8
PARAMS.fracCFL=0.5; % distance by which a marker is allowed to move over a time step, expressed as a fraction of the smallest cell size
PARAMS.R=8.314; % gas constant
PARAMS.etamax=1e25; % maximum viscosity
PARAMS.etamin=1e6; % minimum viscosity
PARAMS.Tsolve=0; % yes (1) or no (0) solve for temperature
% initial temperature profile, polynomial with depth 
% T = a0 + a1*y+a2*y^2+a3*y^3+amp*sin(2*pi*X/lam)
% (make sure it matches the BCs)
PARAMS.a0=0;
PARAMS.a1=0;
PARAMS.a2=0;
PARAMS.a3=1000/(30e3)^3;
PARAMS.amp=0; % amplitude of sinusoidal perturbation
PARAMS.lam=1; % wavelength of sinusoidal perturbation
PARAMS.ynTreset=1; % if ==1, reset T=T0 where im==1 (sticky layer)
PARAMS.T0=0;
% reference values for the constant diffusivity thermal solver
% (kappa = kref / (rhoref*cpref))
PARAMS.rhoref=MAT(2).rho0; 
PARAMS.kref=3;
PARAMS.cpref=1000;

% TOPOGRAPHY EVOLUTION (interface between rock and sticky air/water layer)
PARAMS.Ntopo_markers=1000; % number of markers in marker chain tracking topography
PARAMS.YNSurfaceProcesses=0; % surface processes (diffusion of topography) on or off
% PARAMS.topo_kappa=1e-8; % diffusivity of topography (m^2/s), original
PARAMS.topo_kappa=1e-6; % diffusivity of topography (m^2/s)

% Solver iterations
PARAMS.Npicard_min=10; % minimum number of Picard iterations per time step
PARAMS.Npicard_max=100; % maximum number of Picard iterations per time step
PARAMS.conv_crit_ResL2=1e-9;
PARAMS.pitswitch=0; % number of Picard iterations at which the solver switches to quasi-Newton

% Slip patch time 
% PARAMS.t_start = 3600*24*10; % starts on day 10
PARAMS.t_start = 3000;
PARAMS.t_end = 3600*24*15; 
% BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pressure
PARAMS.p0cell=0; % pressure in the top-left corner of the domain (anchor point)


% flow

% boundary conditions
% entries in BC correspond to
% 1/ rollers? (1=yes, 0=no)
% 2/ type of velocity normal to boundary (0=constant)
% 3/ value of normal velocity 

u0 = 365*24*3600; %normal velocity m/s
visco = 1/(2*MAT(2).pre_diff);
us=-MAT(2).rho0*9.8/(2*visco)*tand(epsilon)*h_ice^2; %anlytical surface vel no slip
% xx=[0:GRID.dx:GRID.x]; %define length vector
yy=[0:GRID.dy:GRID.y]; %define depth vector
zz=min(max(yy-500,0),h_ice+h_air+h_rock);

% non-uniform xx
xx = [];
x0 = 0;
for i = 1:length(GRID.x)
    x_end = GRID.x(i);
    x_i = x0:GRID.dx(i):x_end;
    
    if i > 1
        x_i = x_i(2:end); % avoid duplicate node
    end
    
    xx = [xx x_i];
    x0 = x_end;
end

% Top and bottom boundary conditions
BC.top=[0 3 0];
BC.bot=[0 4 0 0]; %0=no slip; 1=free slip; 2=finite slippery patch; 3=tanhx boundary condition 

% non-uniform grid 
BC.bot_profile=xx*0;
slip_left  = xsize/2 - BC.slip_patch_size/2;
slip_right = xsize/2 + BC.slip_patch_size/2;
nodes_for_slip=find(xx >= slip_left & xx <= slip_right);
w = 0.05*BC.slip_patch_size; % smoothing width (m)
BC.Kx_bot_transient = (1 - tanh((xx - slip_left)/w))/2 + (1 + tanh((xx - slip_right)/w))/2;
BC.Kx_bot_bkgd = ones(1,length(xx));
BC.Kx_bot = BC.Kx_bot_bkgd;
BC.xx = xx;
BC.pressure_tanh= 0.5*(tanh((xx - slip_left)/w) - tanh((xx - slip_right)/w));

BC.smoothing_width=w;

% % uniform grid 
% nodes_for_slip=BC.slip_patch_size/GRID.dx+1;
% slip_l=floor(GRID.x/GRID.dx/2)-(nodes_for_slip-1)/2; % left side of slip patch
% slip_r=floor(GRID.x/GRID.dx/2)+(nodes_for_slip-1)/2; % right side of slip patch
% % BC.Kx_bot_transient=(1-tanh((xx-(GRID.x/2-BC.slip_patch_size/2))/250))/2+(1+tanh((xx-(GRID.x/2+BC.slip_patch_size/2))/250))/2; %normalize tanhx function, shift, factor
% BC.Kx_bot_transient=(1-tanh((xx-(GRID.x/2-BC.slip_patch_size/2))))/2+(1+tanh((xx-(GRID.x/2+BC.slip_patch_size/2))))/2; % step function
% BC.Kx_bot_bkgd=ones(1,floor(GRID.x/GRID.dx)+1); % no slip everywhere
% BC.Kx_bot=BC.Kx_bot_bkgd; % initialize boundary condition
% BC.pressure_tanh=(-tanh((xx-(GRID.x/2+(BC.slip_patch_size/2-BC.slip_patch_size/10)))/350))/2+(tanh((xx-(GRID.x/2-(BC.slip_patch_size/2-BC.slip_patch_size/10)))/350))/2; % normalize tanhx function, shift, factor


BC.top_profile=xx*0;

% Left and right boundary conditions
BC.left=[1 1 u0]; 
BC.right=[1 1 u0]; % 0 3 u0 for open condition

if nonNewtonian
    BC.right_vx=min(us*(1^(n_flow+1)-(1-zz./h_ice).^(n_flow+1)),0); %NON-NEWTONIAN analytical bc no slip
else
    BC.right_vx=min(us*(h_ice-zz)./h_ice.*(2-(h_ice-zz)./h_ice),0); %NEWTONIAN analytical bc no slip
end
BC.left_vx=BC.right_vx; %left velocity profile

PARAMS.BalanceStickyLayer=0; % if set to 1, the code will reset the inflow 
% / outflow BCs to balance the inflow / outflow of sticky layer material,
% and rock separately, based on the position of the sticky layer / air
% interface


% thermal 

% entries in BCtherm correspond to
% 1/ type? (1=Dirichlet, 2=Neumann)
% 2/ value
BCtherm.top=[1 0];
BCtherm.bot=[1 1000];
BCtherm.left=[2 0];
BCtherm.right=[2 0];
