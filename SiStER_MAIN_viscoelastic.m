% SiStER_MAIN_viscoelastic.m
%
% Reproduces viscoelastic results of: 
% Lai, C.-Y., Rines, J. H., Coffey, N. B., Lee, E., Clerc, F., Hilditch, J. P., Gagliardini, O., & L. Stevens. (2026). From blistering to sliding: Ice-sheet stress perturbations from water injection to the bed
% eojinl <at> stanford.edu
% 
% Modified from SiStER_MAIN.m: 
%
% Simple Stokes solver with Exotic Rheologies
%
% Main routine doing initialization, time loop and outputs
%
%
% J.-A. Olive, B.Z. Klein, E. Mittelstaedt, M. Behn, G. Ito, S. Howell
% jaolive <at> ldeo.columbia.edu
% March 2011 - April 2017


% INITIALIZATION

% Input File: loads parameter values, model geometry, boundary conditions
if exist('running_from_SiStER_RUN','var')==0
    InpFil = "SiStER_Input_File_viscoelastic.m";
 
    run(InpFil);
end
% run(InpFil)

% construct grid and initialize marker / node arrays
SiStER_Initialize

%% 
% BEGIN TIME LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time=0;
bkgd_dt=5;
patch_dt=20;
bkgd_min_time=1000;
patch_min_time=((PARAMS.t_end+bkgd_min_time)-(PARAMS.t_start-bkgd_min_time))/patch_dt;
dt_m=bkgd_min_time;
time_v = [];
Nt = 80;

%%
for t=1:Nt
   
    disp(['STARTING ITERATION: ' num2str(t) ' out of ' num2str(Nt)])

    % Update time
    time=time+dt_m;

    % Transient slip patch
    [BC]=get_Kx_from_time(time, PARAMS, BC, GRID); 
    
    % Here we prepare nodal arrays to feed the Stokes solver 
    SiStER_material_props_on_nodes

    if time >= PARAMS.t_start
        min_time=750; %resolve time right after blister/slip patch turn on
        MAT(4).pre_diff=0.5/eta_water; %1e6 for water
        MAT(4).pre_disc=0.5/eta_water; %1e6 for water
    else 
        min_time=bkgd_min_time;
        BC.Kx_bot_p=BC.Kx_bot_bkgd;
    end

    %%% SOLVE STOKES WITH NON-LINEAR RHEOLOGY HERE 
    SiStER_flow_solve

    % Impose pressure
    [BC]=get_decay_from_time(time, BC, PARAMS, MAT, GRID,rho, p);
    pressure=BC.p_input;

    % GET VELOCITIES ON TOPOMARKERS FOR OUTPUT/PLOTTING AVG SURF VX
    % point A (particle pts)
    % Added to original code
    X_point_test = topo_x;
    Y_point_test = topo_y;
    [qd_test,icn_test,jcn_test] = SiStER_locate_markers_in_grid(X_point_test,Y_point_test,x,y,dx,dy);
    [VX_point_test,VY_point_test] = SiStER_get_marker_velocities(qd_test,icn_test,jcn_test,x,y,topo_x,topo_y,vx,vy,dx,dy);
    
    % GET STRAIN RATE FROM CURRENT SOLUTION
    epsIIm=SiStER_interp_shear_nodes_to_markers(epsII_s,x,y,xm,ym,icn,jcn);
    
    % USE STRAIN RATE TO UPDATE STRESSES ON MARKERS
    SiStER_update_marker_stresses;
    
    % BUILD UP PLASTIC STRAIN IN YIELDING AREAS IF PLASTICITY IS ACTIVATED
    if (PARAMS.YNPlas==1) 
        SiStER_update_ep;
    end
  
    % OUTPUT VARIABLES OF INTEREST (prior to rotation & advection)
    if (mod(t,dt_out)==0 && dt_out>0) || t==1 || t==Nt % SAVING SELECTED OUTPUT
        disp('SAVING SELECTED VARIABLES TO OUTPUT FILE')
        filename=strcat("t",num2str(t),"_d",num2str((PARAMS.t_end-PARAMS.t_start)/(3600*24)),"_G",sprintf('%0.1g',MAT(2).G),"_ep",num2str(epsilon),"_mu",sprintf('%0.1g',eta_ice),"_etaw", sprintf('%0.1g',eta_water),"_etar", sprintf('%0.1g',eta_rock),"_p",sprintf('%0.1g', BC.p),'_ell',sprintf('%0.1g',BC.slip_patch_size/2),'_nonlinear_viscoelastic.mat');

        [etam]=SiStER_interp_shear_nodes_to_markers(etas,x,y,xm,ym,icn,jcn); % to visualize viscosity on markers
        save(filename, 'X','Y','vx','vy','p','time','BC','etan','topo_x','topo_y','topo_bot_x','topo_bot_y','sxxOLD','sxyOLD','PARAMS','exx_topo','EXX','EXY','pressure')
    end

    %%
    % SET ADVECTION TIME STEP BASED ON CURRENT FLOW SOLUTION
    cfl_time=SiStER_set_timestep(dx,dy,vx,vy,PARAMS);

    [dt_m] = min(min_time, cfl_time);

    disp(['Chosen dt for next step: ', num2str(dt_m)]);
    disp(['CFL limit: ', num2str(cfl_time)]);

    %%
    % ROTATE ELASTIC STRESSES IN CURRENT FLOW FIELD
    if (PARAMS.YNElast==1) 
        SiStER_rotate_stresses;
    end
    %%
    % EVOLVE TEMPERATURE FIELD THROUGH DIFFUSION
    if PARAMS.Tsolve==1
        SiStER_thermal_update;
    end

    %%
    % MARKER ADVECTION, REMOVAL, AND ADDITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SiStER_move_remove_and_reseed_markers;
    % advect markers in current flow field
    % remove markers if necessary
    % add markers if necessary
    SiStER_update_topography_markers
    % [exx_topo]=SiStER_interp_normal_nodes_to_markers(EXX,xc,yc,topo_x,topo_y,icn,jcn);
    % [sxx_topo]=SiStER_interp_normal_nodes_to_markers(sxxOLD,xc,yc,topo_x,topo_y,icn,jcn);
    % here we do the same for the marker chain that keeps track of topography
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('---------------')
    disp(['END OF ITERATION: ' num2str(t) ' out of ' num2str(Nt) ' - SIMULATION TIME: ' num2str(time/365.25/24/3600/1000) ' kyrs.'])
    disp('--------------------------------')
    disp('--------------------------------')
    fprintf('eta = %0.3g \n',max(max(etan)));
    fprintf('vy = %0.3g \n',vy(end,187));

end 
disp('FIN')
 