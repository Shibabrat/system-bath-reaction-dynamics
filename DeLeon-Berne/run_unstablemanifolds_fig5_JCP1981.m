global eqNum deltaE

% global MASS_A MASS_B EPSILON_S D_X LAMBDA ZETA

% Setting up parameters and global variables
N = 4;          % dimension of phase space
MASS_A = 8.0; MASS_B = 8.0; % De Leon, Marston (1989)
EPSILON_S = 1.0;
D_X = 10.0;

% Fig. 6
ZETA_vec = [0.10, 0.40, 0.50, 0.60, 1.20, 1.60, 2.20, 2.60, 3.00];
LAMBDA = 1.00;

deltaE_vals = 0.02;

for i = 1:length(ZETA_vec)

    parameters = [MASS_A MASS_B EPSILON_S D_X LAMBDA ZETA_vec(i)];

    eqNum = 1;  
    [eqPt] = get_eq_pts_deleonberne(eqNum, parameters);

    % energy of the saddle equilibrium point
    eSaddle = ...
        get_total_energy_deleonberne([eqPt',0,0], parameters); 

    %%% 

    nFam = 100; % use nFam = 10 for low energy

    % first two amplitudes for continuation procedure to get p.o. family
    Ax1  = 2.e-5; % initial amplitude (1 of 2) values to use: 2.e-3
    Ax2  = 2*Ax1; % initial amplitude (2 of 2)

    tic;

    %  get the initial conditions and periods for a family of periodic orbits
    po_fam_file = ['x0_tp_fam_eqPt',num2str(eqNum),'_deleonberne.txt'];
    [po_x0Fam,po_tpFam] = get_POFam_deleonberne(eqNum, Ax1, Ax2, ...
                                nFam, po_fam_file, parameters) ; 

    poFamRuntime = toc;

    x0podata = [po_x0Fam, po_tpFam];
    
    n_mfd_traj = 25;
    frac = 0;
    del = 1e-8;

    paramstarttime = tic;   
    
    % Define branch of the unstable manifold
    stbl = 1; dir = 1;
    
    deltaE = deltaE_vals;

    po_fam_file = ['x0_tp_fam_eqPt',num2str(eqNum),'_deleonberne.txt'];
    eTarget = eSaddle + deltaE; 
    fprintf('Loading the periodic orbit family from data file %s \n',po_fam_file); 
    x0podata = importdata(po_fam_file);

    po_brac_file = ['x0po_T_energyPO_eqPt',num2str(eqNum), ...
                    '_brac',num2str(deltaE),'_deleonberne.txt'];
    tic;
    
    % [x0poTarget,TTarget] = bracket_POEnergy_bp(eTarget, x0podata, po_brac_file);
    [x0poTarget,TTarget] = poBracketEnergy_deleonberne(eTarget, x0podata, ...
                            po_brac_file, parameters);
    poTarE_runtime = toc;

    save(['model_parameters_eqPt',num2str(eqNum), ...
        '_DelE',num2str(deltaE),'_deleonberne.txt'], ...
        'parameters', '-ASCII', '-double');

    %%%

    % target specific periodic orbit
    % Target PO of specific energy with high precision; does not work for the
    % model 

    po_target_file = ['x0po_T_energyPO_eqPt',num2str(eqNum), ...
                        '_DelE',num2str(deltaE),'_deleonberne.txt'];

    [x0_PO, T_PO, e_PO] = poTargetEnergy_deleonberne(x0poTarget, ...
                            eTarget,po_target_file,parameters);


    % data_path = ['./data/UPOs-deltaE',num2str(deltaE),'/x0po_T_energyPO_eqPt', ...
    %     num2str(eqNum),'_DelE',num2str(deltaE),'.txt'];
    data_path = ['./x0po_T_energyPO_eqPt', num2str(eqNum), ...
        '_DelE', num2str(deltaE), '_deleonberne.txt'];
    % data_path = ['./x0po_T_energyPO_eqPt', num2str(eqNum), '_brac', ...
    %     num2str(deltaE), '.txt']

    x0po = importdata(data_path);

    TPOFam = x0po(:,5); 
    ePOFam = x0po(:,6);
    nMed = size(x0po,1);
    tmfd = 12*TPOFam(nMed);

%     manistarttime = tic;    
    
    [xW,x0W] = get_POManiLocal_deleonberne(x0po(nMed,1:4),TPOFam(nMed), ...
        frac,stbl,dir,del,tmfd,n_mfd_traj,parameters);

%     maniendtime = toc;
    
%     maniruntime = maniendtime - manistarttime;

    energyTube = ePOFam(nMed) ;
    title(['Energy of the tube manifold:', num2str(energyTube)]);
    
    paramendtime = toc;
    
    
    % Define branch of the unstable manifold
    stbl = 1; dir = -1;
    
    deltaE = deltaE_vals;

    po_fam_file = ['x0_tp_fam_eqPt',num2str(eqNum),'_deleonberne.txt'];
    eTarget = eSaddle + deltaE; 
    fprintf('Loading the periodic orbit family from data file %s \n',po_fam_file); 
    x0podata = importdata(po_fam_file);

    po_brac_file = ['x0po_T_energyPO_eqPt',num2str(eqNum), ...
                    '_brac',num2str(deltaE),'_deleonberne.txt'];
    tic;
    
    % [x0poTarget,TTarget] = bracket_POEnergy_bp(eTarget, x0podata, po_brac_file);
    [x0poTarget,TTarget] = poBracketEnergy_deleonberne(eTarget, x0podata, ...
                            po_brac_file, parameters);
    poTarE_runtime = toc;

    save(['model_parameters_eqPt',num2str(eqNum), ...
        '_DelE',num2str(deltaE),'_deleonberne.txt'], ...
        'parameters', '-ASCII', '-double');

    %%%

    % target specific periodic orbit
    % Target PO of specific energy with high precision; does not work for the
    % model 

    po_target_file = ['x0po_T_energyPO_eqPt',num2str(eqNum), ...
                        '_DelE',num2str(deltaE),'_deleonberne.txt'];

    [x0_PO, T_PO, e_PO] = poTargetEnergy_deleonberne(x0poTarget, ...
                            eTarget,po_target_file,parameters);


    % data_path = ['./data/UPOs-deltaE',num2str(deltaE),'/x0po_T_energyPO_eqPt', ...
    %     num2str(eqNum),'_DelE',num2str(deltaE),'.txt'];
    data_path = ['./x0po_T_energyPO_eqPt', num2str(eqNum), ...
        '_DelE', num2str(deltaE), '_deleonberne.txt'];
    % data_path = ['./x0po_T_energyPO_eqPt', num2str(eqNum), '_brac', ...
    %     num2str(deltaE), '.txt']

    x0po = importdata(data_path);

    TPOFam = x0po(:,5); 
    ePOFam = x0po(:,6);
    nMed = size(x0po,1);
    tmfd = 12*TPOFam(nMed);

    
    [xW,x0W] = get_POManiLocal_deleonberne(x0po(nMed,1:4),TPOFam(nMed), ...
        frac,stbl,dir,del,tmfd,n_mfd_traj,parameters);


    energyTube = ePOFam(nMed) ;
    title(['Energy of the tube manifold:', num2str(energyTube)]);
    
    paramendtime = toc;
    
    folder_name = ['./fig5/zeta',num2str(ZETA_vec(i))];
%     mkdir('./',folder_name);
%     movefile('./*.txt',folder_name);
end



%%
