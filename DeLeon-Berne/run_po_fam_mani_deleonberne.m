%   SCRIPT to compute periodic orbits for the 2 DoF DeLeon-Berne potential
%--------------------------------------------------------------------------
%   DeLeon-Berne potential energy surface notations:
%
%           Well (stable, EQNUM = 2)    
%
%               Saddle (EQNUM=1)
%
%           Well (stable, EQNUM = 3)    
%
%--------------------------------------------------------------------------
% Shibabrat Naik (20-March-2019)
global eqNum deltaE

% global MASS_A MASS_B EPSILON_S D_X LAMBDA ALPHA

% Setting up parameters and global variables
N = 4;          % dimension of phase space
MASS_A = 8.0; MASS_B = 8.0; % De Leon, Marston (1989)
EPSILON_S = 1.0;
D_X = 10.0;

%Uncoupled system
% ALPHA = 0.0;
% LAMBDA = 1.00;

% Fig. 3-A1
% ALPHA = 0.20;
% LAMBDA = 1.00;

% Fig. 3-B2
% ALPHA = 1.00;
% LAMBDA = 1.5;

% Fig. 3-C2
ALPHA = 2.30;
LAMBDA = 1.95;


parameters = [MASS_A MASS_B EPSILON_S D_X LAMBDA ALPHA];

eqNum = 1;  
[eqPt] = get_eq_pts_deleonberne(eqNum, parameters);

eSaddle = get_total_energy_deleonberne([eqPt',0,0], parameters); % energy of the saddle eq pt

%% 

nFam = 25; % use nFam = 10 for low energy

% first two amplitudes for continuation procedure to get p.o. family
Ax1  = 2.e-5; % initial amplitude (1 of 2) values to use: 2.e-3
Ax2  = 2*Ax1; % initial amplitude (2 of 2)
deltaE = 0.510;

tic;

%  get the initial conditions and periods for a family of periodic orbits
po_fam_file = ['x0_tp_fam_eqPt',num2str(eqNum),'_deleonberne.txt'];
[po_x0Fam,po_tpFam] = get_POFam_deleonberne(eqNum, Ax1, Ax2, ...
                            nFam, po_fam_file, parameters) ; 

poFamRuntime = toc;

x0podata = [po_x0Fam, po_tpFam];


%%

% begins with a family of periodic orbits and steps until crossing the
% initial condition with target energy 
% fileName = 'x0po_T_energy_case1_L41.txt';
% fileName = 'x0po_T.txt';


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

%%

% target specific periodic orbit
% Target PO of specific energy with high precision; does not work for the
% model 

po_target_file = ['x0po_T_energyPO_eqPt',num2str(eqNum), ...
                    '_DelE',num2str(deltaE),'_deleonberne.txt'];
                
[x0_PO, T_PO, e_PO] = poTargetEnergy_deleonberne(x0poTarget, ...
                        eTarget,po_target_file,parameters);


                    
%% Setting parameters for the globalization of manifolds

n_mfd_traj = 25;
% deltaE = 0.125;

frac = 0;

del = 1e-8;

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
tmfd = 10*TPOFam(nMed);

%% % Stable manifold

% stbl = -1; dir = -1;
stbl = -1; dir = 1;

tic;    

[xW,x0W] = get_POManiLocal_deleonberne(x0po(nMed,1:4),TPOFam(nMed),frac, ...
                                            stbl,dir,del,tmfd,n_mfd_traj,parameters);

smaniRuntime = toc

energyTube = ePOFam(nMed) ;
% title(['Energy of the tube manifold:', num2str(energyTube)]);


%% Unstable manifold

% dir = 1;
% tic;    
% 
% [xW,x0W] = get_POManiLocal_deleonberne(x0po(nMed,1:4),TPOFam(nMed),frac, ...
%                                             stbl,dir,del,tmfd,n_mfd_traj,parameters);
% 
% maniRuntime = toc
% 
% energyTube = ePOFam(nMed) ;
% title(['Energy of the tube manifold:', num2str(energyTube)]);


% stbl = 1; dir = -1;
stbl = 1; dir = 1;

tic;    

[xW,x0W] = get_POManiLocal_deleonberne(x0po(nMed,1:4),TPOFam(nMed),frac, ...
                                            stbl,dir,del,tmfd,n_mfd_traj,parameters);

unmaniRuntime = toc

energyTube = ePOFam(nMed) ;
% title(['Energy of the tube manifold:', num2str(energyTube)]);



%%






























