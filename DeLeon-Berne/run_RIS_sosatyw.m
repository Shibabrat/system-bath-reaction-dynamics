% Obtaining iterates of the stable and unstable manifold's first
% intersection with the surface of section
% 

% Load the intersection of the manifolds with the SOS: y = y_w and p_y > 0
deltaE_vals = 0.51;
smani_data_path = ...
    './';
%     '../../data/DeLeon-Berne/reaction-fraction/fig3C2/timeofflight-30/n_mfd_traj10000/';
%     '../../data/DeLeon-Berne/reaction-fraction/stable-mani-top-well/fig3A1/';
unmani_data_path = ...
    './';
%     '../../data/DeLeon-Berne/reaction-fraction/fig3C2/timeofflight-30/n_mfd_traj10000/';
%     '../../data/DeLeon-Berne/reaction-fraction/unstable-mani-top-well/fig3A1/';
% data_path = '../../data/DeLeon-Berne/UPOs-deltaE0.01/';
% data_path = './';

% unmanineg_intersect_sos = importdata([data_path, ...
%     'xeU1_unstable_branch-1_eqPt1_DelE0.01_deleonberne.txt']);
unmanipos_intersect_sos = importdata([unmani_data_path, ...
    'xeU1_unstable_branch1_eqPt1_DelE',num2str(deltaE_vals),'_deleonberne.txt']);
% smanineg_intersect_sos = importdata([data_path, ...
%     'xeU1_stable_branch-1_eqPt1_DelE',num2str(deltaE_vals),'_deleonberne.txt']);
smanipos_intersect_sos = importdata([smani_data_path, ...
    'xeU1_stable_branch1_eqPt1_DelE',num2str(deltaE_vals),'_deleonberne.txt']);

% Integrate the points on this SOS: backward for stable and forward for
% unstable manifold
% Setting up parameters and global variables
N = 4;          % dimension of phase space
MASS_A = 8.0; MASS_B = 8.0; % De Leon, Marston (1989)
EPSILON_S = 1.0;
D_X = 10.0;

% Fig. 3-A1
% ALPHA = 0.20;
% LAMBDA = 1.00;

% Fig. 3-A2
% ALPHA = 1.00;
% LAMBDA = 1.00;

% Fig. 3-B1
% ALPHA = 1.00;
% LAMBDA = 1.30;

% Fig. 3-B2
ALPHA = 1.00;
LAMBDA = 1.5;

% Fig. 3-C1
% ALPHA = 1.00;
% LAMBDA = 2.00;

% Fig. 3-C2
% ALPHA = 2.30;
% LAMBDA = 1.95;

par = [MASS_A MASS_B EPSILON_S D_X LAMBDA ALPHA];
total_energy = EPSILON_S + deltaE_vals;

curr_intersection = smanipos_intersect_sos;

tb = -20;
tf = 0;
OPTIONS = odeset('RelTol',3e-14,'AbsTol',1e-14); % high accuracy
eventSwitch = 'on';
numTimesteps = 500;

smani_iterate = {}; % cell array to store interate of the first intersection

max_iterate = 0;    
j = 1;
while j <= max_iterate
    curr_iterate = [];
    for i = 1:size(curr_intersection,1)
        [x,t,te,xe,ie] = get_traj_sos_deleonberne(curr_intersection(i,:), ...
            tb, tf, OPTIONS, eventSwitch, numTimesteps, par);

        curr_iterate = [curr_iterate; xe(end,:)];
    end
    curr_intersection = curr_iterate;
    smani_iterate{j} = curr_iterate;
    plot(smani_iterate{j}(:,1),smani_iterate{j}(:,3),'-b')
    hold on
    j = j + 1;
end

%%

curr_intersection = unmanipos_intersect_sos;

tb = 20;
tf = 0;
unmani_iterate = {};

j = 1;
while j <= max_iterate
    curr_iterate = [];
    for i = 1:size(curr_intersection,1)
        [x,t,te,xe,ie] = get_traj_sos_deleonberne(curr_intersection(i,:), ...
            tb, tf, OPTIONS, eventSwitch, numTimesteps, par);

        curr_iterate = [curr_iterate; xe(end,:)];
    end
    curr_intersection = curr_iterate;
    unmani_iterate{j} = curr_iterate;
    
    hold on
    j = j + 1;
end



%% %

% domain = [-0.025 0.025 -0.5 0.5];
% domain = [-0.25 0.25 -5 5];
domain = [-0.5 0.5 -7 7];
fimplicit(@(x,px) par(4)*( 1 - exp(-par(5).*x) ).^2 - ...
            exp(-par(6)*par(5).*x) + par(3) + px.^2/(2*par(1)) - ...
            total_energy, domain, '-m','LineWidth',1,'MeshDensity', 1000)
    
        
% axis equal
daspect([(domain(2)-domain(1)) (domain(4)-domain(3)) 1])
xlabel('$x$','Interpreter','Latex')
ylabel('$p_x$','Interpreter','Latex')
xticks([domain(1) 0 domain(2)])
xticklabels({num2str(domain(1)),'0',num2str(domain(2))})
set(gca,'TickDir','out','TickLength',[0.02 0.02]); % The only other option is 'in'
hold on

j = 1;
% plot(unmanineg_intersect_sos(:,1), unmanineg_intersect_sos(:,3),'.m')
plot(unmanipos_intersect_sos(:,1), unmanipos_intersect_sos(:,3),'-k')
% plot(smanineg_intersect_sos(:,1), smanineg_intersect_sos(:,3),'.m')
plot(smanipos_intersect_sos(:,1), smanipos_intersect_sos(:,3),'-b')
% plot(unmani_iterate{j}(:,1),unmani_iterate{j}(:,3),'-k')
% plot(smani_iterate{j}(:,1),smani_iterate{j}(:,3),'-b')



%%




