% Loading and processing the manifold intersection data 
set_plot_settings('publication')
% Setting up parameters and global variables
N = 4;          % dimension of phase space
MASS_A = 8.0; MASS_B = 8.0; % De Leon, Marston (1989)
EPSILON_S = 1.0;
D_X = 10.0;
% ALPHA = 2.30;
% LAMBDA = 1.95;
ALPHA = 1.00;
LAMBDA = 1.5;
par = [MASS_A MASS_B EPSILON_S D_X LAMBDA ALPHA];
total_energy = 1.01;

%%

domain = [-0.5 0.5 -1.5 1.5];
fimplicit(@(x,y) par(4)*( 1 - exp(-par(5)*x) ).^2 + ...
            4*y.^2.*(y.^2 - 1).*exp(-par(6)*par(5)*x) + par(3) - ...
            total_energy, domain, '-k','LineWidth',1, 'MeshDensity', 1000)
    
axis tight
daspect([domain(2)-domain(1) (domain(4)-domain(3)) 1])
xlabel('$x$')
ylabel('$y$')
set(gca,'TickDir','out','TickLength',[0.02 0.02]); % The only other option is 'in'




%%

% domain = [-1 1 -1.5 1.5];
domain = [-1.5 1.5 -4.5 4.5];
fimplicit(@(y,py) 4*y.^2.*(y.^2 - 1) + py.^2/(2*par(2)) + par(3) - ...
            total_energy, domain, '-k','LineWidth',1, ...
            'MeshDensity', 1000)
    
% axis equal
daspect([(domain(2)-domain(1)) (domain(4)-domain(3)) 1])
xlabel('$y$')
ylabel('$p_y$')
xticks([domain(1) 0 domain(2)])
xticklabels({num2str(domain(1)),'0',num2str(domain(2))})
set(gca,'TickDir','out','TickLength',[0.02 0.02]); % The only other option is 'in'





%%


















