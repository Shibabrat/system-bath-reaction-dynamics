function [x0po,t1] = get_PODiffCorr_deleonberne(x0, par)

% parameters = [MASS_A MASS_B EPSILON_S D_X LAMBDA ALPHA];
% global eSaddle
% global OMEGA_X OMEGA_Y DELTA

% set show = 1 to plot successive approximations (default=0)
show = 1 ;
% axesFontName = 'factory';
% axesFontName = 'Times New Roman';
% label_fs = 20; axis_fs = 30; % fontsize for publications 
label_fs = 10; axis_fs = 15; % fontsize for publications 
% set(0,'Defaulttextinterpreter','latex', ...
%     'DefaultAxesFontName', axesFontName, ...
%     'DefaultTextFontName', axesFontName, ...
%     'DefaultAxesFontSize', axis_fs, ...
%     'DefaultTextFontSize', label_fs, ...
%     'Defaultuicontrolfontweight','normal', ...
%     'Defaulttextfontweight','normal', ...
%     'Defaultaxesfontweight','normal');


% tolerances for integration and perpendicular crossing of x-axis
% MAXdxdot1 = 1.e-8 ; RelTol = 3.e-10; AbsTol = 1.e-10; 
% MAXdxdot1 = 1.e-12 ; RelTol = 3.e-14; AbsTol = 1.e-14; 
MAXdxdot1 = 1.e-10 ; RelTol = 3.e-14; AbsTol = 1.e-14; 

MAXattempt = 100;     	% maximum number of attempts before error is declared

% dxdot1 	   = 1;         % to start while loop
dydot1 	   = 1;         % to start while loop
attempt    = 0;         % begin counting number of attempts
% y0(attempt) = 1;

while abs(dydot1) > MAXdxdot1
    
	if attempt > MAXattempt
		ERROR = 'Maximum iterations exceeded' ;
		disp(ERROR) ;
		break
    end
    
    y0 = x0;
    % Find first half-period crossing event
%     MODEL = 'deleonberne2dof';  
    TSPAN = [0 10];        % allow sufficient time for the half-period crossing event         
    OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol, ...
        'Events',@half_period_event); 
%     [tt,xx,t1,xx1,i1] = ode113(MODEL,TSPAN,x0,OPTIONS);
    [tt,xx,t1,xx1,i1] = ode113(@(t, y)deleonberne2dof(t, y, par), ...
                                TSPAN,x0,OPTIONS);
    
    
	x1 = xx1(end,1); 
    y1 = xx1(end,2); 
	dxdot1 = xx1(end,3); 
    dydot1  = xx1(end,4); 
%     plot3(xx(:,1),xx(:,2),xx(:,3),'-r');hold on
%     plot3(x1,y1,dxdot1,'xb');
   
    
    % Compute the state transition matrix from the initial state to
	% the final state at the half-period event crossing
      
    % Events option not necessary anymore
    OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol); 
	[x,t,phi_t1,PHI] = stateTransitMat_deleonberne(x0,t1,OPTIONS,par) ;

	attempt = attempt+1 ;
    
% 	ATTEMPT = sprintf('::poDifCor : iteration %d',attempt) ;
% 	disp(ATTEMPT) ;
     
    if show == 1
        e = get_total_energy_deleonberne(x0,par);
        
        plot3(x(:,1),x(:,2),x(:,3),'.-',x(:,1),x(:,2),-x(:,3),'.-'); 
        hold on;
        m = length(x) ;
        plot3(x(1,1),x(1,2),x(1,3),'b*');
        plot3(x(m,1),x(m,2),x(m,3),'bo');
        
        
%         set(gca,'fontsize',18)
%         axis([min(x(:,1)) max(x(:,1)) min(x(:,2)) max(x(:,2))]);
%         axis equal % destroyes the scale for successive plots
%         set(gca,'fontsize',label_fs)
%         axis([min(x(:,1)) max(x(:,1)) min(x(:,2)) max(x(:,2))]);
%         axis equal % destroyes the scale for successive plots
        xlabel('$x$','interpreter','latex','fontsize',axis_fs);
        ylabel('$y$','interpreter','latex','fontsize',axis_fs);
        zlabel('$p_x$','interpreter','latex','fontsize',axis_fs);
        title(['$\Delta E$ = ',num2str(mean(e) - par(3))], ...
            'interpreter','latex','fontsize',axis_fs);
%         xlim([-15 15])
%         ylim([-15 15])
        pause(0.01) ;
        grid on
        box on

    end


%=========================================================================
% differential correction and convergence test, adjust according to
% the particular problem

    %compute acceleration values for correction term
    dVdx = -2*par(4)*par(5)*exp(-par(5)*x1)*(exp(-par(5)*x1) - 1) - ...
        4*par(6)*par(5)*y1^2*(y1^2 - 1)*exp(-par(6)*par(5)*y1);
    dVdy = 8*y1*(2*y1^2 - 1)*exp(-par(6)*par(5)*x1);
    
    vxdot1 = -dVdx;
    vydot1 = -dVdy;
    
    %correction to the initial x0
%     y0(attempt) = dxdot1/(phi_t1(3,1) - phi_t1(4,1)*(vxdot1/vydot1));	
%     x0(1) = x0(1) - y0(attempt);
    
    %correction to the initial y0
    y0(attempt) = 1/(phi_t1(4,2) - phi_t1(3,2)*vydot1*(1/vxdot1))*dydot1; 
	x0(2) = x0(2) - y0(attempt);
  
   
end

x0po=x0;
t1 = t1(end);

end