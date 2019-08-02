function  [x0_PO, T_PO, ePO] = poTargetEnergy_deleonberne(x0po, ...
                                            energyTarget, po_target_file, ...
                                            par)
%
% poTargetEnergy_deleonberne computes the periodic orbit of target energy using
% bisection method. Using bisection method on the lower and higher energy
% values of the POs to find the PO with the target energy. Use this
% condition to integrate with event function of half-period defined by
% maximum distance from the initial point on the PO
% 
% INPUT
% x0po:     Initial conditions for the periodic orbit with the last two
% initial conditions bracketing (lower and higher than) the target energy
% par: [MASS_A MASS_B EPSILON_S D_X LAMBDA ALPHA]
% 
% OUTPUTS
% x0_PO:    Initial condition of the periodic orbit (P.O.)
% T_PO:     Time period 
% ePO:     Energy of the computed P.O.
% 

%     global eSaddle
    global y0 % shared global variable with ODE integration function
%     label_fs = 10; axis_fs = 15; % small fontsize
    label_fs = 20; axis_fs = 30; % fontsize for publications 
    
    iFam = size(x0po,1);

    energyTol = 1e-10;
    tpTol = 1e-6;
    show = 1;   % for plotting the final PO
    
% bisection method begins here
    iter = 0;
    iterMax = 200;
    a = x0po(end-1,:);
    b = x0po(end,:);
    
    fprintf('Bisection method begins \n');
    while iter < iterMax
%         dx = 0.5*(b(1) - a(1));
%         dy = 0.5*(b(2) - a(2));
        
%         c = [a(1) + dx a(2) + dy 0 0];
        c = 0.5*(a + b); % guess based on midpoint
        
        [x0po_iFam,tfpo_iFam] = get_PODiffCorr_deleonberne(c, par);
        energyPO = get_total_energy_deleonberne(x0po_iFam, par);
%         x0po(iFam,1:N) = x0po_iFam ;
%         T(iFam,1)      = 2*tfpo_iFam ;
        
        c = x0po_iFam;
        iter = iter + 1;
        
        if (abs(get_total_energy_deleonberne(c, par) - ...
                energyTarget) < energyTol) || (iter == iterMax)
            fprintf('Initial condition: %e \t %e \t %e \t %e\n', c);
            fprintf('Energy of the initial condition for PO %e\n', ...
                get_total_energy_deleonberne(c, par));
            x0_PO = c
            T_PO = 2*tfpo_iFam; 
            ePO = get_total_energy_deleonberne(c, par);
            break
        end
        
        if sign( get_total_energy_deleonberne(c, par) - energyTarget ) == ...
            sign ( get_total_energy_deleonberne(a, par) - energyTarget )
            a = c;
        else
            b = c;
        end
        fprintf('Iteration number %d, energy of PO: %f\n',iter, energyPO) ;
        
    end
    fprintf('Iterations completed: %d, error in energy: %e \n', ...
        iter, abs(get_total_energy_deleonberne(c, par) - energyTarget));
    
%     if abs(t1(end) - tfpo_iFam) < tpTol % happy with convergence of bisection
%         fprintf(['Difference in TP between diff. corr and max. distance ',...
%                 'event %e\n'], abs(t1(end) - tfpo_iFam));
%     end
    
    
    % this is the check if the initial condition is really on a P.0.
    % Integrate with maximum distance event and check if the event time is
    % within tolerance of the half-period: this approach didn't work. 
    % Integrating using the same event as used for differential correction
    RelTol = 3.e-14; AbsTol = 1.e-14; 
%     model = 'barbanis2dof';  
    tspan = [0 20]; % allow sufficient time for the half-period crossing event         
%     OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol,'Events','on'); 
    OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol, ...
        'Events',@half_period_event); 
    
    % Integration with maximum distance event for the initial condition
    % obtained using bisection above
%     y0 = c';
% %     [tt,xx,t1,xx1,i1] = ode113(model, tspan, y0, OPTIONS);
%     [tt,xx,t1,xx1,i1] = ode113(@(t, y)deleonberne2dof(t, y, par), ...
%                                 tspan,y0,OPTIONS);
%     
% %     xx1
% %     t1
% %     x0_PO = xx1(end,:);
% %         T_PO = 2*t1(end); 
% %         e_PO = get_energy_points_ball_rolling(xx1(end,:)); 
% %         e_PO = get_energy_points_ball_rolling(y0');
% %         x0po(iFam + (iter+1),1:N) = c ;         
% %         c = xx1(end,:);
% %     abs(t1(end) - tfpo_iFam)
% 
%     
%     
%     % assuming symmetry in the PO, we copy the half orbit
%     po_traj = [[xx(:,1);flip(xx(:,1))], [xx(:,2);flip(xx(:,2))], ...
%                 [xx(:,3);flip(xx(:,3))], [xx(:,4);flip(-xx(:,4))]];
%             
%     if show == 1
%         e = get_total_energy_deleonberne(xx, par);
%         fprintf('Energy of the periodic trajectory %e\n', e);
% %         set(0,'Units','normalized')
%         hFig = figure();
% %         set(gcf,'PaperPositionMode','auto')
% %         set(hFig, 'Position', [0, 0, 0.4, 0.4])
% 
%         plot3([xx(:,1);flip(xx(:,1))], [xx(:,2);flip(xx(:,2))], ...
%                 [xx(:,3);flip(-xx(:,3))],'-','Linewidth',2); 
%         hold on;
%         m = length(xx) ;
%         plot3(xx(1,1), xx(1,2),xx(1,3),'b*','Markersize', 10);
%         plot3(xx(m,1), xx(m,2),xx(m,3),'bx','Markersize', 10);
%         set(gca,'fontsize',label_fs)
% %         axis([min(x(:,1)) max(x(:,1)) min(x(:,2)) max(x(:,2))]);
% %         axis equal % destroyes the scale for successive plots
%         xlabel('$x$','interpreter','latex','fontsize',axis_fs);
%         ylabel('$y$','interpreter','latex','fontsize',axis_fs);
%         zlabel('$v_x$','interpreter','latex','fontsize',axis_fs);
%         grid on
%         box on
% %         view(-90,0)
%         title(['$\Delta E =$ ',num2str(mean(e) - par(3))], ...
%                 'interpreter','latex');
% %         leg1 = legend(['$\Delta E = $',num2str(energyTarget - eSaddle)]);
% %         set(leg1, 'interpreter','latex');
%             
%     end
   
    dum = [c 2*tfpo_iFam energyPO];
    save(po_target_file,'dum','-ascii','-double')
    
%     save('po_traj.txt','po_traj','-ascii','-double');

end





