function xDot = deleonberne2dof(t, x, par)

%   [stateVarsDot,termEvent, direcEvent] = barbanis2dof(t, stateVars, flag)
%   defines the odes for a Hamiltonian system with Barbanis potential
%   and has event option for one of the phase space
%   coordinates (momenta) to cross 0 value.
%
%   [stateVarsDot,termEvent, direcEvent] = ball_rolling2dof(t, stateVars, flag)
%
 
% global OMEGA_X OMEGA_Y DELTA
% par = [MASS_A MASS_B EPSILON_S D_X LAMBDA ALPHA];

% if (nargin < 3 || isempty(flag)) %without event 

    xDot = zeros(length(x),1);

    dVdx = -2*par(4)*par(5)*exp(-par(5)*x(1))*(exp(-par(5)*x(1)) - 1) - ...
        4*par(6)*par(5)*x(2)^2*(x(2)^2 - 1)*exp(-par(6)*par(5)*x(1));
    dVdy = 8*x(2)*(2*x(2)^2 - 1)*exp(-par(6)*par(5)*x(1));
    
    xDot(1) = x(3)/par(1);
    xDot(2) = x(4)/par(2);
    xDot(3) = -dVdx; 
    xDot(4) = -dVdy;
    
    
% else
%    
%     switch lower(flag)      %with event
%         case 'events'
%             if abs(t) > 1e-2
%                 isterminal = 1;  %terminate after waiting for a short time
%             else
%                 isterminal = 0;
%             end
%             
%             direction = 0; %0: all directions of crossing
%             
%             stateVarsDot = stateVars(4);
%     
%             termEvent = isterminal;
%             direcEvent = direction;            
%     end
%     
% end

end