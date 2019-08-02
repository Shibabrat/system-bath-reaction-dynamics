function [x,t,phi_tf,PHI] = stateTransitMat_deleonberne(x0,tf, ...
                            OPTIONS,parameters,fixed_step) 

% function [x,t,phi_tf,PHI] =
% stateTransitionMatrix_boatPR(x0,tf,R,OPTIONS,fixed_step)
%
% Gets state transition matrix, phi_tf = PHI(0,tf), and the trajectory 
% (x,t) for a length of time, tf.  
% 
% In particular, for periodic solutions of % period tf=T, one can obtain 
% the monodromy matrix, PHI(0,T).
%
% Shane Ross (revised 2009-July-07)
% Shibabrat Naik (modified 2019-Mar-23)

    N = length(x0) ;  % <-- N=4 for boat pitch-roll model

    MODEL = 'varEqns_deleonberne' ; % <-- need to change this file accordingly

    if nargin < 5,
        fixed_step = 0 ;
        if nargin < 4,
%             OPTIONS = odeset('RelTol',3e-10,'AbsTol',1e-10);  % lower accuracy
            OPTIONS = odeset('RelTol',3e-14,'AbsTol',1e-14); % high accuracy
        end
    end

    tf = tf(end);
    if fixed_step == 0,
        TSPAN = [ 0  tf ]; 
    else
        TSPAN = 0:tf/(fixed_step-1):tf ;
    end

    PHI_0(1:N^2) 	   = reshape(eye(N),N^2,1); % initial condition for state transition matrix
    PHI_0(1+N^2:N+N^2) = x0;                    % initial condition for trajectory

%     [t,PHI] = ode113(MODEL,TSPAN,PHI_0,OPTIONS);	% integration
    [t,PHI] = ode113(@(t, y)varEqns_deleonberne(t, y, parameters), ...
                            TSPAN,PHI_0,OPTIONS);	% integration

    x = PHI(:,1+N^2:N+N^2); 		   % trajectory from time 0 to tf

    phi_tf = reshape(PHI(length(t),1:N^2),N,N); % state transition matrix, PHI(O,tf)

end