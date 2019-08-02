function [xDot,isterminal,direction] = unstable_intersect_sosatyw(t,x) 

% Function to catch event of crossing y = y_w crossing of a trajectory
% 

    y_w = 1/sqrt(2);
    
    xDot = x(2);
    % The value that we want to be zero
    if x(2) > 0 && x(4) < 0 % unstable manifold's branch in the top well
        xDot = abs(x(2)) - y_w;
    elseif x(2) < 0 && x(4) > 0 % unstable manifold's branch in the bottom well
        xDot = abs(x(2)) - y_w;
    end
    
    if abs(t) > 1e-2 
        isterminal = 1; % Halt integration 
    else
        isterminal = 0; % don't terminate within a short time
    end
    
    % The zero can be approached from either direction
    direction = -1; %unstable manifold's first crossing of sos in the top well 

end