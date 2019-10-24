function [xDot,isterminal,direction] = unstable_intersect_sosatx0(t,x) 

% Function to catch event of crossing x = 0 crossing of a trajectory

    x_w = 0;
    
    xDot = x(1);
    % The value that we want to be zero
    if x(2) > 0 && x(3) < 0 % unstable manifold's branch in the top well
        xDot = abs(x(1)) - x_w;
    elseif x(2) < 0 && x(3) > 0 % unstable manifold's branch in the bottom well
        xDot = abs(x(1)) - x_w;
    end
    
    if abs(t) > 1e-2 
        isterminal = 1; % Halt integration 
    else
        isterminal = 0; % don't terminate within a short time
    end
    
    % The zero can be approached from either direction
    direction = 1; % unstable manifold's 

end