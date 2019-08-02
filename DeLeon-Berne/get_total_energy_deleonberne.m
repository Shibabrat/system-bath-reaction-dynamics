function e = get_total_energy_deleonberne(orbit, parameters)

%   get_total_energy_deleonberne computes the total energy of an input orbit
%   (represented as M x N with M time steps and N = 4, dimension of phase
%   space for the model) for the 2 DoF DeLeon-Berne potential.
% 
%   Orbit can be different initial conditions for the periodic orbit of
%   different energy. When trajectory is input, the output energy is mean.
%

    
    x = orbit(:,1);
    y = orbit(:,2);
    px = orbit(:,3);
    py = orbit(:,4);
    
    e = (1/(2*parameters(1)))*px.^2 + (1/(2*parameters(2)))*py.^2 + ...
            get_potential_energy([x, y],parameters); 
    
%     if length(e) > 1 
%         e = mean(e);
%     end
        
    
end
