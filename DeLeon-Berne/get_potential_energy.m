function pot_energy = get_potential_energy(q,par)

%     par = [MASS_A MASS_B EPSILON_S D_X LAMBDA ALPHA];
    
    pot_energy = par(4)*( 1 - exp(-par(5)*q(:,1)) ).^2 + ...
                    4*q(:,2).^2.*(q(:,2).^2 - 1).*exp(-par(6)*par(5)*q(:,1)) + ... 
                par(3);
                
    
end