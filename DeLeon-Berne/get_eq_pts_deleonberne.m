function [eqPt] = get_eq_pts_deleonberne(eqNum, parameters)

%GET_EQ_PTS_BP solves the saddle center equilibrium point for a system with
%KE + PE. 
%--------------------------------------------------------------------------
%   DeLeon-Berne potential energy surface notations:
%
%           Well (stable, EQNUM = 2)    
%
%               Saddle (EQNUM=1)
%
%           Well (stable, EQNUM = 3)    
%
%--------------------------------------------------------------------------
%   
    
%     global MASS_A MASS_B EPSILON_S D_X LAMBDA ALPHA
    
    %fix the equilibrium point numbering convention here and make a
    %starting guess at the solution
    if 	eqNum == 1 
        x0 = [0; 0];   % EQNUM = 1, saddle  
    elseif 	eqNum == 2 
        eqPt = [0; 1/sqrt(2)];    % EQNUM = 2, stable
        return 
    elseif 	eqNum == 3 
        eqPt = [0; -1/sqrt(2)];    % EQNUM = 3, stable
        return
    end
    
    %%% F(xEq) = 0 at the equilibrium point, solve using in-built function
    
%    options = optimoptions('fsolve','Display','iter'); % Option to display output
%    [eqPt,fval] = fsolve(@func_vec_field_eq_pt,x0,options) % Call solver
    
    [eqPt, fval, ~] = ....
        fsolve(@(x)func_vec_field_eq_pt(x,parameters),x0, ...
        optimset("jacobian","off")); % Call solver
 
end
function F = func_vec_field_eq_pt(x,par)
    
    %%% constant parameters for the surface
%     global MASS_A MASS_B EPSILON_S D_X LAMBDA ALPHA
%     par = [MASS_A MASS_B EPSILON_S D_X LAMBDA ALPHA];

    dVdx = - ( 2*par(4)*par(5)*exp(-par(5)*x(1))*(exp(-par(5)*x(1)) - 1) + ...
        4*par(6)*par(5)*x(2)^2*(x(2)^2 - 1)*exp(-par(6)*par(5)*x(1)) );
    dVdy = 8*x(2)*(2*x(2)^2 - 1)*exp(-par(6)*par(5)*x(1));
    
    F = [ -dVdx;
         -dVdy];
    
end



