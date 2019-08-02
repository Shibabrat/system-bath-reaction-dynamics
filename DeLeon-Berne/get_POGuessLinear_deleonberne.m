function [x0poGuess,TGuess] = ...
    get_POGuessLinear_deleonberne(eqNum,Ax,par)

%   [x0poGuess,TGuess] = get_POGuessLinear_deleonberne(eqNum,Ax) ;

%   par = [MASS_A MASS_B EPSILON_S D_X LAMBDA ALPHA];  
% 
% Uses a small displacement from the equilibrium point (in a direction 
% on the collinear point's center manifold) as a first guess for a planar 
% periodic orbit (called a Lyapunov orbit in th rest. three-body problem).
%
% The initial condition and period are to be used as a first guess for
% a differential correction routine.
%
% Output:
% x0poGuess = initial state on the periodic orbit (first guess)
%           = [ x 0  0 yvel]  , (i.e., perp. to x-axis and in the plane)
% TGuess    = period of periodic orbit (first guess)
%
% Input:
% eqNum = the number of the equilibrium point of interest
% Ax    = nondim. amplitude of periodic orbit (<< 1) 
%
% Shibabrat Naik (21-Mar-19) 
% 
% 

%     global OMEGA_X OMEGA_Y
    x0poGuess  = zeros(4,1);
    
    eqPos = get_eq_pts_deleonberne(eqNum, par);
    eqPt = [eqPos' 0 0]; % phase space location of equil. point

    % Get the eigenvalues and eigenvectors of Jacobian of ODEs at equil. point
    [Es,Eu,Ec,Vs,Vu,Vc] = eqPointEig_deleonberne(eqPt, par);

    l = abs(imag(Ec(1)));
    
    if eqNum == 1
        x0poGuess(1)  = eqPt(1) - Ax;
        x0poGuess(2)  = eqPt(2);
    end
    
    %%% This is where the linearized guess based on center manifold needs
    %%% to be entered.
%     x0poGuess(1)  = ep(1) + (-Ax); %make sure that par(1) and par(2) are equal
%     x0poGuess(2)  = ep(2) + (-Ax)*(l^2 - OMEGA_X^2)/(sqrt(2)*OMEGA_X*OMEGA_Y);
    
    
    TGuess = 2*pi/l ;
    
   
end
