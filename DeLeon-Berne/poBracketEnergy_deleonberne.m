function [x0po,T] = poBracketEnergy_deleonberne(energyTarget,x0podata, ...
                    po_brac_file, par) 

% POBRACKETENERGY_SHIPRP Generates a family of periodic orbits (po)
% given a pair of seed initial conditions from a data file, while targeting
% a specific periodic orbit. This is performed using a scaling factor of
% the numerical continuation step size which is used to obtain initial
% guess for the periodic orbit. The energy of target periodic orbit should
% be higher than the input periodic orbit.
% 
%   [X0PO,T] = POBRACKETENERGY_SHIPRP(ENERGYTARGET, X0PODATA)
%   computes a family of periodic orbit (X0PO and T) while targeting
%   energy, ENERGYTARGET from the starting initial condition of periodic
%   orbit, X0PODATA
%
%
%
% Shane Ross (revised 2.19.04)
% Shibabrat Naik (modified for ball rolling on the surface: 26-Dec-2015)

    global N

    % delt = guessed change in period between successive orbits in family
    delt = -1.e-9 ;   % <==== may need to be changed
%     energyTol = 1e-10; % decrease tolerance for higher target energy
    energyTol = 1e-6;

    N = 4 ; % dimension of phase space
    
    x0po(1:2,1:N)   = x0podata(end-1:end,1:N);
    T(1:2,1)        = x0podata(end-1:end,N+1);
    energyPO(1:2,1) = get_total_energy_deleonberne(x0po(1:2,:),par);
    iFam = 3;
    scaleFactor = 1.25;   %scaling the change in initial guess, usually in [1,2]
    finished = 1;
    
    while finished ~= 0 || iFam > 200, 
        
%         FAMNUM = sprintf('::poFamGet : number %d',iFam) ;
%         disp(FAMNUM) ;
        
        %change in initial guess for next step
        dx  = x0po(iFam-1,1) - x0po(iFam-2,1) ;
        dy  = x0po(iFam-1,2) - x0po(iFam-2,2) ;

        %check p.o. energy and set the initial guess with scale factor
        if energyPO(iFam-1,1) < energyTarget,
            scaleFactor = scaleFactor;
            x0po_g = [ (x0po(iFam-1,1) + scaleFactor*dx) ...
                (x0po(iFam-1,2) + scaleFactor*dy) 0 0] ;
%             x0po_g = [ (x0po(iFam-1,1) + scaleFactor*dx) ...
%                 (x0po(iFam-1,2)) 0 0] ;

            [x0po_iFam,tfpo_iFam] = get_PODiffCorr_deleonberne(x0po_g, par) ; 
            energyPO(iFam,1) = get_total_energy_deleonberne(x0po_iFam, par);
            x0po(iFam,1:N) = x0po_iFam ;
            T(iFam,1)      = 2*tfpo_iFam ;
            iFam = iFam + 1;
            
            %to improve speed of stepping, increase scale factor for very
            %close p.o., this is when target p.o. hasn't been crossed,
            if abs(dx) < 1e-4 && abs(dy) < 1e-4,
                scaleFactor = 2;
            else
                scaleFactor = scaleFactor;
            end
            
        elseif energyPO(iFam-1,1) > energyTarget,
            break;
            scaleFactor = scaleFactor*1e-2;
            x0po_g = [ (x0po(iFam-2,1) + scaleFactor*dx) ...
                (x0po(iFam-2,2) + scaleFactor*dy) 0 0] ;
%             x0po_g = [ (x0po(iFam-2,1) + scaleFactor*dx) ...
%                 (x0po(iFam-2,2)) 0 0] ;

            [x0po_iFam,tfpo_iFam] = get_PODiffCorr_deleonberne(x0po_g, par) ; 
            energyPO(iFam-1,1) = get_total_energy_deleonberne(x0po_iFam, par);
            x0po(iFam-1,1:N) = x0po_iFam ;
            T(iFam-1,1)      = 2*tfpo_iFam ;
        
%         else 
%             finished = 0; % stop the loop when you cross the target energy
        end
        
        if abs(energyTarget - energyPO(iFam-1,1)) > ...
                energyTol*energyPO(iFam-1,1),
            finished = 1;
        else
            finished = 0;
        end
                    
    end

    POENERGYERR = sprintf('Relative error in the po energy from target %e', ...
        abs(energyTarget - energyPO(iFam-1,1))/energyPO(iFam-1,1)) ;
    disp(POENERGYERR) ;
        
    dum = [x0po T energyPO] ;
%     save x0po_T_energyPO.txt -ascii -double dum
    save(po_brac_file,'dum','-ascii','-double');

end
