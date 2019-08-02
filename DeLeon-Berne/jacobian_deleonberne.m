function Df = jacobian_deleonberne(eqPt, par)

     %%% par = [MASS_A MASS_B EPSILON_S D_X LAMBDA ALPHA];

%     global OMEGA_X OMEGA_Y DELTA
    
    syms x y vx vy 
    
    xe = eqPt(1);
    ye = eqPt(2);
    vxe = eqPt(3);
    vye = eqPt(4);
    
    %%% Use vector differentiation 
    f(x,y,vx,vy) = [vx/par(1); 
                    vy/par(2); 
                    2*par(4)*par(5)*exp(-par(5)*x)*(exp(-par(5)*x) - 1) + ...
                    4*par(6)*par(5)*y^2*(y^2 - 1)*exp(-par(6)*par(5)*x); 
                    -8*y*(2*y^2 - 1)*exp(-par(6)*par(5)*x)];
    
    DfMat(x,y,vx,vy) = jacobian(f,[x y vx vy]);
    
    Df = double(DfMat(xe,ye,vxe,vye))
    
end




