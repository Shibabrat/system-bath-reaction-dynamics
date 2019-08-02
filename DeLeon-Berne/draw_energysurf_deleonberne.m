function fs = draw_energysurf_deleonberne(z,lambda,H_val,alpha)

% @author: Víctor J. García-Garrido
% 
% Shibabrat Naik, modified on 2-Apr-2019
% added plot styling
% added epsilon_s to the energy function according to De Leon-Berne,1981.
% added mesh density option for the fimplicit3

    % plot properties
    axesFontName = 'factory';
    % axesFontName = 'Times New Roman';
    axFont = 18;
    textFont = 18;
    labelFont = 25;
    lw = 2;    
    set(0,'Defaulttextinterpreter','latex', ...
        'DefaultAxesFontName', axesFontName, ...
        'DefaultTextFontName', axesFontName, ...
        'DefaultAxesFontSize',axFont, ...
        'DefaultTextFontSize',textFont, ...
        'Defaultuicontrolfontweight','normal', ...
        'Defaulttextfontweight','normal', ...
        'Defaultaxesfontweight','normal');

    % Parameters of the De Leon-Berne Hamiltonian
    mu = 8;             % mass of isomer
    % z interaction parameter
    % lambda Morse range parameter
    D = 10;                 % Dissociation Energy
    bh = 1;                 % barrier height
    wd = 1 / sqrt(2);   % distance from barrier to well
    epsilon_s = 1;
    
    esurf = @(x,y,px) (bh ./ (wd).^4) .* y.^2 .* (y.^2 - 2 .* wd.^2) ...
        .* exp(-z .* lambda .* x) + D .* (1 - exp(-lambda .* x)).^2 + ...
        epsilon_s + (1 ./ (2 .* mu)) .* px.^2  - H_val;
    
    rgb_col = [51/255 153/255 51/255];
%     
%     xi = -0.3; xf = 0.3;
%     yi = -1.2; yf = 1.2;
%     pxi = -5.5; pxf = 5.5;
    
    xi = -1; xf = 1;
    yi = -0.5; yf = 2;
    pxi = -10; pxf = 10;

    fs = fimplicit3(esurf,[xi xf yi yf pxi pxf],...
        'EdgeColor','none','MeshDensity',100,'FaceAlpha',alpha,...
        'FaceColor',rgb_col);
  
    xlabel('$x$','FontSize',labelFont,'Interpreter','Latex');
    ylabel('$y$','FontSize',labelFont,'Interpreter','Latex');
    zlabel('$p_x$','FontSize',labelFont,'Interpreter','Latex');
    light;
    
    % fs.FaceAlpha = alpha;
    fs.AmbientStrength = 1;
    
end

