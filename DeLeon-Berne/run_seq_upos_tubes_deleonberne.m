% Driver script for generating a sequence of unstable periodic orbits and tube
% manifolds for a range of excess energy

global deltaE

% run_po_fam;
% 
% for deltaE = 0.125:0.125:4.875 % change this variable in the script below
%         
%     run_tube_mani;
%    
%     data_path = ['./data/UPOs-deltaE',num2str(deltaE),'/'];
%     [~,~,~] = mkdir('./',data_path);
%     [~,~,~] = movefile('*.txt',data_path);
%     
%     po_fam_file = ['x0_tp_fam_eqPt',num2str(eqNum),'.txt'];
%     [~,~,~] = copyfile([data_path,po_fam_file]);
%     
%     close all;clc
%     
% end

%%

diary on
for deltaE = [0.001]%,10:10:750]
    
    close all;
%     n_mfd_traj = 2000;
    n_mfd_traj = 50;
    stbl = 1; % -1: stable, 1: unstable manifold
    dir = -1;
    
    disp(deltaE);

    % Loading data file of the periodic orbit
    data_path = ['./data/eqPt3/deltaE-',num2str(deltaE), ...
        '/x0po_T_energyPO_eqPt3_DelE',num2str(deltaE),'.txt'];
    
    x0po = importdata(data_path);

    run_tube_mani;

    % Moving the generated text files on tube manifolds 
%     data_path = ['./data/eqPt2/deltaE-',num2str(deltaE), ...
%     '/unstable-U1m-trajs',num2str(n_mfd_traj),'/'];
%     data_path = ['./data/eqPt3/deltaE-',num2str(deltaE), ...
%         '/stable-U1m-trajs',num2str(n_mfd_traj),'/'];
    data_path = ['./data/eqPt3/deltaE-',num2str(deltaE), ...
        '/unstable-U1p-trajs',num2str(n_mfd_traj),'/'];
    
    [~,~,~] = mkdir('./',data_path);
    [~,~,~] = movefile('*.txt',data_path);
    
    
   
end
diary off





