% This is for reading in multiple velocity profiles
% First set the grid number in x and y direction and the number of snapshots 
M = 64;
N = 128;
Tnum = 15;
RE = 40000;
% viscosity and period is used for non-dimensionalize height
viscosity = 1/RE * (17.4e-6/8.9e-4) * 850;
T = sqrt(2*pi);
U = zeros(N,Tnum);
%path = '/home/deike/Desktop/basilisk/wave/wind/parameters/exp/8c_sigma/drift/ak/';
path = '~/wave/ak/ak005_m5/';

cc = colormap(bone(18));
figure('position',[0,0,900,600]);

                                                                                                                                                                                                                                                                                                                                                                                       
for i = 1:Tnum
    % Non dimensionalize by ETA
    if (i == 0)
        ETA = 1;
    else
    ETA = sqrt(2*viscosity*T*i);
    end
    % Set ETA = 1 if not using scaling, otherwise comment it
    % And also change the Y limit
    % ETA = 1;
    % Non dimensionalize by u*
    USTAR = 1/sqrt(2*pi);
    str = sprintf('direct_t%d_max.dat',i);
    filename = [path,str];
    [x,y,u,v] = import_1026(filename); 
    matrixsize = size(x);
    A = zeros(matrixsize(1),2);
    A(:,1) = y/ETA;
    A(:,2) = u/USTAR;
    A_sorted = sortrows(A);
    plot(A_sorted(:,2),A_sorted(:,1),'LineWidth', 2, 'color', cc(i+1,:));
    hold on
end

ax = gca;
ax.YLim = [0. 10];
ax.XLim = [0. 20];
set(gca,'TickLength',[0 0]);
title('Velocity Profile Development', 'FontSize', 14, 'FontWeight', 'normal');
xlabel('Velocity (U/U*)', 'FontSize', 12, 'FontWeight', 'normal');
ylabel({'Height (\zeta/\zeta_{0}) '}, 'FontSize', 12, 'FontWeight', 'normal');
lgd = legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16');
lgd.Location = 'northwest';
lgd.Box = 'off';
title(lgd, 'Time Period', 'FontSize', 10);


