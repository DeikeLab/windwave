% calculating the linlog
RE = 40000;
viscosity = 1/RE * (17.4e-6/8.9e-4) * 850;
USTAR = 1/sqrt(2*pi);
KARMEN = 0.41;
m = 5;
z1 = m*viscosity/USTAR*5;
figure('position',[0,0,600,400]);

z = 0:0.00001:0.2;
y1 = zeros(size(z));
for i = 1:20001
    if (z(i) < z1)
        y1(i) = USTAR*USTAR/viscosity*z(i);
    else
        beta = 2*KARMEN*USTAR/viscosity*(z(i)-z1);
	    alpha = log(beta+sqrt(beta^2+1));
        tanh= (exp(alpha/2)-exp(-alpha/2))/(exp(alpha/2)+exp(-alpha/2));     
        y1(i) = m*USTAR*5 + USTAR/KARMEN*(alpha-tanh);
    end
end

plot(y1, z, 'LineWidth', 1.5);
hold on 

y2 = zeros(size(z));
UFAR = y1(20001);
delta_v = sqrt(1/RE);
delta = 2*delta_v;
y2 = UFAR - UFAR*exp(-z./delta);
plot(y2, z, 'LineWidth', 1.5);
ax = gca;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
%set(gca,'TickLength',[0 0]);
xlabel('Velocity(U/U*)', 'FontSize', 12, 'FontWeight', 'normal');
ylabel({'Height(\zeta/\zeta_{0})'}, 'FontSize', 12, 'FontWeight', 'normal');
lgd = legend('Lin-log','Exponential');
lgd.Location = 'northwest';
lgd.Box = 'off';
title(lgd, 'Velocity Profile', 'FontSize', 10);