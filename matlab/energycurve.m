% This is for plotting the ernergy curve
% If air needs to be counted in, the argument number shall change


function energycurve(t, ke, gpe, dissipation, T, titlename)
    %figure;
    figure('position',[0,0,900,600]);
    total = ke + gpe;
    plot(t, ke./ke(1),'LineWidth',1.5);
    hold on
    plot(t, gpe./gpe(1), 'LineWidth',1.5);
    hold on
    plot(t, (ke+gpe)/(total(1)),'LineWidth',1.5);
    hold on
    dissipe = zeros(size(t));
    dissipe(1) = 0; 
    % calculate the dissipation;
    for i = 2:size(t)
        dissipe(i) = dissipe(i-1) + (t(i)-t(i-1))*dissipation(i)*T;
    end
    plot(t, (ke+gpe+dissipe)/(ke(1)+gpe(1)+dissipe(1)),'LineWidth',1.5);
    hold on
    plot(t,exp(-3.9478*0.001*t*T),'LineWidth',1.5);
    lgd = legend('ke','pe','ke+pe','ke+pe+dissipation','theoretical');
    lgd.Location = 'northwest';
    ax = gca;
    ax.XLim = [0. 16.];
    ax.YLim = [0.8 4.5];
    ax.LineWidth = 1.5;
    ax.FontSize = 12;
    ax.FontWeight = 'normal';
    set(gca,'TickLength',[0 0]);
    xlabel('t/T','FontSize',12);
    ylabel('E/E0','FontSize',12);
    title(titlename, 'FontWeight', 'normal');
end