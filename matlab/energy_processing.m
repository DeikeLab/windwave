% This is for processing the energy 
% arg1 is where the batch is, arg2 is the individual case name
% For example 
% arg1 = /scratch/network/jiarongw/parameter/ 
% arg2 = m2B0Ustar1ak0.05
function [] = energy_processing(arg1, arg2)
    path = arg1;
    filename = [path, arg2, '/budgetWaterwind.dat'];

    [t,ke,gpe,dissipation] = readinenergy(filename);
    T = sqrt(0.5*pi);
    domain = 1; 
    % domain = 1 means doubled domain therefore extra processing to gpe
    if (domain)
        gpe = gpe;
    end

    energycurve(t, ke, gpe, dissipation, T,'');
    plotname = [path, 'energy/', arg2, '.jpg'];
    saveas(gcf, plotname);
