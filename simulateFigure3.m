%%This Matlab script generates Figure 3 in the paper:
%
%Emil Björnson, Özlem Tuğfe Demir, and Luca Sanguinetti, 
%"A Primer on Near-Field Beamforming for Arrays and Reconfigurable 
%Intelligent Surfaces,"  Asilomar Conference on Signals, Systems, and
%Computers, Virtual conference, October-November 2021.
%
%Download article: https://arxiv.org/pdf/2110.06661.pdf
%
%This is version 1.0 (Last edited: 2021-10-14)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%paper as described above.

close all
clear;

%Carrier frequency
f_c = 3e9;

%Wavelength
lambda = 3e8/f_c;

%This factor can be used to change the antenna size. We need to
%have a value above 0.6 to satisfy the condition d_F >= 1.2D
scalefactor = 2;


%Diagonal of the receive antenna
D = scalefactor*lambda;

%Compute the Fraunhofer distance 
fraunhoferDistance = 2*D^2/lambda;

%Determine the range of distances to be considered
zRange = logspace(-1,1,200)*fraunhoferDistance;


%Prepare to save results
G_exact = zeros(length(zRange),1);
G_approx = zeros(length(zRange),1);


%Go through all distances
for n = 1:length(zRange)
    
    %Extract distance to transmitter
    z = zRange(n);
    
    %Define the integrand E(x,y) in (5) and its absolute value square (we
    %have removed E_0/sqrt(4*pi) since it will cancel out)
    E_fun = @(x,y) sqrt((z.*(x.^2+z.^2))./((x.^2+y.^2+z.^2).^(5/2))).*exp(-1j*(2*pi).*sqrt((x.^2+y.^2+z.^2))/lambda);
    E2_fun = @(x,y) (z.*(x.^2+z.^2))./((x.^2+y.^2+z.^2).^(5/2));
    
    %Compute the exact antenna gain in (8) numerically
    numerator_exact = abs(integral2(E_fun, -D/sqrt(8),+D/sqrt(8),-D/sqrt(8),+D/sqrt(8))).^2;
    denominator_exact = integral2(E2_fun, -D/sqrt(8),+D/sqrt(8),-D/sqrt(8),+D/sqrt(8));
    G_exact(n) = (2/D^2)*numerator_exact/denominator_exact;

    
    %Compute the approximate antenna gain using (9)
    G_approx(n) = ((8*z)/fraunhoferDistance)^2*( fresnelc(sqrt(fraunhoferDistance/(8*z))).^2 + fresnels(sqrt(fraunhoferDistance/(8*z))).^2).^2;
    

end


set(groot,'defaultAxesTickLabelInterpreter','latex');  


%% Plot the simulation results
figure;
hold on; box on; grid on;
plot(zRange/fraunhoferDistance, G_exact,'k-', 'Linewidth', 2);
plot(zRange/fraunhoferDistance, G_approx,'r-.', 'Linewidth', 2);
set(gca,'XScale','log');
xticks([0.1 0.3 1 10])
xticklabels({'$0.1 d_F$','$0.3 d_F$','$d_F$','$10 d_F$'})
ylim([0 1]);
legend({'Exact', 'Fresnel approximation'},'Interpreter','Latex','Location','SouthEast');
xlabel('Propagation distance ($z$)','Interpreter','Latex');
ylabel('Normalized antenna gain','Interpreter','Latex');
set(gca,'fontsize',18);

