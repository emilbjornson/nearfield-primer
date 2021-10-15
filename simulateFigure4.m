%%This Matlab script generates Figure 4 in the paper:
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

%Changing this factor will not affect the result since everything becomes
%scaled accordingly
scalefactor = 1/4;

%Diagonal of each receive antenna
D_antenna = scalefactor*lambda;

%Number of receive antennas per dimension
Ndim = 25;

%Total number of receive antennas
N = Ndim^2;

%Diagonal of entire array
D_array = Ndim*D_antenna;

%Center points of the antennas in x and y dimensions
gridPoints = (-(Ndim-1)/2:1:(Ndim-1)/2)*D_antenna/sqrt(2);

relativeRange = logspace(0,3,200);


%Compute the Fraunhofer distance
fraunhoferDistanceAntenna = 2*D_antenna^2/lambda;

%Determine the range of distances to be considered
zRange = relativeRange*fraunhoferDistanceAntenna;


%Prepare to save results
G_exact = zeros(length(relativeRange),1);
G_approx = zeros(length(relativeRange),1);
G_bound = zeros(length(relativeRange),1);

%Go through all distances
for m = 1:length(zRange)
        
    disp(['Distance ' num2str(m) ' out of ' num2str(length(zRange))]);
        
    %Extract distance to transmitter
    z = zRange(m);
        
    %Define the integrand E(x,y) in (5) and its absolute value square
    %(we have removed E_0/sqrt(4*pi) since it will cancel out)
    E_fun = @(x,y) sqrt((z.*(x.^2+z.^2))./((x.^2+y.^2+z.^2).^(5/2))).*exp(-1j*(2*pi).*sqrt((x.^2+y.^2+z.^2))/lambda);
    E2_fun = @(x,y) (z.*(x.^2+z.^2))./((x.^2+y.^2+z.^2).^(5/2));
        
    %Define the integrand described in the text on Page 4 and its absolute
    %value square (we have removed E_0/sqrt(4*pi) since it will cancel out)
    E_fun2 = @(x,y) sqrt(1./((x.^2+y.^2+z.^2))).*exp(-1j*(2*pi).*sqrt((x.^2+y.^2+z.^2))/lambda);
    E2_fun2 = @(x,y) (1./((x.^2+y.^2+z.^2)));
        
        
    numerator_exact = zeros(Ndim,Ndim);
    numerator_approx = zeros(Ndim,Ndim);
        
    denominator_exact =  integral2(E2_fun, -D_antenna/sqrt(8),D_antenna/sqrt(8),-D_antenna/sqrt(8),D_antenna/sqrt(8));
    denominator_approx =  integral2(E2_fun2, -D_antenna/sqrt(8),D_antenna/sqrt(8),-D_antenna/sqrt(8),D_antenna/sqrt(8));
        
    for xdim = 1:Ndim
            
        for ydim = 1:Ndim
                
            numerator_exact(xdim,ydim) = abs(integral2(E_fun, gridPoints(xdim)-D_antenna/sqrt(8),gridPoints(xdim)+D_antenna/sqrt(8),gridPoints(ydim)-D_antenna/sqrt(8),gridPoints(ydim)+D_antenna/sqrt(8))).^2;
                
            numerator_approx(xdim,ydim) =  abs(integral2(E_fun2, gridPoints(xdim)-D_antenna/sqrt(8),gridPoints(xdim)+D_antenna/sqrt(8),gridPoints(ydim)-D_antenna/sqrt(8),gridPoints(ydim)+D_antenna/sqrt(8))).^2;
                
        end
            
    end
        
    G_exact(m) = (2/(N*D_antenna^2))*sum(numerator_exact(:))/denominator_exact;
        
    G_approx(m) = (2/(N*D_antenna^2))*sum(numerator_approx(:))/denominator_approx;
        
    alpha_z = D_antenna^2/(8*z^2);
    numerator_bound = N*alpha_z/(2*(N*alpha_z+1)*sqrt(2*N*alpha_z+1)) + atan(N*alpha_z/sqrt(2*N*alpha_z+1));
    denominator_bound = alpha_z/(2*(alpha_z+1)*sqrt(2*alpha_z+1)) + atan(alpha_z/sqrt(2*alpha_z+1));
        
    G_bound(m) = numerator_bound/(N*denominator_bound);
    
end

set(groot,'defaultAxesTickLabelInterpreter','latex');

%% Plot the simulation results
figure;
hold on; box on; grid on;
plot(relativeRange, G_approx,'b--', 'Linewidth', 2);
plot(relativeRange, G_bound,'r-.', 'Linewidth', 2);
plot(relativeRange, G_exact,'k-', 'Linewidth', 2);
set(gca,'XScale','log');
xticks([1 10 100 1000])
xticklabels({'$d_F$','$10 d_F$','$100 d_F$','$1000 d_F$'})
ylim([0 1]);
legend({'Scalar-field approximation','Upper bound','Exact'},'Interpreter','Latex','Location','SouthEast');
xlabel('Propagation distance ($z$)','Interpreter','Latex');
ylabel('Normalized antenna array gain','Interpreter','Latex');
set(gca,'fontsize',18);
