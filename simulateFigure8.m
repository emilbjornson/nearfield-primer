%%This Matlab script generates Figure 8 in the paper:
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

%Element spacing in fraction of wavelengths
scalefactor = 1/4;

%Diagonal of each RIS element
D_element = scalefactor*lambda;

%Number of RIS elements per dimension
Ndim = 100;

%Total number of RIS elements
N = Ndim^2;

%Compute the Fraunhofer distance
fraunhoferDistanceElement = 2*D_element^2/lambda;

%Björnson distance in number of Fraunhofer distance of a single element
distance_B = (2*D_element*Ndim)/(fraunhoferDistanceElement);

%Define the range of focusing points along the z-axis
relativeRange = [distance_B, N/10, N/5];



%Define the range of points along the horizontal axis
relativeRangeX = sort([ linspace(0.001,0.1,10), -linspace(0.001,0.1,10),  -logspace(log10(50),-1,101) 0 logspace(-1,log10(50),101), 0.886*relativeRange*sqrt(lambda/fraunhoferDistanceElement/N), -0.886*relativeRange*sqrt(lambda/fraunhoferDistanceElement/N) ]);

%Determine the range of distances to be considered
xRange = relativeRangeX*fraunhoferDistanceElement;
zRange = relativeRange*fraunhoferDistanceElement;

%Incident angle
theta_i = pi/6;

%Prepare to save results
channelVectors = zeros(N,length(zRange),length(xRange));

%Diagonal of entire surface
D_RIS = Ndim*D_element;


%Center points of the elements in x and y dimensions
gridPoints = (-(Ndim-1)/2:1:(Ndim-1)/2)*D_element/sqrt(2);

%Go through all distances
for mx = 1:length(xRange)
    
    for m = 1:length(zRange)
        disp(['Distance x ' num2str(mx) ' out of ' num2str(length(xRange))]);
        disp(['Distance z ' num2str(m) ' out of ' num2str(length(zRange))]);
        
        %Extract the z-distance from RIS to receiver
        z = zRange(m);
        
        %Extract the x-distance from RIS to receiver
        x_r = xRange(mx);
        
        %Define the integrand E_r(x,y) in (31) and its absolute value square
        %(we have removed E_i/(4*pi*rho)*exp(-1j*2*pi*rho/lambda)
        %since it will cancel out)
        E_fun = @(x,y) sqrt((z.*((x-x_r).^2+z.^2))./(((x-x_r).^2+y.^2+z.^2).^(5/2))).*exp(-1j*(2*pi).*(sqrt((x-x_r).^2+y.^2+z.^2)+sin(theta_i)*x)/lambda);
        E2_fun = @(x,y) (z.*((x-x_r).^2+z.^2))./(((x-x_r).^2+y.^2+z.^2).^(5/2));
        
        %Compute all the terms in (30)
        numerator_exact = zeros(Ndim,Ndim);
        denominator_exact =  integral2(E2_fun, -D_element/sqrt(8),D_element/sqrt(8),-D_element/sqrt(8),D_element/sqrt(8));
        
        
        for xdim = 1:Ndim
            
            for ydim = 1:Ndim
                
                numerator_exact(xdim,ydim) = integral2(E_fun, gridPoints(xdim)-D_element/sqrt(8),gridPoints(xdim)+D_element/sqrt(8),gridPoints(ydim)-D_element/sqrt(8),gridPoints(ydim)+D_element/sqrt(8));
                
            end
            
        end
        
        channelVectors(:,m,mx) = sqrt((2/D_element^2)/(N^2*denominator_exact))*numerator_exact(:);
        
    end
end


XfocusIndex = find(xRange==0);

%Determine the RIS gains with different phase shifts
RISGain = zeros(length(relativeRangeX),length(relativeRange));

for n = 1:length(relativeRange)
    
    phaseshifts = exp(-1j*angle(channelVectors(:,n,XfocusIndex)));
    
    for mx = 1:length(xRange)
        RISGain(mx,n) = abs(sum(channelVectors(:,n,mx).*phaseshifts)).^2;
    end
end


set(groot,'defaultAxesTickLabelInterpreter','latex');


%% Plot the simulation results
figure;
hold on; box on; grid on;
plot(relativeRangeX, RISGain(:,1),'k-', 'Linewidth', 2);
plot(relativeRangeX, RISGain(:,2),'b--', 'Linewidth', 2);
plot(relativeRangeX, RISGain(:,3),'r-.', 'Linewidth', 2);
xlim([-50 50])
xticks([-50 -25 -10 10 25 50])
xticklabels({'$-50 d_F$', '$-25 d_F$','$-10 d_F$','$10 d_F$','$25 d_F$', '$50 d_F$'})
ylim([0 1]);
yticks([0 0.25 0.5 0.75 1]);
legend({'Focus at $z=d_B$','Focus at $z=d_{FA}/10$','Focus at $z=d_{FA}/5$'},'Interpreter','Latex','Location','Best');
xlabel('$x_r$','Interpreter','Latex');
ylabel('Normalized RIS gain','Interpreter','Latex');
set(gca,'fontsize',18);
