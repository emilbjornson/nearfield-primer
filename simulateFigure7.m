%%This Matlab script generates Figure 7 in the paper:
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

%Define the range of points along the horizontal axis
relativeRange = sort([logspace(1,5,300) 400 1000 10000]);


%Compute the Fraunhofer distance
fraunhoferDistanceElement = 2*D_element^2/lambda;

%Determine the range of distances to be considered
zRange = relativeRange*fraunhoferDistanceElement;

%Incident angle
theta_i = pi/6;


%Prepare to save results
channelVectors = zeros(N,length(relativeRange));

%Diagonal of entire surface
D_RIS = Ndim*D_element;


%Center points of the elements in x and y dimensions
gridPoints = (-(Ndim-1)/2:1:(Ndim-1)/2)*D_element/sqrt(2);



%Go through all distances
for m = 1:length(zRange)
        
    disp(['Distance ' num2str(m) ' out of ' num2str(length(zRange))]);
        
    %Extract distance from RIS to receiver
    z = zRange(m);
    
    %x-axis value of the receiver
    x_r = 0;
    
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
        
    channelVectors(:,m) = sqrt((2/D_element^2)/(N^2*denominator_exact))*numerator_exact(:);
        
    
end

%Determine the channel for the far-field receiver
E0_fun = @(x,y) exp(-1j*(2*pi).*(sin(theta_i)*x)/lambda);
numerator_exact0 = zeros(Ndim,Ndim);


for xdim = 1:Ndim
    
    for ydim = 1:Ndim
        
        numerator_exact0(xdim,ydim) = integral2(E0_fun, gridPoints(xdim)-D_element/sqrt(8),gridPoints(xdim)+D_element/sqrt(8),gridPoints(ydim)-D_element/sqrt(8),gridPoints(ydim)+D_element/sqrt(8));
        
    end
    
end

channelVectors0 = sqrt((2/D_element^2)^2/N^2)*numerator_exact0(:);

%Determine which points the RIS focuses at
focusIndex = [find(relativeRange==400) find(relativeRange==1000) find(relativeRange==10000)]; 


%Determine the RIS gains with different phase shifts
RISGain = zeros(length(relativeRange),length(focusIndex));

for n = 1:length(focusIndex)
    
    phaseshifts = exp(1j*angle(channelVectors(:,focusIndex(n))));
    
    RISGain(:,n) = abs(channelVectors'*phaseshifts).^2;
    
end

%Determine the phase-shifts for focusing at far-field
phaseshifts0 = exp(-1j*angle(channelVectors0));


%The following is for optimized RIS phase-shifts
optimizedPhaseShiftGains = (sum(abs(channelVectors),1)').^2;

RISGainFarfield = abs(sum(channelVectors.*phaseshifts0,1)').^2;
    

set(groot,'defaultAxesTickLabelInterpreter','latex');


%% Plot the simulation results
figure;
hold on; box on; grid on;
plot(relativeRange, optimizedPhaseShiftGains,'k:', 'Linewidth', 2);
plot(relativeRange, RISGain(:,1),'k-', 'Linewidth', 2);
plot(relativeRange, RISGain(:,2),'b--', 'Linewidth', 2);
plot(relativeRange, RISGainFarfield,'r-.', 'Linewidth', 2);
set(gca,'XScale','log');
xticks([10 100 1000 10000 100000])
xticklabels({'$10 d_F$','$10^2 d_F$','$10^3 d_F$','$10^4 d_F$','$10^5 d_F$'})
ylim([0 1]);
yticks([0 0.25 0.5 0.75 1]);
legend({'Maximum gain', 'Focus at $z=d_B$','Focus at $z=d_{FA}/10$','Focus at $z=\infty$'},'Interpreter','Latex','Location','Best');
xlabel('Propagation distance ($z$)','Interpreter','Latex');
ylabel('Normalized RIS gain','Interpreter','Latex');
set(gca,'fontsize',18);


