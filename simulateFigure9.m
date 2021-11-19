%%This Matlab script generates Figure 9 in the paper:
%
%Emil Björnson, Özlem Tuğfe Demir, and Luca Sanguinetti, 
%"A Primer on Near-Field Beamforming for Arrays and Reconfigurable 
%Intelligent Surfaces,"  Asilomar Conference on Signals, Systems, and
%Computers, Virtual conference, October-November 2021.
%
%Download article: https://arxiv.org/pdf/2110.06661.pdf
%
%This is version 1.1 (Last edited: 2021-11-19)
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

%Define the range of points along the z-axis
relativeRange = sort([logspace(1,5,50) distance_B N/10]);



%Define the range of points along the horizontal axis
relativeRangeX = sort([  -logspace(log10(4000),-2,26) 0 logspace(-2,log10(4000),26)  ]);

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
focusIndex = [find(relativeRange==distance_B) find(relativeRange==N/10) ]; 

XfocusIndex = find(xRange==0);

%Determine the RIS gains with different phase shifts
RISGain = zeros(length(relativeRangeX),length(relativeRange),length(focusIndex)+1);

for n = 1:length(focusIndex)
    
    phaseshifts = exp(1j*angle(channelVectors(:,focusIndex(n),XfocusIndex)));
    
    for mx = 1:length(xRange)
        RISGain(mx,:,n) = abs(channelVectors(:,:,mx)'*(phaseshifts)).^2;
    end
end

%Determine the phase-shifts for focusing at far-field
phaseshifts0 = exp(1j*angle(channelVectors0));

for mx = 1:length(xRange)
    RISGain(mx,:,end) = abs(channelVectors(:,:,mx)'*(phaseshifts0)).^2;
end

set(groot,'defaultAxesTickLabelInterpreter','latex');


%% Plot the simulation results
[X,Y] = meshgrid(relativeRangeX,relativeRange);

% Heatmaps
%(a)
figure;
hold on; box on; grid on;
contourf(X,Y,RISGain(:,:,1).',100,'LineColor','none')
set(gca,'YScale','log');
xlim([-40/fraunhoferDistanceElement 40/fraunhoferDistanceElement])
xticks([  -2000 -1000 0  1000 2000 ])
xticklabels({'$-2\cdot10^3 d_F$','$-10^3 d_F$','$0 d_F$','$10^3 d_F$','$2\cdot10^3 d_F$'})
yticks([10 100 1000 10000 100000])
yticklabels({'$10 d_F$','$10^2 d_F$','$10^3 d_F$','$10^4 d_F$','$10^5 d_F$'})
xlabel('$x_r$','Interpreter','Latex');
ylabel('$z_r$','Interpreter','Latex');
plot(linspace(-D_element*sqrt(N/8)/fraunhoferDistanceElement, D_element*sqrt(N/8)/fraunhoferDistanceElement,2), [1 1],'r', 'Linewidth', 4)
set(gca,'fontsize',24);
colorbar
annotation('arrow', [.51 .55], [.45 .46],'Color','white');
axes('pos',[.56 .17 .3 .3])
contourf(X,Y,RISGain(:,:,1).',100,'LineColor','none')
set(gca,'YScale','log');
set(gca,'fontsize',18);
xlim([-1.77*distance_B*sqrt(lambda/fraunhoferDistanceElement/N) 1.77*distance_B*sqrt(lambda/fraunhoferDistanceElement/N)])
ylim([200 1000])
xticks([-1.77/2*distance_B*sqrt(lambda/fraunhoferDistanceElement/N) 0  1.77/2*distance_B*sqrt(lambda/fraunhoferDistanceElement/N) ])
xticklabels({'$-10 d_F$', '$0 d_F$','$10 d_F$'})
yticks([N*distance_B/(N+10*distance_B) distance_B N*distance_B/(N-10*distance_B) ])
yticklabels({'$286 d_F$','$400 d_F$', '$667 d_F$'})
rectangle('Position',[-1.77/2*distance_B*sqrt(lambda/fraunhoferDistanceElement/N) N*distance_B/(N+10*distance_B) 1.77*distance_B*sqrt(lambda/fraunhoferDistanceElement/N) N*distance_B/(N-10*distance_B)-N*distance_B/(N+10*distance_B)],'Linewidth', 2)
set(gca,'XColor',[1 1 1]); 
set(gca,'YColor',[1 1 1]); 


%(b)
figure;
hold on; box on; grid on;
contourf(X,Y,RISGain(:,:,2).',100,'LineColor','none')
set(gca,'YScale','log');
xlim([-40/fraunhoferDistanceElement 40/fraunhoferDistanceElement])
xticks([  -2000 -1000 0  1000 2000 ])
xticklabels({'$-2\cdot10^3 d_F$','$-10^3 d_F$','$0 d_F$','$10^3 d_F$','$2\cdot10^3 d_F$'})
yticks([10 100 1000 10000 100000])
yticklabels({'$10 d_F$','$10^2 d_F$','$10^3 d_F$','$10^4 d_F$','$10^5 d_F$'})
xlabel('$x_r$','Interpreter','Latex');
ylabel('$z_r$','Interpreter','Latex');
plot(linspace(-D_element*sqrt(N/8)/fraunhoferDistanceElement, D_element*sqrt(N/8)/fraunhoferDistanceElement,2), [10 10],'r', 'Linewidth', 4)
set(gca,'fontsize',24);
colorbar
annotation('arrow', [.51 .55], [.51 .47],'Color','white');
axes('pos',[.56 .17 .3 .3])
contourf(X,Y,RISGain(:,:,2).',100,'LineColor','none')
set(gca,'YScale','log');
set(gca,'fontsize',18);
xlim([-1.77*(N/10)*sqrt(lambda/fraunhoferDistanceElement/N) 1.77*(N/10)*sqrt(lambda/fraunhoferDistanceElement/N)])
ylim([400 20000])
xticks([-1.77/2*(N/10)*sqrt(lambda/fraunhoferDistanceElement/N) 0  1.77/2*(N/10)*sqrt(lambda/fraunhoferDistanceElement/N) ])
xticklabels({'$-25 d_F$', '$0 d_F$','$25 d_F$'})
yticks([N*(N/10)/(N+10*(N/10)) (N/10) N ])
yticklabels({'$50 d_F$','$10^3 d_F = d_{FA}/10$', '$10^4 d_F = d_{FA}$'})
rectangle('Position',[-1.77/2*(N/10)*sqrt(lambda/fraunhoferDistanceElement/N) N*(N/10)/(N+10*(N/10)) 1.77*(N/10)*sqrt(lambda/fraunhoferDistanceElement/N) 3*N-N*(N/10)/(N+10*(N/10))],'Linewidth', 2)
set(gca,'XColor',[1 1 1]); 
set(gca,'YColor',[1 1 1]); 


%(c)
figure;
hold on; box on; grid on;
contourf(X,Y,RISGain(:,:,3).',100,'LineColor','none')
set(gca,'YScale','log');
xlim([-40/fraunhoferDistanceElement 40/fraunhoferDistanceElement])
xticks([  -2000 -1000 0  1000 2000 ])
xticklabels({'$-2\cdot10^3 d_F$','$-10^3 d_F$','$0 d_F$','$10^3 d_F$','$2\cdot10^3 d_F$'})
yticks([10 100 1000 10000 100000])
yticklabels({'$10 d_F$','$10^2 d_F$','$10^3 d_F$','$10^4 d_F$','$10^5 d_F$'})
xlabel('$x_r$','Interpreter','Latex');
ylabel('$z_r$','Interpreter','Latex');
plot(linspace(-D_element*sqrt(N/8)/fraunhoferDistanceElement, D_element*sqrt(N/8)/fraunhoferDistanceElement,2), [10 10],'r', 'Linewidth', 4)
set(gca,'fontsize',24);
colorbar

