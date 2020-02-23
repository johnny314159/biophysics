%PHYS 319 Assignment 3 - Johnny Liu
%Feb 22, 2020
%Q2

close all, clear all
axisLim =[0 3 0 1;0 0 0 0;0 2 0 1;0 0 0 0;0 1 0 1];

for E = 0:2:4
    x = 0:0.001:3;
   
    Z = 1+(exp(-2))+2*x*(1+exp(-(2+(-E))))+(x.^2)*(1+exp(-(2+2*(-E)))); 
    p0 = (1+exp(-2))./Z;
    p1 = ((2*x)*(1+exp(-(2+(-E))))./Z);
    p2 = ((x.^2)*(1+exp(-(2+2*(-E))))./Z);
    
    h = figure;
    plot(x,p0,'--','LineWidth',1)
    hold on
    plot(x,p1,'-.','LineWidth',1)
    plot(x,p2,'LineWidth',1)
    hold off
    xlabel('Dimensionless Concentration','Interpreter','latex','FontSize',12)
    ylabel('Probability','Interpreter','latex','FontSize',12)
    title(sprintf('Probability of Dimoglobin Binding States for $$\\varepsilon = %d k_{b}T$$',E),'Interpreter','latex','FontSize',12);
    axis(axisLim((E+1),:))
    set(gca,'TickLabelInterpreter','latex')
    legend('$p_{0}$','$p_{1}$','$p_{2}$','Interpreter','latex','Location','best')
    saveas(h,sprintf('PHYS319_HW3_Q2_FIG%d.png',E))
    clear p0 p1 p2
end

%%
%4A
clear all
close all
numSteps = 1000;
numFrames = numSteps+1;
coord=zeros(3,numFrames);                                               %3xnumSteps+1 because the 1st coordinate is (0,0,0)
dice=zeros(1,numSteps);

for n = 1:numSteps
    dice(n) = randi(6);
    if dice(n) ==1
        coord(1,n+1:end) = coord(1,n) +1;
    elseif dice(n) ==2
        coord(1,n+1:end) = coord(1,n) -1;
    elseif dice(n) ==3
        coord(2,n+1:end) = coord(2,n) +1;
    elseif dice(n) ==4
        coord(2,(n+1):end) = coord(2,n) -1;
    elseif dice(n) ==5
        coord(3,n+1:end) = coord(3,n) +1;
    else coord(3,n+1:end) = coord(3,n) -1;
    end
   
end
H = figure;
plot3(coord(1,:),coord(2,:),coord(3,:))
xlabel('x (pixel)','Interpreter','latex','FontSize',12)
ylabel('y (pixel)','Interpreter','latex','FontSize',12)
zlabel('z (pixel)','Interpreter','latex','FontSize',12)
title(sprintf('3D Random Walk Trajectory for N = %d Steps (Step Size =1)',numSteps),'Interpreter','latex','FontSize',12)
saveas(H,'4A.png')

%4B
numerator = 0;

for tau = 1:(numFrames-1)
    for i = (tau+1):numFrames
        numerator = numerator + (coord(1,i)-coord(1,i-tau))^2 + (coord(2,i)-coord(2,i-tau))^2 + (coord(3,i)-coord(3,i-tau))^2;
    end
    MSD(1,tau) = (numerator)/(numFrames-tau);
    numerator = 0;
end
tauX = 1:numFrames-1;
polyCo = polyfit(tauX,MSD,1);
J = figure;
plot(tauX,MSD)
hold on
plot(tauX,(polyval(polyCo,tauX)))
ylabel('Mean squared displacement (pixels$^{2}$)','Interpreter','latex','FontSize',12)
xlabel('Time lag $\tau$ (time unit)','Interpreter','latex','FontSize',12)
title('Mean Squared Displacement at Varying Time Lags (Step Size = 1)','Interpreter','latex','FontSize',12)
legend('MSD','1st-order fit','Interpreter','latex','location','best','FontSize',12)
hold off
gtext(sprintf('Diffusion Coefficient = %d pixels$^{2}$/units time',(polyCo(1,1)/6)),'Interpreter','latex')
saveas(J,'4B.png')

%4C
for j =1:50:numSteps
    coord_50(1,1+((j-1)/50)) = coord(1,j);
    coord_50(2,1+((j-1)/50)) = coord(2,j);
    coord_50(3,1+((j-1)/50)) = coord(3,j);
end

L = figure;
plot3(coord_50(1,:),coord_50(2,:),coord_50(3,:))
xlabel('x (pixel)','Interpreter','latex','FontSize',12)
ylabel('y (pixel)','Interpreter','latex','FontSize',12)
zlabel('z (pixel)','Interpreter','latex','FontSize',12)
title(sprintf('3D Random Walk Trajectory for N = %d Steps (Step Size = 50)',numSteps),'Interpreter','latex','FontSize',12)
saveas(L,'4C_walk.png')

numFrames_50 = length(coord_50);
num_50 = 0;
for tau = 1:(numFrames_50-1)
    for i = (tau+1):numFrames_50
        num_50 = num_50 + (coord(1,i)-coord(1,i-tau))^2 + (coord(2,i)-coord(2,i-tau))^2 + (coord(3,i)-coord(3,i-tau))^2;
    end
    mSD_50(1,tau) = (num_50)/(numFrames_50-tau);
    num_50 = 0;
end
X_50 = 1:numFrames_50-1;
polyCo_50 = polyfit(X_50,mSD_50,1);
M = figure;
plot(X_50,mSD_50,'o')
hold on
plot(X_50,polyval(polyCo_50,X_50))
ylabel('Mean squared displacement (pixels$^{2}$)','Interpreter','latex','FontSize',12)
xlabel('Time lag $\tau$ (time unit)','Interpreter','latex','FontSize',12)
title('Mean Squared Displacement at Varying Time Lags (Step Size = 50)','Interpreter','latex','FontSize',12)
legend('MSD','1st-order fit','Interpreter','latex','location','best','FontSize',12)
gtext(sprintf('Diffusion Coefficient = %d pixels$^{2}$/units time',(polyCo_50(1,1)/6)),'Interpreter','latex')
hold off
saveas(M,'4C_MSD.png')

%4D
aVG = zeros(1,200);
for N = 101:200 %# of steps
    dice = zeros(1,N);
    dist = zeros(1,100);
    for k = 1:100 %# of walks
        XYZ = zeros(3,N+1);
        for n=1:N %current step
            dice(n) = randi(6);
            if dice(n) ==1
                XYZ(1,n+1:end) = XYZ(1,n) +1;
            elseif dice(n) ==2
                XYZ(1,n+1:end) = XYZ(1,n) -1;
            elseif dice(n) ==3
                XYZ(2,n+1:end) = XYZ(2,n) +1;
            elseif dice(n) ==4
                XYZ(2,(n+1):end) = XYZ(2,n)-1;
            elseif dice(n) ==5
                XYZ(3,n+1:end) = XYZ(3,n) +1;
            else XYZ(3,n+1:end) = XYZ(3,n) -1;
            end
        end
        dist(k)= ((XYZ(1,end))^2 + (XYZ(2,end))^2 + (XYZ(3,end))^2)^(0.5);  
    end
    AVG(1,N)=sum(dist)/100;
end
AVG_use = aVG(1,[101:200]);
NAVG = 101:200;
N = figure
plot(NAVG,AVG_use)
xlabel('Number of steps per walk (N)','Interpreter','latex','FontSize',12)
ylabel('Mean displacement (pixels)','Interpreter','latex','FontSize',12)
title('Average End-to-end Distance for 3D Random Walks between 101 and 200 Steps (n=100)','Interpreter','latex','FontSize',12)
saveas(N,'4D.png')