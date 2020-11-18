%Brian Lui
%bl569
%MAE 5730
%HW due 9/12/18
%Mass on Spring 3D- #9

clear all;
close all;

%toggleable items
p.k = 2;         %spring constant
p.m = 1;
p.g = 9.8;

dur = 10;
npoints = 1001;
tspan = linspace(0, dur, npoints);

%initial condition stored as a matrix, the rows each represent a state, 
%the columns represent x, xDot, y, yDot, z, zDot at time = 0
ic = [0,0,0,0,0,1;      %point in xy, oval in xz, yz
   1,0,0,0,0,0;...      %line in xy (1 oscillates but the other doesn't)
   1,0,1,0,0,0;...      %line in xy (in phase). Also a line in xz and yz (in phase)
   1,0,0,1,0,0;...      %circle in xy, oval in yz, line in xz (check scale!)
   1,0,0,0,0,1;...      %oval in xz, line in xy, line in yz
   3,5,1,2,4,6];        %oval in all 3 directions


for i = 1:length(ic(:,1))
    state0 = ic(i,:);
    f = @(t,state) rhs(state, p);
    [tArray, stateArray] = ode45(f,tspan, state0);
    
    xPlot = stateArray(:,1);
    yPlot = stateArray(:,3);
    zPlot = stateArray(:,5);
    
%     %visualizing how x,y, and z changes over time
%     figure(2)
%     subplot(3,1,1);
%     plot(tArray, xPlot);
%     title('x vs t');
%     subplot(3,1,2);
%     plot(tArray, yPlot);
%     title('y vs t');
%     subplot(3,1,3);
%     plot(tArray, zPlot);
%     title('z vs t');
    
    figure(i);
    subplot(3,1,1);
    plot(xPlot, yPlot);
    title('xy plane');
    xlabel('x')
    ylabel('y')
    
    subplot(3,1,2);
    plot(xPlot, zPlot);
    title('xz')
    xlabel('x')
    ylabel('z');
    
    subplot(3,1,3);
    plot(yPlot, zPlot);
    title('yz');    
    xlabel('y')
    ylabel('z')
    
end

    
function stateDot = rhs(state,p)
    x = state(1); xDot = state(2);
    y = state(3); yDot = state(4);
    z = state(5); zDot = state(6);
    k = p.k;
    m = p.m;
    g = p.g;
    
    xDoubleDot = -k/m * x;
    yDoubleDot = -k/m*y;
    zDoubleDot = -k/m*z - g;
    stateDot = [xDot; xDoubleDot; yDot; yDoubleDot; zDot; zDoubleDot];
    
end

