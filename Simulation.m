%%% ME354 - ASSIGNMENT 
%%% Vibrations
%%% Abhinav Maheshwari

clc;
clear all;
close all;

m1 = 1;
m2 = 1;
k = 1;
l0 = 1;
alpha = 0.1;
beta = 1/500;

%%% INITIAL CONDITIONS
ch = 'case5';
switch ch
    case 'case1'
        % CASE 1
        x1_0 = alpha;
        x2_0 = alpha + l0 + .25*l0;
        y1_0 = 0;
        y2_0 = 0;
        
        u1_0 = 1/500;
        v1_0 = 1/500;
        u2_0 = 1/500;
        v2_0 = 1/500;
    case 'case2'
        %CASE 2
        x1_0 = alpha;
        x2_0 = alpha + l0;
        y1_0 = 0;
        y2_0 = 0;
        
        u1_0 = beta;
        v1_0 = 1/500;
        u2_0 = beta;
        v2_0 = 1/500;
    case 'case3'
        %CASE 3
        x1_0 = alpha;
        x2_0 = alpha + l0 + .25*l0;
        y1_0 = 0;
        y2_0 = 0;
        
        u1_0 = beta;
        v1_0 = 1/500;
        u2_0 = beta;
        v2_0 = 1/500;
    case 'case4'
        %CASE 4;
        x1_0 = alpha;
        x2_0 = alpha + l0 + .25*l0;
        y1_0 = 0;
        y2_0 = 0;
        
        u1_0 = 0;
        v1_0 = -beta;
        u2_0 = 0;
        v2_0 = beta;
    case 'case5'
        % CASE 5
        x1_0 = alpha;
        x2_0 = alpha + l0;
        y1_0 = 0;
        y2_0 = 0;
        
        u1_0 = 0;
        v1_0 = -beta;
        u2_0 = 0;
        v2_0 = beta;
end
tspan = [0,2000];
X0 = [x1_0 y1_0 x2_0 y2_0 u1_0 v1_0 u2_0 v2_0]; %%INITIAL CONDITIONS
[T, X_sol] = ode45(@(t,x)odefun(m1,m2,k,l0, x, t), tspan, X0);


%%% SIMULATION

x1 = X_sol(:,1);
y1 = X_sol(:,2);
x2 = X_sol(:,3);
y2 =  X_sol(:, 4);
%%% we know that Mass = density*volume
%%% therefore, radius = constant*mass^(1/3)
constant = 0.1;
Rad_1 = constant*m1^(1/3);
Rad_2 = constant*m2^(1/3);
theta = linspace(0, 2*pi, 100);
figure(1);
grid on;
grid minor;
daspect([1 1 1]);
xlabel('$\textit{\textbf{x-axis}}$','Interpreter', 'latex', 'FontSize',11)
ylabel('$\textit{\textbf{y-axis}}$','Interpreter', 'latex', 'FontSize',11)
xlim([-5,5]);
ylim([-5,5]);
hold on;
title('$\textbf{ME354 - GROUP \#07}$','Interpreter', 'latex','FontSize',18);
mass_1 = fill(Rad_1*cos(theta)+x1(1), Rad_1*sin(theta)+y1(1), 'R', DisplayName='Mass m1');
mass_2 = fill(Rad_2*cos(theta)+x2(1), Rad_2*sin(theta)+y2(1), 'B', DisplayName='Mass m2');
legend;
% v = VideoWriter('Group7.avi');
% open(v);
pause(0.1);
for i = 1:length(T)
%     hold on
    hg = hggroup;
    mass_1.XData = Rad_1*cos(theta) +x1(i);
    mass_1.YData = Rad_1*sin(theta) + y1(i);
    mass_2.XData = Rad_2*cos(theta) + x2(i);
    mass_2.YData = Rad_2*sin(theta) + y2(i);
    xa = x1(i); ya = y1(i); xb = x2(i); yb = y2(i); ne = 10; a = 1; ro = 0.1;
    [xs_,ys_] = spring(xa,ya,xb,yb,ne,a,ro); plot(xs_,ys_,'LineWidth',0.2, 'Parent', hg, 'Color', 'black');
%     hold off
%     pause(5);
    drawnow limitrate;
%     frame = getframe(gcf);
%     writeVideo(v, frame);
    if i<length(T)
        delete(hg);
    end
end
% close(v);

function [xs ys] = spring(xa,ya,xb,yb,varargin)
persistent ne Li_2 ei b

if nargin > 4 % calculating some fixed spring parameters only once time
    [ne a r0] = varargin{1:3};                  % ne: number of coils - a = natural length - r0 = natural radius
    Li_2 =  (a/(4*ne))^2 + r0^2;                % (large of a quarter of coil)^2
    ei = 0:(2*ne+1);                            % vector of longitudinal positions
    j = 0:2*ne-1; b = [0 (-ones(1,2*ne)).^j 0]; % vector of transversal positions
end
R = [xb yb] - [xa ya]; mod_R = norm(R); % relative position between "end_B" and "end_A"
L_2 = (mod_R/(4*ne))^2; % (actual longitudinal extensión of a coil )^2
if L_2 > Li_2
   error('Spring:TooEnlargement', ...
   'Initial conditions cause pulling the spring beyond its maximum large. \n Try reducing these conditions.')
else
    r = sqrt(Li_2 - L_2);   %actual radius
end
c = r*b;    % vector of transversal positions
u1 = R/mod_R; u2 = [-u1(2) u1(1)]; % unitary longitudinal and transversal vectors 
xs = xa + u1(1)*(mod_R/(2*ne+1)).*ei + u2(1)*c; % horizontal coordinates
ys = ya + u1(2)*(mod_R/(2*ne+1)).*ei + u2(2)*c; % vertical coordinates
end

function [x_dot] = odefun(m1, m2, k, l0, X, t)
    x_dot = zeros(8,1);
    x_dot(1,1) = X(5);
    x_dot(2,1) = X(6);
    x_dot(3,1) = X(7);
    x_dot(4,1) = X(8);
    x_dot(5,1) = (1/m1)*k*(sqrt((X(1) - X(3))^2 + (X(2)-X(4))^2) - l0)*(X(3) - X(1))/(sqrt((X(1) - X(3))^2 + (X(2)-X(4))^2));
    x_dot(6,1) = x_dot(5)*(X(4)-X(2))/(X(3) - X(1));
    x_dot(7,1) = -(m1/m2)*x_dot(5,1);
    x_dot(8,1) = -(m1/m2)*x_dot(6,1);
end

%%% THANK YOU!