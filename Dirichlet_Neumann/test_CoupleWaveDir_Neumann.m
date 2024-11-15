% This script is created to test the coupled wave equation with the
% function CoupledWaveDir_Neumann. 
% Particularly, we deal with the exact solution  
% y(x,t) = sin(pi*x)*(t-T)
% lambda(x,t) = alpha*pi^2*sin(pi*x)*(t-T)
% y hat(x,t) = (alpha*pi^4+1)*sin(pi*x)*(t-T)
% z(x) = -alpha*pi^2*sin(pi*x)

% We simulate the solution with each numerical test for some different
% values of T 
% Moreover, we compute the l2 norm error and determine the order of
% convergence
% We use the coupledwavedir1 function with the block structure proposed by
% Laurence Halpern
clc;  
close all;
clear all;

% Fixed parameter
a = 0; 
b = 1;
T0 = 0;
alpha = 1; 
gamma = 0;
mu = 1;
error = [];
error_control = [];
error_max = [];
relative = [];
%N = [51 101 201 401];
%M = [51 101 201 401];
N = [21 41 81 161];
M = [21 41 81 161];

for k = 1:length(N)
    T = 1;
    n = N(k);
    m = M(k);
    dx = (b-a)/(n-1);
    dt = T/(m-1);
    bl = zeros(1,m);
    br = zeros(1,m);
    x = (a:dx:b);
    t = (0:dt:T);

% Define the y hat solution 
yhat = zeros(n,m);
z = zeros(n,1); %gamma*y(x,T) - lambda_t(x,T)

      
      % Test for example y(x,t) = sin(pi*x)*(t-T)
        y0 = -sin(pi*x)*T; %y(x,0)
        y0t = sin(pi*x); %y_t(x,0)
        %y0t = 100*rand(length(x));
             
for i = 1:length(x)
        % Test for example y(x,t) = sin(pi*x)*(t-T)
        yhat(i,:) = (alpha*pi^4+1)*sin(pi*x(i))*(t-T);
end
     % Test for example y(x,t) = sin(pi*x)*(t-T)
    z = -alpha*pi^2*sin(pi*x);
    z = z';   
    % Define the exact solution y 
yexact = zeros(n,m);
Lexact = zeros(n,m);
for i = 1:length(x)
        % Exact solution for example sin(pi*x)*(t-T)
            yexact(i,:) = sin(pi*x(i))*(t-T);
            Lexact(i,:) = alpha*pi^2*sin(pi*x(i))*(t-T);
end
lam_T = Lexact(:,end);
lam_dom1 = alpha*pi^2*sin(pi*x);
temp = 1;
[Y,Lambda] = CoupledWaveDir_Neumann(y0,y0t,z,lam_dom1',yhat,lam_T,mu,bl,br,a,b,T0,T,n,m,alpha,gamma);
% for j=1:length(t)
%     subplot(2,1,1);
%     % Simulate the y
%     plot(x,Y(:,j),x,yexact(:,j));
%     legend('Y','Yexact');
%     %axis([0 1 -1 1]);
%     % Simulate the control
%     subplot(2,1,2);
%     plot(x,Lambda(:,j),x,Lexact(:,j));
%     legend('Lambda','Lambdaexact');
%    % axis([0 1 -5 5]);
%     title(['time = ', num2str(t(j))]);
%     pause(0.5) 
% end
% [XX,TT] = meshgrid(x,t);
% subplot(2,2,1);
% surf(Y);
% title('Y');
% xlabel('t');
% ylabel('x');
% subplot(2,2,2);
% surf(yexact);
% title('yexact');
% xlabel('t');
% ylabel('x');
% subplot(2,2,3);
% surf(Lambda);
% title('Lambda');
% xlabel('t');
% ylabel('x');
% subplot(2,2,4);
% surf(Lexact);
% title('exact Lambda');
% xlabel('t');
% ylabel('x');
   error = [error (dt*dx)^(1/2)*norm(Y-yexact)];
   error_max = [error_max max(max(abs(Y-yexact)))];
   error_control = [error_control (dt*dx)^(1/2)*norm(Lambda-Lexact)];
end
% Plot the order of convergence
figure
loglog(N,N.^(-2),N,error,'-x',N,error_max,'-+',N,error_control,'-o','LineWidth',2.0);
legend({'2x','error','infinity error','control error'},'FontSize',20);
title('The convergence plot with dt = 0.8dx','FontSize',20);
set(gca,'FontSize',20);

