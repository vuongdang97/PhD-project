% This script is created to test the coupled wave equation with the
% function CoupledWaveDir. 
% Particularly, we deal with the manufactured solution y(x,t) = sin(pi*x)*(t-T)
% And the original and dual solution denoted by y and lambda.
% The dual solution of 
% lambda(x,t) = alpha*pi^2*sin(pi*x)*(t-T)
% and target state
% y hat(x,t) = (alpha*pi^4+1)*sin(pi*x)*(t-T)
% z(x) = -alpha*pi^2*sin(pi*x)  (gamma = 0)

% We simulate the solution with each numerical test for some different
% values of T 
% Moreover, we compute the l2 norm error of y and determine the order of
% convergence
clc;  
close all;
clear all;

% Fixed parameter
a = 0; 
b = 1;
T0 = 0;
TT = 1:1:50;
alpha = 1; 
gamma = 0;
mu = 1;
error = [];
error_control = [];
error_max = [];
relative = [];
N = [51 101 201 401];
M = [51 101 201 401];


for l = 1:length(N)
    T = 1;
    n = N(l);
    m = M(l);
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
        y0 = -sin(pi*x)*T; 
        y0t = sin(pi*x);
            
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
[Y,Lambda] = CoupledWaveDir_FVM1(y0,y0t,z,yhat,lam_T,mu,bl,br,a,b,T0,T,n,m,alpha,gamma);
%    myVideo = VideoWriter('myfile.avi');
%    uncompressedVideo = VideoWriter('myfile.avi', 'Uncompressed AVI');
%    myVideo.FrameRate = 10;
%    open(myVideo);
for j=1:length(t)
    subplot(2,1,1);
    % Simulate the y
    plot(x,Y(:,j));
    title('Y','FontSize',30);
    axis([0 1 -1 1]);
    % Simulate the control
    subplot(2,1,2);
    plot(x,Lambda(:,j));
    title('\Lambda','FontSize',30);
   % axis([0 1 -5 5]);
    %title(['time = ', num2str(t(j))]);
    pause(0.5) 
    frame=getframe(gcf);
    writeVideo(myVideo,frame);
end
%    close(myVideo);

   error = [error (dt*dx)^(1/2)*norm(Y-yexact)];
   error_max = [error_max max(max(abs(Y-yexact)))];
   error_control = [error_control (dt*dx)^(1/2)*norm(Lambda-Lexact)];
end



% Plot the order of convergence
figure
loglog(N,N.^(-2),N,error,'-x',N,error_max,'-+',N,error_control,'-o','LineWidth',2.0);
legend({'2x','||y_d - y_{exact}||_{l^{2}}','||y_d - y_{exact}||_{l^{\infty}}','||\lambda_d - \lambda_{exact}||_{l^{2}}'},'FontSize',30);
xlabel('N','FontSize',30);
ylabel('error','FontSize',30);
%title('The convergence plot with dt = 0.8dx','FontSize',20);
set(gca,'FontSize',30);
% [XX,TT] = meshgrid(x,t);
% figure;
% subplot(2,2,1);
% surf(XX,TT,Y);
% xlabel('t','FontSize',20);
% ylabel('x','FontSize',20);
% zlabel('value','FontSize',20);
% title('Discrete solution Y','Fontsize',20);
% subplot(2,2,2);
% surf(XX,TT,yexact);
% xlabel('t','FontSize',20);
% ylabel('x','FontSize',20);
% zlabel('value','FontSize',20);
% title('Exact solution Y','FontSize',20);
% subplot(2,2,3);
% surf(XX,TT,Lambda);
% xlabel('t','FontSize',20);
% ylabel('x','FontSize',20);
% zlabel('value','FontSize',20);
% title('Discrete control','FontSize',20);
% subplot(2,2,4);
% surf(XX,TT,Lexact);
% xlabel('t','FontSize',20);
% ylabel('x','FontSize',20);
% zlabel('value','FontSize',20);
% title('Exact control','FontSize',20);
