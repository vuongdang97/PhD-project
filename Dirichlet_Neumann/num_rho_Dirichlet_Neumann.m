function rho = num_rho_Dirichlet_Neumann(xi,theta,T,n,m,a,b)
% This function is to compute the numerical convergence factor of Dirichlet
% Neumann algorithm with toy model alpha = 1, gamma = 0, given T and no
% overlapping delta = 0
% Test the computation with sine and cosine of matrices
% m is the number of mesh point in time
% n is the number of mesh point in space
% (a,b) is space interval
% Fixed parameter
global alpha gamma delta
alpha = 1;
gamma = 0;
delta = 0;
mu = 1;
error_FVM = [];
error_control = [];
error_max = [];
itermax = 4;
test = 1;
T1 = T/2;
TOL = 10^(-50);
dx = (b-a)/(n-1);
dt = (b-a)/(m-1);
x = (a:dx:b);
t = (0:dt:T);
% Time domain 1
t1 = (0:dt:T1);
% Time domain 2
t2 = (T1:dt:T);

% Boundary data for domain 1
bl1 = zeros(1,length(t1));
br1 = zeros(1,length(t1));

% Boundary data for domain 2
bl2 = zeros(1,length(t2)); 
br2 = zeros(1,length(t2));

% Define necessary data 
yhat1 = zeros(n,length(t1));
yhat2 = zeros(n,length(t2));
z2 = zeros(length(x),1);
y01 = zeros(1,length(x)); 
y01t = zeros(1,length(x));

% Test for example y(x,t) = sin(pi*x)*(t-T)
% z2 = -alpha*pi^2*sin(10*pi*x); % domain 2

s1 = dt/dx;

% Create the matrix C
e = ones(n-2,1); 
C = spdiags([s1/2*e -s1*e s1/2*e], -1:1, n-2, n-2);


% % Initial data for domain 1      
% % Test for example y(x,t) = sin(pi*x)(t-T)
%        y01 = -sin(10*pi*x)*T; 
%       y01t = sin(10*pi*x);
% for i = 1:length(x)
% % Test for example y(x,t) = sin(pi*x)*(t-T)
%          yhat1(i,:) = (alpha*pi^4+1)*sin(10*pi*x(i))*(t1-T);
%          yhat2(i,:) = (alpha*pi^4+1)*sin(10*pi*x(i))*(t2-T);
% end
    
% Define the exact solution y and lambda on two domain
yexact1 = zeros(n,length(t1));
yexact2 = zeros(n,length(t2));
Lexact1 = zeros(n,length(t1));
Lexact2 = zeros(n,length(t2));

% yexact1 = Y(:,1:length(t1));
% yexact2 = Y(:,end-length(t2)+1:end);
% Lexact1 = Lambda(:,1:length(t1));
% Lexact2 = Lambda(:,end-length(t2)+1:end);
%for i = 1:length(x)        
% Exact solution for example sin(pi*x)*(t-T)
%           yexact1(i,:) = sin(10*pi*x(i))*(t1-T);
%           yexact2(i,:) = sin(pi*x(i))*(t2-T);
%           Lexact1(i,:) = alpha*pi^2*sin(pi*x(i))*(t1-T);
%           Lexact2(i,:) = alpha*pi^2*sin(pi*x(i))*(t2-T);
%end

% Define the discrete solution
Y1_old = zeros(n,length(t1));
Lambda1_old = zeros(n,length(t1)); 
Y1_new = zeros(n,length(t1));
Lambda1_new = zeros(n,length(t1));
Y2_old = zeros(n,length(t2));
Lambda2_old = zeros(n,length(t2)); 
Y2_new = zeros(n,length(t2));
Lambda2_new = zeros(n,length(t2));

% Data lambda_t(M+1)
lam_T2 = zeros(n,1); % the final value of lambda for domain 2
iter = 0;
IT = [];
%if(temp == 1)
g_lam = sin(k*x);
%else
%g_lam = rand(1,length(x));
%end
%lam_T1 = zeros(1, length(x));
g_y = zeros(1,length(x));
%y_dom2 = sin(k*x);
g_y = g_y';

% lam_T1 = zeros(1,length(x));
% z1 = zeros(1,length(x));
% z1 = z1';
y02 = zeros(1,length(x));
% y02t = (10*pi*x);
% Solving the time DDM system
%tic
while(test > TOL && iter < itermax)
   
    % Solve the equation on domain 1
    [Y1_new,Lambda1_new] = CoupledWaveDir_Dirichlet(y01,y01t,g_y',yhat1,g_lam,mu,bl1,br1,a,b,0,T1,n,length(t1),alpha,gamma); % y01,y01t are given, we update y_dom2 and lam_T1
%     subplot(2,1,1)
%     plot(x,y_dom2);
%     hold on
%     subplot(2,1,2)
%     plot(x,lam_T1);
%     hold on
%     pause();
     % Update data for domain 2
    y02t(2:n-1) = (eye(n-2)+dt/dx*C)*Y1_new(2:n-1,end)-Y1_new(2:n-1,end-1);
    y02t(2:n-1) = 1/dt*y02t(2:n-1);
    % Applying paper of Laurence to update the data for domain two
      lam_t(2:n-1,1) = (eye(n-2)+dt/dx*C)*Lambda1_new(2:n-1,end)-Lambda1_new(2:n-1,end-1); 
      lam_t(2:n-1,1) = 1/dt*lam_t(2:n-1,1);
     % Solve the equation on domain 2
     [Y2_new,Lambda2_new] = CoupledWaveDir_Neumann(y02,y02t,z2,lam_t,yhat2,lam_T2,mu,bl1,br1,a,b,T1,T,n,length(t2),alpha,gamma); % z2 and lam_T2 are give, we update y02 and y02t         
     % Update data for domain one
    g_lam = (1-theta)*g_lam + theta*Lambda2_new(2:n-1,1);
    
    % Applying paper of Laurence to update the data for domain one
    g_y(2:n-1) = (1-theta)*g_y(2:n-1) + theta*Y2_new(2:n-1,1);

    % Compute the error
     test = norm(Y1_new-Y1_old) + norm(Y2_new-Y2_old) + norm(Lambda1_new-Lambda1_old) + norm(Lambda2_new - Lambda2_old);
     test = (dt*dx)^(1/2)*test;
    
   
    
    % Increase the iteration
     iter = iter + 1;
     IT = [IT iter];
     error_FVM = [error_FVM (dt*dx)^(1/2)*(norm(Y1_new-yexact1)+norm(Y2_new-yexact2))];
     error_max = [error_max max(max(abs(Y1_new-yexact1)))+max(max(abs(Y2_new-yexact2)))];
     error_control = [error_control (dt*dx)^(1/2)*(norm(Lambda1_new-Lexact1)+norm(Lambda2_new-Lexact2))];
%             if(iter)
%             figure('units','normalized','outerposition',[0 0 1 1]);           
%             subplot(2,1,1);
%             [TT1,XX] = meshgrid(t1,x);
%              [TT2,XX] = meshgrid(t2,x);
%             surf(TT1,XX,Y1_new);
%             xlabel('t','FontSize',20);
%             ylabel('x','FontSize',20);
%             hold on
%           
%             %subplot(2,2,2);
%             surf(TT2,XX,Y2_new);
%             title(['$(y^{n}_{1},y^{n}_{2})$, n = ',num2str(iter)],'fontsize',20,'interpreter','latex');          
%             set(gca,'FontSize',20);
%             %zlim([-0.5 0.5]);
%               pause();
%             subplot(2,1,2);
%             surf(TT1,XX,Lambda1_new);
%             hold on
%             %subplot(2,2,4);
%             surf(TT2,XX,Lambda2_new);
%             xlabel('t','FontSize',20);
%             ylabel('x','FontSize',20);
%             title(['$(\lambda^{n}_{1},\lambda^{n}_{2})$, n = ', num2str(iter)],'fontsize',20,'interpreter','latex');          
%             set(gca,'FontSize',20);
%         %    zlim([-2 2]);
%               pause();
%             %print -depsc iterates_DN_6.eps
            %  end
 % Update new solution Y and Lambda
    Y1_old = Y1_new; 
    Lambda1_old = Lambda1_new; 
    Y2_old = Y2_new; 
    Lambda2_old = Lambda2_new;
end
%rho = mean(error_FVM(2:end)./error_FVM(1:end-1));
if(iter >3)
    rho = error_FVM(4)/error_FVM(3);
else
    rho = error_FVM(end)/error_FVM(end-1);
end
%     if(temp == 1)
%         semilogy(IT(1:end),error_FVM(1:end),'-*','DisplayName',['k = ',num2str(k)],'LineWidth',2.0);
%         hold on
%         legend show 
%         set(gca,'FontSize',20);
%         xlabel('iteration','FontSize',20);
%         ylabel('error','FontSize',20);
%     else
%         semilogy(IT(1:end),error_FVM(1:end),'-*','DisplayName',['rand'],'LineWidth',2.0);
%         hold on 
%         legend show
%         set(gca,'FontSize',20);
%         xlabel('iteration','FontSize',20);
%         ylabel('error','FontSize',20);
%     end

