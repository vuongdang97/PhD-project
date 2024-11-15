function rho = num_rho_Schwarz(xi,theta,T,n,m,a,b)
% This function is write to compute the numerical convergence factor of
% Schwarz method with toy model alpha = 1, gamma = 0, given T and no
% overlapping delta = 0
% Test the computation with sine and cosine of matrices
% m is the number of mesh point in time
% n is the number of mesh point in space
% (a,b) is space interval
global alpha gamma delta
alpha = 1;
gamma = 0;
delta = 0;

error_FVM = []; % l2 error in y
itermax = 4; % Maximum iteration
test = 1; % termination of algorithm
TOL = 10^(-100); % tolerance
dx = (b-a)/(n-1); % mesh size in space
dt = (b-a)/(m-1); % mesh size in time
T1 = T/2;
mu = 1;
x = (a:dx:b); % mesh point in space
% Time domain 1
t1 = (0:dt:T1+delta); % mesh point in domain 1
% Time domain 2
t2 = (T1:dt:T); % mesh point in domain 2
t = [t1 t2(2:end)];
% Boundary data for domain 1
bl1 = zeros(1,length(t1)); % boundary data on domain 1
br1 = zeros(1,length(t1));

% Boundary data for domain 2
bl2 = zeros(1,length(t2)); % boundary data on domain 2
br2 = zeros(1,length(t2));
bl = zeros(1,length(t));
br = zeros(1,length(t));
% Define necessary data
yhat1 = zeros(n,length(t1)); % target y on domain 1
yhat2 = zeros(n,length(t2)); % target y on domain 2
z2 = zeros(1,length(x)); % final data on domain 2
y01 = zeros(1,length(x)); % initial data on domain 1
y01t = zeros(1,length(x));
s1 = dt/dx;
%z2 = -alpha*pi^2*sin(pi*x); % domain 2
% Create the matrix C
e = ones(n-2,1);
C = spdiags([s1/2*e -s1*e s1/2*e], -1:1, n-2, n-2);
y01 = zeros(length(x),1);
y01t = zeros(length(x),1);
yhat1 = zeros(length(x),length(t1));
yhat2 = zeros(length(x),length(t2));
yhat = zeros(length(x),length(t));
%             y01 = -sin(pi*x)*T;
%             y01t = sin(pi*x);
%             for i = 1:length(x)
%                 % % Test for example y(x,t) = sin(pi*x)*(t-T)
%                 yhat1(i,:) = (alpha*pi^4+1)*sin(pi*x(i))*(t1-T);
%                 yhat2(i,:) = (alpha*pi^4+1)*sin(pi*x(i))*(t2-T);
%                 yhat(i,:) = (alpha*pi^4+1)*sin(pi*x(i))*(t-T);
% 
%             end
 % Define the exact solution y and lambda on two domain
    yexact1 = zeros(n,length(t1));
    yexact2 = zeros(n,length(t2));
    Lexact1 = zeros(n,length(t1));
    Lexact2 = zeros(n,length(t2));
    lam_T2 = zeros(n,1); % the final value of lambda for domain 2
    temp = 1;
     % Define the discrete solution
    Y1_old = zeros(n,length(t1));
    Lambda1_old = zeros(n,length(t1));
    Y1_new = zeros(n,length(t1));
    Lambda1_new = zeros(n,length(t1));
    Y2_old = zeros(n,length(t2));
    Lambda2_old = zeros(n,length(t2));
    Y2_new = zeros(n,length(t2));
    Lambda2_new = zeros(n,length(t2));
    iter = 0;
    IT = [];
    z1 = zeros(1,length(x));
     % Compute the exact solution in a discrete way
    [Y,Lambda] = CoupledWaveDir_FVM1(y01,y01t,z2',yhat,lam_T2,mu,bl,br,a,b,0,T,n,length(t),alpha,gamma);
    Y1 = Y(:,1:length(t1));
    Lambda1 = Lambda(:,1:length(t1));
    Y2 = Y(:,length(t1):end);
    Lambda2 = Lambda(:,length(t1):end);        
      % Initialize the data on domain 1

     lam_T1 = sin(xi*x);
     lam_T1(1) = 0;
     lam_T1(end) = 0;
    z1 = zeros(1,length(x));
    z1 = z1';
    error_FVM_theo = [];
    % Solving the time DDM system
     % We solve the alternate time DDM
        test = 1;
        while(test > TOL && iter < itermax) % stopping criterion
            % Solve the equation on domain 1
            [Y1_new,Lambda1_new] = CoupledWaveDir_FVM1(y01,y01t,z1,yhat1,lam_T1,mu,bl1,br1,a,b,0,T1+delta,n,length(t1),alpha,gamma); % y01,y01t are given, we update z1 and lam_T1
            
            % Update data for domain 2
            y02 = Y1_new(:,end);
            
            % Applying paper of Gander, Laurence, Nataf to update the data for domain two
            y02t(2:n-1) = (Y1_new(2:n-1,end)-Y1_new(2:n-1,end-1))/dt + 1/dx*C*Y1_new(2:n-1,end) + dt/(2*alpha)*Lambda1_new(2:n-1,end);
            y02t(1) = 0;
            y02t(n) = 0;
            
            % Solve the equation on domain 2
            [Y2_new,Lambda2_new] = CoupledWaveDir_FVM1(y02,y02t,z2',yhat2,lam_T2,mu,bl2,br2,a,b,T1,T,n,length(t2),alpha,gamma); % z2 and lam_T2 are give, we update y02 and y02t
            % Compute the error
            test = norm(Y1_new-Y1_old) + norm(Y2_new-Y2_old) + norm(Lambda1_new-Lambda1_old) + norm(Lambda2_new - Lambda2_old);
            test = (dt*dx)^(1/2)*test;
            
            % Update data for domain one
            lam_T1 = (1-theta)*lam_T1+theta*(Lambda2_new(:,1))';
            
            % Applying paper of Gander, Laurence, Nataf to update the data for domain one
            z1(2:n-1) = (1-theta)*((Lambda1_old(2:n-1,end)-Lambda1_old(2:n-1,end-1))/dt +1/dx*C*Lambda1_old(2:n-1,end) - dt/2*(Y1_old(2:n-1,end)-yhat1(2:n-1,end)))...
                +theta*((Lambda2_new(2:n-1,2)-Lambda2_new(2:n-1,1))/dt - 1/dx*C*Lambda2_new(2:n-1,1) + dt/2*(Y2_new(2:n-1,1)-yhat2(2:n-1,1)));
            z1 = -z1;
            % Update new solution Y and Lambda
            Y1_old = Y1_new;
            Lambda1_old = Lambda1_new;
            Y2_old = Y2_new;
            Lambda2_old = Lambda2_new;
            % Increase the iteration
            iter = iter + 1;
            IT = [IT iter];
            error_FVM = [error_FVM (dt*dx)^(1/2)*(norm(Y1_new-Y1)+norm(Y2_new-Y2))];
             
%             if(iter == 20)
%             figure('units','normalized','outerposition',[0 0 1 1]);           
%             subplot(2,1,1);
%             [TT1,XX] = meshgrid(t1,x);
%             [TT2,XX] = meshgrid(t2,x);
%             surf(TT1,XX,Y1_new);
%             xlabel('t','FontSize',20);
%             ylabel('x','FontSize',20);
%             hold on
%             %subplot(2,2,2);
%             surf(TT2,XX,Y2_new);
%             title(['$(y^{20}_{1},y^{20}_{2})$'],'fontsize',20,'interpreter','latex');          
%             set(gca,'FontSize',20);
%             subplot(2,1,2);
%             surf(TT1,XX,Lambda1_new);
%             hold on
%             %subplot(2,2,4);
%             surf(TT2,XX,Lambda2_old);
%             xlabel('t','FontSize',20);
%             ylabel('x','FontSize',20);
%             title(['$(\lambda^{20}_{1},\lambda^{19}_{2})$'],'fontsize',20,'interpreter','latex');          
%             set(gca,'FontSize',20);
%             print -depsc iteratesrand_5.eps
%             end
%             
%             Lambda2_old = Lambda2_new;             
        end

            rho = error_FVM(4)/error_FVM(3);
        %end
        % IT

%             semilogy(IT(1:end),error_FVM(1:end),'-*','LineWidth',2.0);
%             hold on
%         set(gca,'FontSize',20);
%         xlabel('iteration','FontSize',20);
%         ylabel('error','FontSize',20);
%     legend('k = 1', 'k = 20', 'k = 40', 'rand');
         %title('Iterative error with sine initial guess','FontSize',20);
        % error_FVM = [];
%          error_FVM_theo = [];
%          print -depsc convnumk.eps


    
         
         