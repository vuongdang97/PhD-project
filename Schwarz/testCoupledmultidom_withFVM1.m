% This script is created to test the time DDM for wave equation on multisubdomain.
% Particularly, we solve the Coupled Wave equation on two domains by finite
% volume method proposed in paper of Laurence Halpern, Martin Jacob Gander
% and Férédic Nataf.
% We simulate the solution at each iteration and after the termination of
% the algorithm
% Moreover, we compute the l2 error and determine the order of convergence
% Remark: The function we use is the coupledwavedir1 which organize the
% block matrix structure proposed by Laurence Halpern

clc;
close all;
clear all;




% Fixed parameter
a = 0;
b = 1;
T = 1;
T1 = T/2;
L = b-a; % the length of the interval (a,b)
M = 26;
N = 51; % number of mesh point in space
h = L/(N-1); % dx
%dt = h; % determine dt with fixed ratio
kk = 1:1:1;


delta = 0;
mu = 1;

itermax = 200;
test = 1;
TOL = 10^(-10);

alpha = 1;
number_subdomain = [2 4 8 16];
error_FVM = cell(length(number_subdomain),1);
error_control = cell(length(number_subdomain),1);
error_max = cell(length(number_subdomain),1);
IT = cell(length(number_subdomain),1);

for k = 1:length(N)
    n = N(k);
    dx = (b-a)/(n-1);
    m = M(k);
    for q = 1:length(number_subdomain)
        gamma = 0;
        sub = number_subdomain(q);
        %% Compute the solution on the full time domain
        x = (a:dx:b);
        dt = T/(16*(m-1));
        ratio = dt/h
        t = (0:dt:T);
        delta = 0;
        
        % Define the y hat solution
        yhat_full = zeros(n,length(t));
        z_full = zeros(n,1); %(gamma*y(x,T) - lambda_t(x,T))/gamma
        
        % Test for example y(x,t) = sin(pi*x)(t-T)
        %y0_full = -sin(pi*x)*T;
        %y0t_full = sin(pi*x);
        y0_full = zeros(n,1);
        y0t_full = zeros(n,1);
        
%          for i = 1:length(x)
% %             % Test for example y(x,t) = sin(pi*x)*(t-T)
%             % yhat_full(i,:) = (alpha*pi^4+1)*sin(pi*x(i))*(t-T);
% %             % yhat(i,:) = 0;
%          end
        
        % Test for example y(x,t) = sin(pi*x)*(t-T)
        %z_full = -alpha*pi^2*sin(pi*x);
        z_full = zeros(length(x),1);
        %z_full = z_full';
        % Define the exact solution y
        yexact = zeros(n,length(t));
        Lexact = zeros(n,length(t));
%         for i = 1:length(x)
%             % Exact solution for example sin(pi*x)*(t-T)
%             yexact(i,:) = sin(pi*x(i))*(t-T);
%             Lexact(i,:) = alpha*pi^2*sin(pi*x(i))*(t-T);
%         end
        lam_T_full = Lexact(:,end);
        bl_full = zeros(1,length(t));
        br_full = zeros(1,length(t));
        T0 = 0;
        gamma = 0;
        [Y,Lambda] = CoupledWaveDir_FVM1(y0_full,y0t_full,z_full,yhat_full,lam_T_full,mu,bl_full,br_full,a,b,T0,T,n,16*(m-1)+1,alpha,gamma);
        
        %%
        % Create the solution y, lambda
        Y_old = cell(sub,1);
        Y_new = cell(sub,1);
        Lambda_old = cell(sub,1);
        Lambda_new = cell(sub,1);
        
        
        
        % Create mesh point in time
        t_mesh = cell(sub+1,1);
        for l = 1:sub+1
            t_mesh{l} = (l-1)*T/sub;
        end
        
        % Create mesh in time
        mesh_t = cell(sub,1);
        for l = 1:sub
            if(l==sub)
                mesh_t{l} = (t_mesh{l}:dt:t_mesh{l+1});
            else
                mesh_t{l} = (t_mesh{l}:dt:t_mesh{l+1}+delta);
            end
        end
        
        % Create boundary data
        bl = cell(sub,1);
        br = cell(sub,1);
        for l = 1:sub
            bl{l} = zeros(1,length(mesh_t{l}));
            br{l} = zeros(1,length(mesh_t{l}));
        end
        
        % Create initial data for y
        y0 = cell(sub,1);
        y0t = cell(sub,1);
        
        % Create final data for lambda
        lam_T = cell(sub,1);
        z = cell(sub,1);
        
        % Create internal data y hat
        yhat = cell(sub,1);
        x = (a:dx:b);
        t = (0:dt:T);
        
        % Test for example y(x,t) = sin(pi*x)*(t-T)
        %z{sub} = (-alpha*pi^2*sin(pi*x))'; %fixed data
        z{sub} = zeros(n,1);
        
        s1 = dt/dx;
        
        % Create the matrix C
        e = ones(n-2,1);
        C = spdiags([s1/2*e -s1*e s1/2*e], -1:1, n-2, n-2);
        
        
        % Test for example y(x,t) = sin(pi*x)*(t-T)
        % fixed initial data
        %y0{1} = (-sin(pi*x)*T)';
        %y0t{1} = (sin(pi*x))';
        y0{1} = zeros(n,1);
        y0t{1} = zeros(n,1);
        for l = 1:sub
            for i = 1:length(x)
                % % % Test for example y(x,t) = sin(pi*x)*(t-T)
                %yhat{l}(i,:) = (alpha*pi^4+1)*sin(pi*x(i))*(mesh_t{l}-T);
                yhat{l}(i,:) = zeros(1,length(mesh_t{l}));
            end
        end
        % Define the exact solution y and lambda on multisubdomain
        y_exact = cell(sub,1);
        lambda_exact = cell(sub,1);
        
        for l = 1:sub
            for i = 1:length(x)
                % % Exact solution for example sin(pi*x)*(t-T)
                %y_exact{l}(i,:) = sin(pi*x(i))*(mesh_t{l}-T);
                %  lambda_exact{l}(i,:) = alpha*pi^2*sin(pi*x(i))*(mesh_t{l}-T);
                y_exact{l}(i,:) = Y(i,t_mesh{l}*sub/T+1:t_mesh{l+1}*sub/T);
                lambda_exact{l}(i,:) = Lambda(i,t_mesh{l}*sub/T+1:t_mesh{l+1}*sub/T);
            end
            % Create Y_old and Lambda_old
            Y_old{l} = zeros(n,length(mesh_t{l}));
            Lambda_old{l} = zeros(n,length(mesh_t{l}));
        end
        
        % Data lambda_t(M+1)
        lam_T{sub} = zeros(n,1); % the final value of lambda for domain 2
        iter = 0;
        IT{q} = [];
        rng('default')
        for l = 1:sub-1
            lam_T{l} = sin(k*pi*x);
            z{l} = zeros(1,length(x))';
            y0{l+1} = zeros(1,length(x));
            y0t{l+1} = zeros(1,length(x));
        end
        
        % Solving the time DDM system
        while(test > TOL && iter < itermax)
            test = 0;
            err_FVM = 0;
            err_max = 0;
            err_con = 0;
            % Solve the equation on each subdomain
            for l = 1:sub
                if(l == sub)
                    [Y_new{l},Lambda_new{l}] = CoupledWaveDir_FVM1(y0{l},y0t{l},z{l},yhat{l},lam_T{l},mu,bl{l},br{l},a,b,t_mesh{l},t_mesh{l+1},n,length(mesh_t{l}),alpha,gamma);
                    test = test + norm(Y_new{l}-Y_old{l}) + norm(Lambda_new{l}-Lambda_old{l});
                else
                    [Y_new{l},Lambda_new{l}] = CoupledWaveDir_FVM1(y0{l},y0t{l},z{l},yhat{l},lam_T{l},mu,bl{l},br{l},a,b,t_mesh{l},t_mesh{l+1}+delta,n,length(mesh_t{l}),alpha,gamma);
                    test = test + norm(Y_new{l}-Y_old{l}) + norm(Lambda_new{l}-Lambda_old{l});
                end
            end
            test = (dt*dx)^(1/2)*test;
            
            % Update new solution Y and Lambda and new data
            for l = 1:sub
                Y_old{l} = Y_new{l};
                Lambda_old{l} = Lambda_new{l};
                if(l == 1)
                    lam_T{l} = Lambda_old{l+1}(:,delta/dt+1);
                    z{l}(2:n-1) = Lambda_old{l+1}(2:n-1,delta/dt+2)/dt - 1/dt*(eye(size(C))+dt/dx*C)*Lambda_old{l+1}(2:n-1,delta/dt+1) - gamma*Y_old{l+1}(2:n-1,delta/dt+1) + dt/2*(Y_old{l+1}(2:n-1,delta/dt+1)-yhat{l}(2:n-1,end));
                    z{l} = -z{l};
                elseif(l == sub)
                    y0{l}(2:n-1) = Y_old{l-1}(2:n-1,end-delta/dt);
                    y0t{l}(2:n-1) = (Y_old{l-1}(2:n-1,end-delta/dt)-Y_old{l-1}(2:n-1,end-delta/dt-1))/dt + 1/dx*C*Y_old{l-1}(2:n-1,end-delta/dt) + dt/(2*alpha)*Lambda_old{l-1}(2:n-1,end-delta/dt);
                else
                    lam_T{l} = Lambda_old{l+1}(:,delta/dt+1);
                    z{l}(2:n-1) = Lambda_old{l+1}(2:n-1,delta/dt+2)/dt - 1/dt*(eye(size(C))+dt/dx*C)*Lambda_old{l+1}(2:n-1,delta/dt+1) - gamma*Y_old{l+1}(2:n-1,delta/dt+1) + dt/2*(Y_old{l+1}(2:n-1,delta/dt+1)-yhat{l}(2:n-1,end));
                    z{l} = -z{l};
                    y0{l}(2:n-1) = Y_old{l-1}(2:n-1,end-delta/dt);
                    y0t{l}(2:n-1) = (Y_old{l-1}(2:n-1,end-delta/dt)-Y_old{l-1}(2:n-1,end-delta/dt-1))/dt + 1/dx*C*Y_old{l-1}(2:n-1,end-delta/dt) + dt/(2*alpha)*Lambda_old{l-1}(2:n-1,end-delta/dt);
                end
                err_FVM = err_FVM + (dt*dx)^(1/2)*norm(Y_old{l}-y_exact{l});
                err_max = err_max + max(max(abs(Y_old{l}-y_exact{l})));
                err_con = err_con + (dt*dx)^(1/2)*(norm(Lambda_old{l}-lambda_exact{l}));
            end
            % Increase the iteration
            iter = iter + 1;
            IT{q} = [IT{q} iter];
            error_FVM{q} = [error_FVM{q} err_FVM];
            error_max{q} = [error_max{q} err_max];
            error_control{q} = [error_control{q} err_con];
%             % Plot the solution on domain 1
%             for l = 1:sub
%                 for j=1:length(mesh_t{l})
%                     % Plot the solution Y
%                     subplot(2,1,1);
%                     plot(x,Y_new{l}(:,j),'-x b',x,y_exact{l}(:,j));
%                     axis([0 1 -2 1]);
%                     % Plot the control u
%                     subplot(2,1,2);
%                     plot(x,Lambda_new{l}(:,j),'-x b',x,lambda_exact{l}(:,j));
%                     legend('discrete solution', 'exact solution');
%                     %axis([0 1 -5 1]);
%                     title(['time = ', num2str(mesh_t{l}(j)),'it = ', num2str(iter)]);
%                     pause(0.1)
%                 end
%                 pause
%             end
        end
        
        % Renew the iteration and the stopping criteria
        test = 1;
        iter = 0;
        
        error = (dt*dx)^(1/2)*norm(Y-yexact);
        semilogy(IT{q}(1:end),error_FVM{q}(1:end),'-*');
        %legend('sub = 2', 'sub = 4','sub = 8', 'sub = 16','FontSize',30);
        set(gca,'FontSize',30);
        xlabel('iteration','FontSize',30);
        ylabel('l^2 error','FontSize',30);
        %title('Iterative error with l^2 norm','FontSize',30);
        hold on
    end
end

% semilogy(IT(1:end),error_FVM(1:end),'-*');
% legend('error');

% % Plot the order of convergence
% loglog(N,N.^(-2),N,error_FVM{1},'-x',N,error_max{1},'-+',N,error_control{1},'-o');
% legend({'2x','error','infinity error','control error'},'FontSize',20);
% title('The convergence plot with dt = dx','FontSize',20);