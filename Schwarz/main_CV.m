% This script is used to compute the convergence factor of Schwarz method
% theoretical and numerical
clear all;
close all;
clc;

TT= 1:0.01:1;
theta11 = 0:0.01:1;
a = 0;
b = 1;
n = 51; % number mesh point in space
m = 51; % number mesh point in time
kk = 1/(b-a):1/(b-a):(n-1)/(b-a);
vec_theo_rho = [];
vec_num_rho = [];
eigenvalue1 = zeros(length(TT),length(theta11));
eigenvalue2 = zeros(length(TT),length(theta11));
for l = 1:length(TT)
     T = TT(l);
    for j = 1:length(theta11)
        theta = theta11(j);
       for i = 1:length(kk)
            xi = kk(i)*pi
            % Theoretical convergence factor
            rho_theo = rho_Schwarz(xi,theta,T);
            vec_theo_rho = [vec_theo_rho rho_theo];
            % Numerical convergence factor
            %rho_num = num_rho_Schwarz(xi,theta,T,n,m,a,b);
            %vec_num_rho = [vec_num_rho rho_num];
        end
        eigenvalue1(l,j) =  max(abs(vec_theo_rho));
        %eigenvalue2(l,j) = max(abs(vec_num_rho));
        vec_theo_rho = [];
        %vec_num_rho = [];
        %figure;
        %semilogy(kk,vec_theo_rho,'DisplayName', ['theory'],'LineWidth',2.0);
        %hold on
        %legend show
        %set(gca,'FontSize',20);
        %semilogy(kk,vec_num_rho,'--','DisplayName', ['numerics'],'LineWidth',2.0);
        %legend show
        %xlabel('k','FontSize',20);
        %ylabel('convergence factor','FontSize',20);
       % print -depsc Comparison_N=201.eps
        
        % %     hold on
        % %     pause();
        % %      vec_theo_rho = [];
        % %     vec_num_rho = [];
    end

end
     plot(theta11,eigenvalue1,'LineWidth',2.0);
     hold on
     plot(theta11(end-1),eigenvalue1(end-1),'o','LineWidth',3.0);
     hold on
     %plot(theta11,1.0+0.*eigenvalue1,'--','LineWidth',2.0);
%     hold on
%     plot(TT/2,ones(size(TT)),'LineWidth',2.0);
%     legend show
     set(gca,'FontSize',20);
%    % title('Convergence factor respect to T','FontSize',20);
%     legend('theory','numerics','rho = 1','FontSize',20);
     %xlabel('\theta','FontSize',20);
     %ylabel('Convergence factor','FontSize',20);
%     % eigenvalue1 = [];

% [TTT,Theta] = meshgrid(TT,theta11);
%     subplot(2,1,1)
%     surf(TTT,Theta,eigenvalue1','DisplayName',['theory ']);
%     hold on
%     %surf(TTT,Theta,ones(length(TT),length(theta11))','DisplayName',['rho = 1 ']);
%     legend show
%     set(gca,'FontSize',20);
%     %title('Convergence factor respect to T','FontSize',20);
%     xlabel('T','FontSize',20);
%     ylabel('theta','FontSize',20);
%     zlabel('Convergence factor','FontSize',20);
%     subplot(2,1,2)
%     surf(TTT,Theta,eigenvalue2','DisplayName',['numerics ' ]);
%     hold on
%     %surf(TTT,Theta,ones(length(TT),length(theta11))','DisplayName',['rho = 1 ']);
%     legend show
%     set(gca,'FontSize',20);
%     %title('Convergence factor respect to T','FontSize',20);
%     xlabel('T','FontSize',20);
%     ylabel('theta','FontSize',20);
%     zlabel('Convergence factor','FontSize',20);





%     plot(TT,eigenvalue1,TT,eigenvalue2,TT,ones(size(TT)),'LineWidth',2.0);
%     hold on
%     legend show
%     set(gca,'FontSize',20);
%     title('Convergence factor respect to T','FontSize',20);
%     legend('theory','numerics','rho = 1','FontSize',20);
%     xlabel('T','FontSize',20);
%     ylabel('Convergence factor','FontSize',20);
    % eigenvalue1 = [];
