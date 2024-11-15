% This script is used to compute the convergence factor of Dirichlet Neumann method
% theoretical and numerical
clear all; 
close all; 
clc; 

TT = [100 200];
%theta11 = 1/394780:0.05:1/394780;
theta11 = 1:0.01:1;
a = 0; 
b = 1;
%Xi = [pi 2*pi 5*pi 9*pi];
n = 11; % number mesh point in space
m = 11; % number mesh point in time
Xi = pi/(b-a):pi/(b-a):(n-1)*pi/(b-a);
%Xi = [pi/(b-a) 20*pi/(b-a) 40*pi/(b-a) 100];
vec_theo_rho1 = [];
vec_theo_rho2 = [];
vec_num_rho = [];
vec_theo_rhoL = [];
eigenvalue1 = zeros(length(TT),length(theta11)); 
eigenvalue2 = zeros(length(TT),length(theta11));
temp = 1;
for l = 1:length(TT)
     T = TT(l)
    for j = 1:length(theta11)
        theta = theta11(j);
        for i = 1:length(Xi)
            xi = Xi(i);
            % Theoretical convergence factor
            [rho_theo1,rho_theo2] = rho_Dirichlet_Neumann(xi,theta,T);
            %rho_theoL = rho_Dirichlet_NeumannL(xi,T,theta);
            vec_theo_rho1 = [vec_theo_rho1 abs(rho_theo1)];
            %vec_theo_rho2 = [vec_theo_rho2 abs(rho_theo2)];
            %vec_theo_rhoL = [vec_theo_rhoL abs(rho_theoL)];
            % Numerical convergence factor
            rho_num =  num_rho_Dirichlet_Neumann(xi,theta,T,n,m,a,b);
            vec_num_rho = [vec_num_rho abs(rho_num)];
        end
       % eigenvalue1(l,j) =  max(abs(vec_theo_rho));
        %eigenvalue2(l,j) = max(abs(vec_num_rho));
        %vec_theo_rho = [];
        %vec_num_rho = [];
%         figure
         semilogy((b-a)*Xi/pi,vec_theo_rho1,'DisplayName', ['theory T = ',num2str(T)],'LineWidth',2.0);
         hold on
        % semilogy((b-a)*Xi/pi,vec_num_rho2,'DisplayName', ['numerics'],'LineWidth',2.0);
        % hold on
        %semilogy(Xi/pi,vec_theo_rhoL,'DisplayName', ['theoretical_L theta = ',num2str(theta)],'LineWidth',2.0);
        %hold on
            semilogy((b-a)*Xi/pi,vec_num_rho,'--','DisplayName', ['numerics T = ',num2str(T)],'LineWidth',2.0);
            hold on;
             legend show
           % title('Convergence factor of Dirichlet Neumann algorithm respect to frequency','FontSize',20);
            % semilogy(kk,vec_theo_rho,kk,vec_num_rho,'LineWidth',2.0);
             xlabel('k','FontSize',20);
             ylabel('Convergence factor','FontSize',20);
             set(gca,'FontSize',20)
        vec_theo_rho1 = [];
        vec_num_rho = [];
    end
%     semilogy(TT,eigenvalue1,'DisplayName',['theory   ',],'LineWidth',2.0);
%     hold on
%     semilogy(TT,eigenvalue2,'DisplayName',['numerics ', ],'LineWidth',2.0);
%     hold on
%     legend show
%     set(gca,'FontSize',20);
%     %title('Convergence factor respect to T','FontSize',20);
%     xlabel('T','FontSize',20);
%     ylabel('Convergence factor','FontSize',20);
    %eigenvalue1 = [];
    %eigenvalue2 = [];
end
% hold on
% semilogy(TT,ones(size(TT)),'DisplayName',['rho = 1'],'LineWidth',2.0);
% legend show
%     semilogy(TT,eigenvalue2,'LineWidth',2.0);
%     hold on
%     semilogy(TT,ones(size(TT)),'LineWidth',2.0);
%     legend show
%     set(gca,'FontSize',20);
%     title('Convergence factor respect to T','FontSize',20);
%     legend('numerics','rho= 1','FontSize',20);
%     xlabel('T','FontSize',20);
%     ylabel('Convergence factor','FontSize',20);


% [TTT,Theta] = meshgrid(TT,theta11);
%     subplot(2,1,1)
%     surf(TTT/2,Theta,eigenvalue1','DisplayName',['theory ']);
%     hold on
%     %surf(TTT,Theta,ones(length(TT),length(theta11))','DisplayName',['rho = 1 ']);
%     legend show
%     set(gca,'FontSize',20);
%     %title('Convergence factor respect to T','FontSize',20);
%     xlabel('T1','FontSize',20);
%     ylabel('theta','FontSize',20);
%     zlabel('Convergence factor','FontSize',20);
%     subplot(2,1,2)
%     surf(TTT/2,Theta,eigenvalue2','DisplayName',['numeric ' ]);
%     hold on
%     %surf(TTT,Theta,ones(length(TT),length(theta11))','DisplayName',['rho = 1 ']);
%     legend show
%     set(gca,'FontSize',20);
%     %title('Convergence factor respect to T','FontSize',20);
%     xlabel('T1','FontSize',20);
%     ylabel('theta','FontSize',20);
%     zlabel('Convergence factor','FontSize',20);
