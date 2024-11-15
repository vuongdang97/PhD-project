% This script is run the test code for theta.



%TT = 2:1:8;
TT = 0.2:0.2:30;
%theta = 1;
a = 0;
b = 1;
n = 51; % number mesh point in space
m = 51; % number mesh point in time
kk = 1/(b-a):1/(b-a):(n-1)/(b-a);
vec_theo_rho = [];
vec_num_rho = [];
eigenvalue1 = zeros(length(TT),length(kk));
eigenvalue2 = zeros(length(TT),length(kk));
eigenvalue3 = []
f = [];
% k = pi;
for l = 1:length(TT)
     T = TT(l)
     T1 = T/2;    
     for i = 1:length(kk)
           k = kk(i);
           xi = k*pi;
          % k = xi;
           nu1 = sqrt((xi^2+sqrt(xi^4+1))/2);
           nu2 = 1/(2*nu1);
           if(i == 1)
               f = [f 1/2*sinh(2*nu2*T1)+1/(4*(nu1^2+nu2^2))*(cosh(2*nu2*T1)+1)];
           end
           phi = 1/(4*(nu1^2+nu2^2))*(nu1^2*sinh(2*nu2*T1)^2-nu2^2*sin(2*nu1*T1)^2);
           d = (nu1^2*cosh(nu2*T1)^2+nu2^2*cos(nu1*T1)^2)/(nu1^2+nu2^2);
           lambda1 = -1/d^2*(sqrt(phi)+1/(4*(nu1^2+nu2^2))*(cosh(2*nu2*T1)-cos(2*nu1*T1)))^2;
           lambda2 = -1/d^2*(sqrt(phi)-1/(4*(nu1^2+nu2^2))*(cosh(2*nu2*T1)-cos(2*nu1*T1)))^2;
           eigenvalue1(l,i) = abs(lambda1);
     end
     eigenvalue3 = [eigenvalue3 max(abs(eigenvalue1(l,:)))];
end
% [TTT,kkk] = meshgrid(TT/2,kk);
% surf(kkk,TTT,eigenvalue1');
% xlabel('k','FontSize',30);
% ylabel('T_{1}','FontSize',30);
% zlabel('\rho','FontSize',30);
% set(gca,'FontSize',30);
plot(TT/2,eigenvalue3,TT/2,f,'LineWidth',2.0);
xlabel('T_{1}','FontSize',30);
ylabel('value','FontSize',30);
legend('\rho','f(T_{1})','FontSize',30);
set(gca,'FontSize',30);
   
