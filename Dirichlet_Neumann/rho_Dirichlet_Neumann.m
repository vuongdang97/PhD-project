function [rho1,rho2] = rho_Dirichlet_Neumann(xi,theta,T)
% This function is write to compute the theoretical convergence factor of
% Dirichlet Neumann method with toy model alpha = 1, gamma = 0, given T and no
% overlapping delta = 0
% Test the computation with sine and cosine of matrices
     global alpha gamma delta
     alpha = 1;
     gamma = 0;
     delta = 0;
     nu1 = sqrt((xi^2 + sqrt(xi^4+1))/2);
     nu2 = 1/(2*nu1); 
     T1 = T/2;
     % Test the computation with sine and cosine of matrices
    A1 = zeros(2,2);
    A1(1,1) = sin(nu1*T1)*sinh(nu2*T1);
    A1(1,2) = -1/(nu1^2+nu2^2)*(nu1*cos(nu1*T1)*sinh(nu2*T1)-nu2*sin(nu1*T1)*cosh(nu2*T1));
    A1(2,1) = cos(nu1*T1)*cosh(nu2*T1);
    A1(2,2) = 1/(nu1^2+nu2^2)*(nu1*sin(nu1*T1)*cosh(nu2*T1)+nu2*cos(nu1*T1)*sinh(nu2*T1));
    
    A2 = zeros(2,2);
    A2(1,1) = cos(nu1*(T-T1))*cosh(nu2*(T-T1));
    A2(1,2) = 1/(nu1^2+nu2^2)*(nu1*sin(nu1*(T-T1))*cosh(nu2*(T-T1))+nu2*cos(nu1*(T-T1))*sinh(nu2*(T-T1)));
    A2(2,1) = -sin(nu1*(T-T1))*sinh(nu2*(T-T1));
    A2(2,2) = 1/(nu1^2+nu2^2)*(nu1*cos(nu1*T1)*sinh(nu2*T1)-nu2*sin(nu1*T1)*cosh(nu2*T1));
    
    B1 = zeros(2,2);
    B1(1,1) = nu1*cos(nu1*T1)*sinh(nu2*T1)+nu2*sin(nu1*T1)*cosh(nu2*T1);
    B1(1,2) = sin(nu1*T1)*sinh(nu2*T1);
    B1(2,1) = -(nu1*sin(nu1*T1)*cosh(nu2*T1)-nu2*cos(nu1*T1)*sinh(nu2*T1));
    B1(2,2) = cos(nu1*T1)*cosh(nu2*T1);
    
    B2 = zeros(2,2);
    B2(1,1) = (nu1*sin(nu1*(T-T1))*cosh(nu2*(T-T1))-nu2*cos(nu1*(T-T1))*sinh(nu2*(T-T1)));
    B2(1,2) = -cos(nu1*(T-T1))*cosh(nu2*(T-T1));
    B2(2,1) = nu1*cos(nu1*(T-T1))*sinh(nu2*(T-T1))+nu2*sin(nu1*(T-T1))*cosh(nu2*(T-T1));
    B2(2,2) = sin(nu1*(T-T1))*sinh(nu2*(T-T1));
    
    Theta = [theta 0; 0 theta];
    MDN = inv(A1)*A2*inv(B2)*B1;
    MDN = eye(2)-Theta*(eye(2)-MDN);
    [VDN,DDN] = eig(MDN);
    DDN1 = diag(DDN);
    rho1 = max(abs(DDN1));
    rho2 = min(abs(DDN1));
    