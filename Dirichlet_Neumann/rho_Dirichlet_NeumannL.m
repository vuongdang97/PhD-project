function [lamp,lamm,rho] = rho_Dirichlet_NeumannL(xi,T,teta)
% This function   computes the theoretical convergence factor of
% Dirichlet Neumann method with toy model alpha = 1, gamma = 0, given T and no
% overlapping delta = 0
     global alpha gamma delta
     alpha = 1;
     gamma = 0;
     delta = 0;
     nu1 = sqrt((xi.^2 + sqrt(xi.^4+1))/2);
     nu2 = 1./(2*nu1);
     mnu=sqrt(nu1.^2+nu2.^2);
     C=cos(nu1*T);
     S=sin(nu1*T);
     Ch=cosh(nu2*T);
     Sh=sinh(nu2*T);
     u=(Ch.^2-C.^2)/2;
     v=nu1.^2.*Ch.^2+nu2.^2.*C.^2;
     %w=nu1.^2.*Sh.^2-nu2.^2.*S.^2;
     w=v-mnu.^2;
     l1=u.^2+(v.^2+w.^2)/2;
     dis=mnu.^4.*(u.^2+(v+w).^2/4);
     den= S.*C.*nu2.^2-  Sh.*Ch.*nu1.^2;
     den=den.*mnu.^2;
     lamp0= (l1+sqrt(dis))./den;
     lamm0=(l1-sqrt(dis))./den;
     lamp=(1-teta)+teta*lamp0;
     lamm=(1-teta)+teta*lamm0;
     A=[lamp;lamm];
     rho=max(abs(A));