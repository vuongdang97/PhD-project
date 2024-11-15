function [Y,Lambda] = CoupledWaveDir_FVM1(y0,y0t,z,yhat,lam_T,mu,bl,br,a,b,T0,T,n,m,alpha,gamma)
% WAVE solves the coupled wave equation with zero Dirichlet boundary conditions
%      [Y,Lambda]=CoupledWaveDir_FVM(y0,y0t,z,yhat,lam_T,mu,bl,br,a,b,T0,T,n,m,alpha,gamma) solves 
%      the coupled wave equation on the domain [a,b]x[T0,T]  with initial condition
%      y0 and y0t and final data lam(T) and z, zero Dirichlet boundary condition bl and br
%      using a finite volume method proposed in the paper of Laurence
%      Halpern, Martin J Gander and Frérédic Nataf
%      Alpha and gamma are two parameters
%      Remark: This function use the block structure proposed by Laurence
%      Halpern
%      n is number of mesh point in space
%      m is number of mesh point in time

Y = zeros(n,m); % coefficient of the basis with dim = n, m: time step
Lambda =  zeros(n,m); % coefficient of the dual solution with dim = n, m:time step
% Update the boundary data for Y and Lambda
Y(1,:) = bl; 
Y(n,:) = br;
Lambda(1,:) = bl;
Lambda(n,:) = br;
GrandVector = sparse(2*(n-2)*m,1); 
% Create a giant vector, we arrange the vector as follow [Y(2:n-1,1),Y(2:n-1,2),...,Y(2:n-1,M),Y(2:n-1,M+1),Lambda(2:n-1,1),Lambda(2:n-1,2),...,Lambda(2:n-1,M),Lambda(2:n-1,M+1)]
dx = (b-a)/(n-1);
dt = (T-T0)/(m-1);
x = (a:dx:b)';
t = (T0:dt:T);
s = dt^2/dx^2;
s1 = dt/dx;

% Create the matrix A
e = ones(n-2,1); 
A = spdiags([s*e 2*(1-s)*e s*e], -1:1, n-2, n-2);
C = spdiags([s1/2*e -s1*e s1/2*e], -1:1, n-2, n-2);

% Create the iterative matrix for the big system
B = sparse(2*(n-2)*m,2*(n-2)*m);
Yhat = sparse(2*(n-2)*m,1);
Z = sparse(2*(n-2)*m,1);
Yhat((n-2)+1:2*(n-2),1)= dt*y0t(2:n-1);
Yhat(1:n-2,1) = y0(2:n-1);
B((n-2)+1:2*(n-2),1:(n-2)) = eye(n-2)+dt/dx*C;
B((n-2)+1:2*(n-2),(m-1)*(n-2)+(n-2)+1:m*(n-2)+(n-2)) = 1/alpha*(dt)^2/2*eye(n-2);
for i = 1:m
    k1 = i; 
    k2 = i+m;
    if(k1 > 2)
       B((k1-1)*(n-2)+1:k1*(n-2),(k1-2)*(n-2)+m*(n-2)+1:(k1-1)*(n-2)+m*(n-2)) = (dt)^2/alpha*eye(n-2); 
       B((k1-1)*(n-2)+1:k1*(n-2),(k1-2)*(n-2)+1:(k1-1)*(n-2)) = A; 
       B((k1-1)*(n-2)+1:k1*(n-2),(k1-3)*(n-2)+1:(k1-2)*(n-2)) = -eye(n-2);
     end
    if(k2 <2*m-1)
       B((k2-1)*(n-2)+1:k2*(n-2),(k2-m)*(n-2)+1:(k2-m+1)*(n-2)) = -(dt)^2*eye(n-2);
       B((k2-1)*(n-2)+1:k2*(n-2),k2*(n-2)+1:(k2+1)*(n-2)) = A; 
       B((k2-1)*(n-2)+1:k2*(n-2),(k2+1)*(n-2)+1:(k2+2)*(n-2)) = -eye(n-2);
       % Create the yhat vector
       Yhat((k2-1)*(n-2)+1:k2*(n-2),1)= (dt)^2*yhat(2:n-1,i+1);
    end
    if(k2 == 2*m-1)
         B((k2-1)*(n-2)+1:k2*(n-2),(k2-m)*(n-2)+1:(k2-m+1)*(n-2)) = -(dt)^2/2*eye(n-2) - dt*gamma*eye(n-2);
         B((k2-1)*(n-2)+1:k2*(n-2),k2*(n-2)+1:(k2+1)*(n-2)) = dt/dx*C+eye(n-2); 
    end
end
% Create the final state vector
Z((2*m-1)*(n-2)+1:2*m*(n-2)) = lam_T(2:n-1);

Z((2*m-2)*(n-2)+1:(2*m-1)*(n-2)) = dt*z(2:n-1) + (dt)^2/2*yhat(2:n-1,m); % This is holds for the case gamma = 0;
%Z((2*m-2)*(n-2)+1:(2*m-1)*(n-2)) = gamma*dt*z(2:n-1) + (dt)^2/2*yhat(2:n-1,m); % This is holds for the case gamma \neq 0;


% Solve the problem numerically
M = speye(size(B))-B;
F = Yhat + Z;
GrandVector = M\F;

%Exact the discrete solution for big Vector

for i = 1:m
    Y(2:n-1,i) = GrandVector((i-1)*(n-2)+1:i*(n-2),1);
    Lambda(2:n-1,i) = GrandVector((i+m-1)*(n-2)+1:(i+m)*(n-2),1);
end