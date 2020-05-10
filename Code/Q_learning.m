%% A
%Constants
n = 3;
m = 1;
l = (n+m)*(n+m+1)/2;

A = [0 1 0; 0 0 1; 0 0 0]; %uknown
B = [0; 0; 1]; %uknown
M = eye(n); %from the cost criterion

%Finding Z
%Choosing initial x
x(:,1) = [1; 2; 3];
Z= zeros(l,l);
u = zeros(1,l);

while rank(Z) < l
    for i = 1:10
       %selecting random integer input between -10 and 10
       u(i) = randi([-10 10],1,1);

       Z(i,:) = [x(1,i)^2 x(2,i)^2 x(3,i)^2 u(i)^2 2*x(1,i)*x(2,i) 2*x(1,i)*x(3,i) 2*x(1,i)*u(i) 2*x(2,i)*x(3,i) 2*x(2,i)*u(i) 2*x(3,i)*u(i)];
        
       %testing input for new x
       x(:,i+1) = A*x(:,i) + B*u(i);
    end
end

Zinv = inv(Z);
%Run the iterative algorithm to estimate H and K
ro = 1;
[Hest, Kest] = EstimateQ(ro,Zinv,x,u,n,m,l,M)

%% B

syms p11 p12 p13 p22 p23 p33

P = [p11 p12 p13
     p12 p22 p23
     p13 p23 p33];

%Ricatti P =
M + A'*P*A - A'*P*B*inv(ro+B'*P*B)*B'*P*A 

%From ricatti
p = [1 0 0
     0 2 0
     0 0 3]

%Expected Theoretical optimal gain
R = 1;
Kexpected = -inv(ro+B'*p*B)*B'*p*A

%Simulation on Matlab
simtime = 40;

umat = zeros(1,simtime);
umat(1:l) = u;
for i = (l+1):simtime
    umat(i) = Kest * x(:,i);
    x(:,i+1) = A*x(:,i) + B*umat(i);
end

%Plotting simulation
%States
figure;
hold on
plot(x(1,:))
plot(x(2,:))
plot(x(3,:))
title('States Simulation')
xlabel('Time samples')
legend('x1','x2','x3')
grid on
hold off

%Input
figure;
plot(umat)
title('Input u Simulation')
xlabel('Time Samples')
grid on

%Simulation on Simulink
usim = [[0:9]' u'];