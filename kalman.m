clear all

N=1000 %sampling timew

dt=0.001%del t

t=dt*(1:N)%vector for time

F=[1 dt

    0 1];

G=[-1/2*dt^2

    -dt];

H=[1 0];

Q=[0 0

    0 0]% no noise assumed

u=9.80665;

I=eye(2);

y0=100; %initial position

v0=0; %initial velocity



%state vectors

xt=zeros(2,N);

xt(:,1)=[y0;v0]%first state are assumed as initial position and velocity

%to predict true states are to be considered

for k=2:N

    xt(:,k)=F*xt(:,k-1)+G*u; 

end



%noisy measurement from the actual states

R=4;%error variance

v=sqrt(R)*randn(1,N) %random error in velocity measure,ent

z=H*xt+v;

x=zeros(2,N)

x(:,1)=[10 0]%guess values

P=[50 0

    0 0.01]

for k=2:N

    x(:,k)=F*x(:,k-1)+G*u;

    %covariance matrix prediction

    P=F*P*F'+Q;

    %kalman gain matrix

    K=P*H'/(H*P*H'+R);

    x(:,k)=x(:,k)+K*(z(k)-H*x(:,k));

    %update covariance matrix

    P=(I-K*H)*P;

end

%plotting the states

figure(1);

subplot(211);

plot(t,z,'g-',t,x(1,:),'b--','Linewidth',2);

hold on

plot(t,xt(1,:),'r:','Linewidth',1.5);

legend('Measured','Estimated','True');

subplot(212);

plot(t,x(2,:),'Linewidth',2);

hold on;

plot(t,xt(2,:),'r:','Linewidth',1.5);

legend('Estimated','True');

figure(2);

subplot(211);

plot(t,x(1,:)-xt(1,:),'b--','Linewidth',2);



legend('Measured','Estimated','True');

subplot(212);

plot(t,x(2,:)-xt(2,:),'Linewidth',2);



legend('Estimated','True')