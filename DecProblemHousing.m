%% Set parameters
clc
clear 
lambda = 0.2;
gamma = 3.5;
eta = 1.3;
p = 8.0;
rho = 0.045;
r = 0.015;
d = 0.4;
hmin = 0.5
alpha = 0.2;
% grids
Na = 600;
amin=0.0;
amax=3.;
agrid = linspace(amin,amax,Na)';
Da = agrid(2)-agrid(1);

w =0.015;
z = [1 5];
ygrid = w*z;
Ny=length(ygrid);

% create matrices for the states to make some stuff easier
aa = agrid * ones(1,Ny);
yy = ones(Na,1)*ygrid;

% set the utility functions and it's derivative
util = @(c,gamma) c.^(1-gamma)/(1-gamma);
uprime = @(c,gamma) c.^(-gamma);
uprimeinv = @(dV,gamma) dV.^(-1/gamma);

func = @(h,eta) -alpha*exp(-eta*h) + exp(0) - r*p*h;
fprime =@(h,eta) alpha*eta*exp(-eta*h);
fprimeinv =@(dV,eta) -1/eta * log(r*p/(eta*alpha));

h = min(fprimeinv(r*p,eta), aa./(d*p));
h = h.*(h>hmin);
bc = @(c,f,y,a) y + f +r*a - c;
%%
% numerical parameters
maxit = 30;
crit = 10^(-6);
Delta = 100;

% preallocate some variables
dVf = zeros(Na,Ny);
dVb = zeros(Na,Ny);
dV0 = zeros(Na,Ny);
cf = zeros(Na,Ny);
cb = zeros(Na,Ny);

adotf = zeros(Na,Ny);
adotb = zeros(Na,Ny);
If = false(Na,Ny);
Ib = false(Na,Ny);
I0 = false(Na,Ny);

% initial guess (present value of staying put forever)
V0 = util(yy + func(h,eta) + r.*aa,gamma)/rho;
Vnew = V0;
tic
%%
for n=1:maxit
    %disp(sprintf('Starting iteration %d',n))
    V = Vnew;

    dVf(1:Na-1,:) = (V(2:Na,:) - V(1:Na-1,:))/Da;
    dVb(2:Na,:) = (V(2:Na,:) - V(1:Na-1,:))/Da;

    % End point corrections, only the first is important
    dVb(1,:) = uprime(ygrid + r*amin + func(h(1,:),eta),gamma);
    dVf(Na,:) = uprime(ygrid + r*amax + func(h(Na,:),eta),gamma); 
    
    cf = uprimeinv(dVf,gamma) ;
    cb = uprimeinv(dVb,gamma) ;
        
    adotf = bc(cf,func(h,eta),yy,aa);
    adotb = bc(cb,func(h,eta),yy,aa);
 
    Hf = util(cf,gamma)  + dVf.*adotf;
    Hb = util(cb,gamma)  + dVb.*adotb;
    
    c0 = yy +func(h,eta) +r.*aa;
    Ineither = (1-(adotf>0)) .* (1-(adotb<0));
    Iunique = (adotb<0).*(1-(adotf>0)) + (1-(adotb<0)).*(adotf>0);
    Iboth = (adotb<0).*(adotf>0);
    Ib = Iunique.*(adotb<0) + Iboth.*(Hb>=Hf);
    If = Iunique.*(adotf>0) + Iboth.*(Hf>=Hb);
    
    Ib(1,:) = false;
    I0 = Ineither;
    I0 = (1-If-Ib);
    
    
    dV0 = uprime(c0,gamma) ;%+ fprime(h0,eta);
      
    dVupwind = (dVf.*If + dVb.*Ib + dV0.*I0);
    c = uprimeinv(dVupwind,gamma);
    c = cf.*If + cb.*Ib + c0.*I0;
    adot = bc(c,func(h,eta),yy,aa);
    u = util(c,gamma) ;
    
    % Construct the A matrixIb
    Xvec = - Ib.*adotb/Da;
    Yvec = Ib.*adotb/Da - If.*adotf/Da - lambda;
    Zvec = If.*adotf/Da;
    
    
    A1block = spdiags(Yvec(:,1),0,Na,Na) + spdiags(Xvec(2:Na,1),-1,Na,Na) + spdiags([0;Zvec(1:Na-1,1)],1,Na,Na);
    A2block = spdiags(Yvec(:,2),0,Na,Na) + spdiags(Xvec(2:Na,2),-1,Na,Na) + spdiags([0;Zvec(1:Na-1,2)],1,Na,Na);
    lambdablock = lambda*speye(Na,Na);
    A = [A1block,lambdablock; lambdablock, A2block];
    
    B = (rho + 1/Delta)*speye(2*Na) - A;
    
    ustack = [u(:,1); u(:,2)];
    Vstack = [V(:,1); V(:,2)];
    
    b = ustack + Vstack/Delta;
    Vstack = B\b ;
    Vnew = [Vstack(1:Na), Vstack(Na+1:2*Na)];
    
    diff = max(max(abs(Vnew - V)));
    if diff<crit
        fprintf('Value function converged on iteration %d, distance %f \n',n,diff);
        break
    end
end
%V=Vnew;

It = If + Ib;


toc
%% Distribution - Old way of finding the stationary distribution
AT= A';
tempvec = zeros(Na*2,1);

% Need to hack one value so that it's not singular
ihack = 500;
tempvec(ihack) = 0.1;
row = zeros(1,Na*2);
row(ihack) = 1;
AT(ihack,:) = row;

gstack = AT\tempvec;
gmass = ones(1,2*Na)*gstack*Da;
gstack = gstack/gmass;

g = [gstack(1:Na), gstack(Na+1:2*Na)];

%% Distribution - iterative
% start with uniform
gstack = ones(2*Na,1);
gmass = ones(1,2*Na)*gstack*Da;
gstack = gstack./gmass;
gnew = gmass;
N = 500;
dt = 10;
for i=1:N
    gnew= (speye(2*Na) - AT*dt)\gstack;
    dist = max(abs(gnew-gstack));
    gstack = gnew;
    if dist < crit
        break
    end
end

g2 = [gstack(1:Na), gstack(Na+1:2*Na)];



%% Next lets plot the results

figure(1)
subplot(3,2,1)
plot(agrid,V)
title("Value function")

subplot(3,2,2)
plot(agrid,adot)
title("Savings adot")
hold on
plot(agrid,zeros(Na))
hold off

subplot(3,2,3)
plot(agrid,[c c-func(h,eta)])
title("Consumption")

subplot(3,2,4)
plot(agrid,h)
title("Housing")


subplot(3,2,5)
plot(agrid,g)
title("Distribution using old method")

subplot(3,2,6)
plot(agrid,g2)
title("Distribution using iterative method")

    
    
    
