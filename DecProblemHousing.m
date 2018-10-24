%% Set parameters
clc
clear 
lambda = 0.02;
gamma = 3.5;
eta = 0.5;
chi = 1.0;
rho = 0.018;
r = 0.01;

% grids
Na = 550;
amin=0.0;
amax=5.0;
agrid = linspace(amin,amax,Na)';
Da = agrid(2)-agrid(1);

w =0.01;
z = [1 3];
ygrid = w*z;
Ny=length(ygrid);

p=1.0;
% create matrices for the states to make some stuff easier
aa = agrid * ones(1,Ny);
yy = ones(Na,1)*ygrid;

% set the utility functions and it's derivative
util = @(c,gamma) c.^(1-gamma)/(1-gamma);
uprime = @(c,gamma) c.^(-gamma);
uprimeinv = @(dV,gamma) dV.^(-1/gamma);

f = @(h,eta) chi*h.^(1-eta)/(1-eta);
fprime = @(h,eta) chi*h.^(-eta);
fprimeinv = @(dV,eta) 1/chi*dV.^(-1/eta);
bc = @(c,y,a) y - c + r*a ;


% numerical parameters
maxit = 1;
crit = 10^(-6);
Delta = 1000;

% preallocate some variables
dVf = zeros(Na,Ny);
dVb = zeros(Na,Ny);
dV0 = zeros(Na,Ny);
cf = zeros(Na,Ny);
cb = zeros(Na,Ny);
c0 = zeros(Na,Ny);
adotf = zeros(Na,Ny);
adotb = zeros(Na,Ny);
If = false(Na,Ny);
Ib = false(Na,Ny);
I0 = false(Na,Ny);
V   = zeros(Na,Ny);

% initial guess (present value of staying put forever)
V0 =util(r.*aa + yy,gamma)/rho*0.1;
Vnew = V0;

%%
for n=1:maxit
    %disp(sprintf('Starting iteration %d',n))
    V = Vnew;
   
    dVf(1:Na-1,:) = (V(2:Na,:) - V(1:Na-1,:))/Da;
    dVb(2:Na,:) = (V(2:Na,:) - V(1:Na-1,:))/Da;

    % End point corrections, only the first is important
    dVb(1,:) = uprime(ygrid + r*amin,gamma);
    dVf(Na,:) = uprime(ygrid + r*amax,gamma); 

    cf = uprimeinv(dVf,gamma) ;
    cb = uprimeinv(dVb,gamma) ;

    hf = fprimeinv(p*dVf,eta);
    hb = fprimeinv(p*dVb,eta);

    adotf = bc(cf+p*hf,yy,aa);
    adotb = bc(cb+p*hb,yy,aa);

    c0 = yy + r*aa;
    dV0 = uprime(c0,gamma);


   
    
    Hf = util(cf,gamma) + f(hf,eta) + dVf.*adotf
    Hb = util(cb,gamma) + f(hb,eta) + dVb.*adotf
    Ineither = (1-(adotf>0)) .* (1-(adotb<0));
    Iunique = (adotb<0).*(1-(adotf>0)) + (1-(adotb<0)).*(adotf>0);
    Iboth = (adotb<0).*(adotf>0);
    Ib = Iunique.*(adotb<0) + Iboth.*(Hb>=Hf);
    If = Iunique.*(adotf>0) + Iboth.*(Hf>=Hb);
    I0 = Ineither;
    I0 = (1-If-Ib);
    
   
    dVupwind = (dVf.*If + dVb.*Ib + dV0.*I0);
    c = uprimeinv(dVupwind,gamma);
    h = fprimeinv(p*dVupwind,eta);
    adot = bc(c + p*h,yy ,aa);
    u = util(c,gamma) + f(h,eta);
    
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
    disp([n,diff])
    if diff<crit
        disp(sprintf('Value function converged on iteration %d',n));
        break
    end
end
V=Vnew;

It = If + Ib;



%% Distribution
AT= A';
tempvec = zeros(Na*2,1);

% Need to hack one value so that it's not singular
ihack = 1;
tempvec(ihack) = 0.1;
row = zeros(1,Na*2);
row(ihack) = 1;
AT(ihack,:) = row;

gstack = AT\tempvec;
gmass = ones(1,2*Na)*gstack*Da;
gstack = gstack/gmass;

g = [gstack(1:Na), gstack(Na+1:2*Na)];
%% Next lets plot the results

figure(1)
subplot(2,2,1)
plot(agrid,V)
title("Value function")

subplot(2,2,2)
plot(agrid,adot)
title("Savings adot")
hold on
plot(agrid,zeros(Na))
hold off

subplot(2,4,5:6)
plot(agrid,g)
title("Distribution")



    
    
    
