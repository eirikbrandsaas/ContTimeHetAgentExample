%% Set parameters
clc
clear 
lambda = 0.02;
gamma = 2.0;
eta = 1.5;
p = 1.0;
rho = 0.018;
r = 0.01;

% grids
Na = 200;
amin=0.0;
amax=5.0;
agrid = linspace(amin,amax,Na)';
Da = agrid(2)-agrid(1);

w =0.01;
z = [1 3];
ygrid = w*z;
Ny=length(ygrid);

% create matrices for the states to make some stuff easier
aa = agrid * ones(1,Ny);
yy = ones(Na,1)*ygrid;

% set the utility functions and it's derivative
util = @(c,gamma) c.^(1-gamma)/(1-gamma);
uprime = @(c,gamma) c.^(-gamma);
uprimeinv = @(dV,gamma) dV.^(-1/gamma);

% next we find h as a function of c
f = @(h,eta) h.^(1-eta)/(1-eta);
fprime = @(h,eta) h.^(-eta);
fprimeinv = @(dV,eta) dV.^(-1/eta);

findh = @(c,gamma,eta) fprimeinv(p*uprime(c,gamma),eta);
bc = @(c,y,a) y - c - p*findh(c,gamma,eta) + r*a ;
cgrid = linspace(0,max(ygrid(2) + r*amax)*1.05,5000);

c0=zeros(Na,Ny);
% next, find out the total expenditure if we stay put
temp = @(c,inc,assets) bc(c,inc,assets);
for ia=1:Na
    for iy=1:2
        expenses = bc(cgrid,ygrid(iy),agrid(ia));
        [~,ic]=min(abs(expenses));
        
        inc = ygrid(iy);
        assets = agrid(ia);
        
        c0(ia,iy)=cgrid(ic);
        h0(ia,iy) = findh(c0(ia,iy),gamma,eta);
    end
end


%%
% numerical parameters
maxit = 10;
crit = 10^(-4);
Delta = 1000;

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
V0 = util(c0,gamma)/rho + f(h0,eta)/rho;
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
    cf(1,:) = min(c0(1,:)-0.000001,cf(1,:));
    cb(1,:) = min(c0(1,:)-0.000001,cb(1,:));
    hf = findh(cf,gamma,eta);
    hb = findh(cb,gamma,eta);
    
    
    adotf = bc(cf,yy,aa);
    adotb = bc(cb,yy,aa);
 
    Hf = util(cf,gamma) + f(hf,eta) + dVf.*adotf;
    Hb = util(cb,gamma) + f(hb,eta) + dVb.*adotb;
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
    h = findh(c,gamma,eta);
    adot = bc(c,yy,aa);
    u = util(c,gamma) + f(h,eta) ;
    
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
        diff
        break
    end
end
%V=Vnew;

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

subplot(2,2,3)
plot(agrid,g)
title("Distribution")

subplot(2,2,4)
plot(agrid,[c findh(c,gamma,eta)])
title("Distribution")


    
    
    
