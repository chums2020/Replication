clear all;
close all;
%%
% Parameter values 
alpha = (1/3)*0.85;
gamma = (2/3)*0.85;
beta = 0.96 ;
delta = 0.08;
ce = 1
cf = 0;
lamda = 1/10;

L = 1; %Aggregate supply of labor is inelastic and fixed at 1
s1 = 1
%% 
%Distribution of firm's productivity, s
load establishment_dist.txt;
datazupper = establishment_dist(:,1);
datahs = establishment_dist(:,2);
datazupper = [1; datazupper]
datahs = [0; datahs]
dataHs = cumsum(datahs)


plot(log(datazupper),dataHs, 'o');
hold on;

n_grid = logspace(0,4, 100)  % log-spaced grid with 100 points on the interval [1,10000]
n_cdf = interp1(datazupper,dataHs,n_grid,'linear'); % CDF of s is linear interpolation across the data points
plot(log(n_grid), n_cdf)
xlabel('Number of employees (log scale)')
ylabel('Cumulative Distribution of Establishments')
legend({'Model' 'Data'});
%% 
%Use number of employees, n, to back out firm's productivity, s
for i = 1:length(n_grid)
s(i) = s1*(n_grid(i)/n_grid(1))^(1-gamma- alpha)  %s1 =1 from normalization
end
% s ranges between [1, 3.98]

%density
h(1) = n_cdf(1)
for i = 2:length(n_cdf)
h(i) = n_cdf(i) - n_cdf(i-1)
end

%% 
%tau_t is a vector of different levels of distortion (Uncorrelated case)
%tau_s is the corresponding subsidy s.t. net effect on SS capital is zero
tau_t = [0.1, 0.2, 0.3, 0.4]
tau_s = []
for i = 1:length(tau_t)
tau_s(i) = 1-(2-(1-tau_t(i))^(1/(1-alpha-gamma)))^(1-alpha-gamma) 
end

TAU = [0,0; tau_t',tau_s']
%tau_s = [-0.0632, -0.0898, -0.1017, -0.1068]
%% 
%Case: No distortion, tau = TAU(1,:)
r = 1/beta - (1-delta); %from consumer optimization
rho = (1-lamda)/(1+r-delta)
tau = TAU(1,:) 

syms w;
for i = 1:length(s)
    for j = 1:length(tau)
       k(i,j) = (alpha/r)^((1-gamma)/(1-gamma-alpha))*(gamma/w)^(gamma/(1-gamma-alpha))*(s(i)*(1-tau(j)))^(1/(1-gamma-alpha))
       n(i,j) = ((1-tau(j))*s(i)*gamma/w)^(1/(1-gamma))*k(i,j)^(alpha/(1-gamma))
       pi(i,j) = (1-tau(j))*s(i)*k(i,j)^alpha*n(i,j)^gamma - w*n(i,j) - r*k(i,j) -cf
       W(i,j) = pi(i,j)/(1-rho)
    end
end

%We = (1/2)*h*(W(:,1)+W(:,2))
%Use secant method to find wss such that We == ce
%Verify that all firms have non-negative discounted profit

w0=1.8
w1=1.85
iterw = w0
iterw = [iterw;w1]
iterfw = double(subs((1/2)*h*(W(:,1)+W(:,2))-1,w,w0))
iterfw = [iterfw;double(subs((1/2)*h*(W(:,1)+W(:,2))-1,w,w1))]
error = 1
while error >=1e-5  
w2=w1-(w1-w0)*double(subs((1/2)*h*(W(:,1)+W(:,2))-1,w,w1))/double((subs((1/2)*h*(W(:,1)+W(:,2))-1,w,w1))-double(subs((1/2)*h*(W(:,1)+W(:,2))-1,w,w0)))
error = abs(double(subs((1/2)*h*(W(:,1)+W(:,2))-1,w,w2)))
iterw=[iterw;w2] 
iterfw = [iterfw;double(subs((1/2)*h*(W(:,1)+W(:,2))-1,w,w2))]
w0 = w1
w1= w2
end
%output wage in steady state, wss, and check free entry condition is
%satisfied
wss = iterw(length(iterw))
We = iterfw(length(iterfw)) + ce
%Check if all firms' discounted present values of profit is positive; otherwise,
%some firms are not producing, and they should be excluded when calculating
%We
Wss = double(subs(W,w, wss))
[nrows, ncols] = size(Wss);
for r = 1:nrows
    for c = 1:ncols
        if Wss(r,c) < 0
            Positive(r,c)  = 1 
        else Positive(r,c) = 0
        end
    end
end
sum(Positive)
%% 
%Case: No distortion, tau = TAU(1,:)
%Stead state aggregate
nss = double(subs(n,w, wss))
kss = double(subs(k,w, wss))
for i = 1:length(s)
    for j = 1:length(tau)
        yss(i,j) = s(i)*kss(i,j)^alpha*nss(i,j)^gamma
    end
end
E(1) = sqrt(2*lamda/(h*(nss(:,1)+nss(:,2))))
for i = 1:length(h)
    mu(i,1) = E(1)*(1/lamda)*h(i)/2
    mu(i,2) = E(1)*(1/lamda)*h(i)/2
end
Kss(1) = (h/2)*(kss(:,1) + kss(:,2));
Yss(1) = (h/2)*(yss(:,1) + yss(:,2));
TFP(1) = Yss(1)/(Kss(1)^gamma) ; %Note aggregate labor is fixed at 1
%% 
%Case: tax rate is 10%, tau = TAU(2,:)
%Solve wage in steady state, wss, from free entry condition, We - cf =0
%where ce=1

r = (1/beta)-(1-delta)
tau = TAU(2,:)
clearvars k n pi W w0 w1 w2 iterw  iterfw error wss We Wss
%express firm's policies, k, n, pi, W, in terms of wage, w
syms w;
for i = 1:length(s)
    for j = 1:length(tau)
       k(i,j) = (alpha/r)^((1-gamma)/(1-gamma-alpha))*(gamma/w)^(gamma/(1-gamma-alpha))*(s(i)*(1-tau(j)))^(1/(1-gamma-alpha))
       n(i,j) = ((1-tau(j))*s(i)*gamma/w)^(1/(1-gamma))*k(i,j)^(alpha/(1-gamma))
       pi(i,j) = (1-tau(j))*s(i)*k(i,j)^alpha*n(i,j)^gamma - w*n(i,j) - r*k(i,j) -cf
       W(i,j) = pi(i,j)/(1-rho)
    end
end

syms w;

%We = (1/2)*h*(W(:,1)+W(:,2))
%Use secant method to find wss such that We == ce
%Verify that all firms have non-negative discounted profit

w0=1.8
w1=1.85
iterw = w0
iterw = [iterw;w1]
iterfw = double(subs((1/2)*h*(W(:,1)+W(:,2))-1,w,w0))
iterfw = [iterfw;double(subs((1/2)*h*(W(:,1)+W(:,2))-1,w,w1))]
error = 1
while error >=1e-5  
w2=w1-(w1-w0)*double(subs((1/2)*h*(W(:,1)+W(:,2))-1,w,w1))/double((subs((1/2)*h*(W(:,1)+W(:,2))-1,w,w1))-double(subs((1/2)*h*(W(:,1)+W(:,2))-1,w,w0)))
error = abs(double(subs((1/2)*h*(W(:,1)+W(:,2))-1,w,w2)))
iterw=[iterw;w2] 
iterfw = [iterfw;double(subs((1/2)*h*(W(:,1)+W(:,2))-1,w,w2))]
w0 = w1
w1= w2
end
%output wage in steady state, wss, and check free entry condition is
%satisfied
wss = iterw(length(iterw))
We = iterfw(length(iterfw)) + ce
%Check if all firms' discounted present values of profit is positive; otherwise,
%some firms are not producing, and they should be excluded when calculating
%We
Wss = double(subs(W,w, wss))
[nrows, ncols] = size(Wss);
for r = 1:nrows
    for c = 1:ncols
        if Wss(r,c) < 0
            Positive(r,c)  = 1 
        else Positive(r,c) = 0
        end
    end
end
sum(Positive)

%% 
%Case: tax rate is 10%, tau = TAU(2, :)
%Stead state aggregate
nss = double(subs(n,w, wss))
kss = double(subs(k,w, wss))
for i = 1:length(s)
    for j = 1:length(tau)
        yss(i,j) = s(i)*kss(i,j)^alpha*nss(i,j)^gamma
    end
end
E(2) = sqrt(2*lamda/(h*(nss(:,1)+nss(:,2))))
for i = 1:length(h)
    mu(i,1) = E(2)*(1/lamda)*h(i)/2
    mu(i,2) = E(2)*(1/lamda)*h(i)/2
end
Kss(2) = (h/2)*(kss(:,1) + kss(:,2));
Yss(2) = (h/2)*(yss(:,1) + yss(:,2));
TFP(2) = Yss(2)/(Kss(2)^gamma) ; %Note aggregate labor is fixed at 1
Ys_ss(2) = (h/2)*( yss(:,2)); %output share of firms that are receiving a subsidy;
S(2) = h*yss(:,2)*(-tau(2))

%% 
%Effect of distortion -- uncorrelated case
Relative_TFP = TFP(2)/TFP(1);
Relative_Y = Yss(2)/Yss(1);
Relative_E = E(2)/E(1);
YsOverY = Ys_ss(2)/Yss(2);
SOverY = S(2)/Yss(2);

%% 


