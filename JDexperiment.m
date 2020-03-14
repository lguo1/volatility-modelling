%% Main
clf
N = 500;
mu = 0;
sigma = 1;
a = 0;
b = 2;
lambda = 2;
x = 1:N-1;
y = MJD_Ito(3, mu, sigma, a, b, lambda, N);
R = abs(diff(log(y)))';
figure(1)
plot(x, R, 'k')
hold on 
title("MJD Experiment")
hold on 
sigmas = zeros(N-1,1);
Ns = zeros(N-1,1);
ks = 0:16;
for n = 1:N-1
    N_probs = pdfNs(R(n), mu, sigma, a, b, lambda, ks);
    Ns(n) = exptN(N_probs, ks);
    sigmas(n) = sqrt(sigma^2+ b^2*Ns(n));
end
plot(x, sigmas, 'b')
%% Function
function price = MJD_Euler(price0, mu, sigma, nu, tau, lambda, N)
    price = zeros(N,1);
    price(1) = price0;
    W = makedist('Normal');
    J = makedist('Normal', 'mu', nu, 'sigma', tau);
    
    for t = 1:N-1
        Nt = poissrnd(lambda);
        jumps = 1;
        for k = 1:Nt
            jumps = jumps*exp(random(J));
        end
        price(t+1) = price0*exp(mu*t+sigma*random(W))*jumps;
    end
end

function price = MJD_Ito(price0, mu, sigma, nu, tau, lambda, N)
    price = zeros(N,1);
    price(1) = price0;
    W = makedist('Normal');
    J = makedist('Normal', 'mu', nu, 'sigma', tau);
    jumps = 1;
    for t = 1:N-1
        Nt = poissrnd(lambda);
        for k = 1:Nt
            jumps = jumps*exp(random(J));
        end
        price(t+1) = price0*exp(mu*t+sigma*random(W))*jumps;
    end
end

function expectation = exptN(N_probs, ks)
% p8 (12)
    expectation = sum(N_probs.*ks);
end

function probs = pdfNs(Rn, mu, sigma, nu, tau, lambda, ks)
% p8 (12)
    probs = ks;
    for k = ks
        probs(k+1) = pdf('Normal', Rn, mu+k*nu, sqrt(sigma^2+k*tau^2))*lambda^k/factorial(k); 
    end
    probs = probs/sum(probs);
end 
