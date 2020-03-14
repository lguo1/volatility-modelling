%% Initial Parameters and Data
clear
data = read("^GSPC 00-19.csv",16,16);
dates = data.Time;
R = diff(log(data.Close))';
N = length(dates);
dt = 1/N;
sigma = 0.104283431591000*sqrt(dt);
mu = 0.106496177725054*dt;
nu = .1*sigma;
lambda = 0.2;
k_bound = 5;
b_bounds = 1;
iterMax = 100;

%% Estimate
clf
[mu, sigma, nu, tau, lambda] = main(R, mu, sigma, nu, lambda, k_bound, b_bounds, iterMax);

%% Checks
price = MJD_Ito(data.Close(1), mu, sigma, nu, tau, lambda, N);
R_predicted = diff(log(price))';
plot(1:N-1, R, 'k')
hold on 
plot(1:N-1, R_predicted, 'b:')
sigmaE = sqrt(sigma^2+lambda*(nu^2+tau^2))
sigmaS = std(R_predicted)
sigmaA = std(R)
muE = mu +lambda*nu
muS = mean(R_predicted)
muA = mean(R)
title("log-return checks")
legend("actual log returns","simulated log returns")
xlabel("days")
ylabel("log-returns")


%% Functions
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

function [mu, sigma, nu, tau, lambda] = main(R, mu, sigma, nu, lambda, k_bound, b_bound, iterMax)
    N = length(R);
    R_bar = mean(R);
    for i = 1:iterMax
        [b, lambda_tilde, R_tilde, alpha_bar, sigma_square] = find_mins(R, N, R_bar, mu, sigma, nu, lambda, k_bound, b_bound);
        alpha_bar_frac = 1/alpha_bar -1; 
        mu = R_tilde - nu*alpha_bar_frac/b^2; %p10 (17) checked
        nu = b^2*(R_bar - R_tilde)/(b^2*lambda - alpha_bar_frac); %p10 (17) checked, lambda is lambda_tilde in paper
        sigma = sqrt(sigma_square);
        lambda = lambda_tilde;
    end 
    tau = b*sigma;
end

function [min_b, min_lambda_tilde, min_R_tilde, min_alpha_bar, min_sigma_square] = find_mins(R, N, R_bar, mu, sigma, nu, lambda, k_bound, b_bound)
    bs = linspace(0, b_bound, 20);
    old_min = inf;
    for i = 1:20
        [value, lambda_tilde, R_tilde, alpha_bar, sigma_tilde_square] = q_b_square(R, N, R_bar, mu, sigma, nu, bs(i), lambda, k_bound);
        if value < old_min
            old_min = value;
            min_b = bs(i);
            min_lambda_tilde = lambda_tilde;
            min_R_tilde = R_tilde;
            min_alpha_bar = alpha_bar;
            min_sigma_square = sigma_tilde_square;
        end
    end
end

function [value, lambda_tilde, R_tilde, alpha_bar, sigma_tilde_square] = q_b_square(R, N, R_bar, mu, sigma, nu, b, lambda, k_bound)
    alphas = zeros(1,N);
    Ns = zeros(1,N);
    logNBs = zeros(1,N);
    ks = 0:k_bound;
    for n = 1:N
        N_probs = pdfNs(R(n), mu, sigma, nu, b*sigma, lambda, ks);
        alphas(n) = alpha(b, N_probs, ks);
        logNBs(n) = logNB(b, N_probs, ks);
        Ns(n) = exptN(N_probs, ks);
    end 
    lambda_tilde = mean(Ns); % p7 (8)
    R_tilde = sum(alphas.*R)/sum(alphas); % p10 (18.5)
    alpha_bar = mean(alphas); % p10 (18.5)
    sigma_tilde_square = sum(alphas.*(R-R_tilde).^2)/N - (R_tilde-R_bar)^2/(b^2*lambda_tilde - (1 -alpha_bar)/alpha_bar);% p10 (18)
    value = log(sigma_tilde_square) + mean(logNBs);
end 

function expectation = logNB(b, N_probs, ks)
% p10 (19)
    expectation = sum(N_probs.*log(1+ks.*b^2));
end

function expectation = alpha(b, N_probs, ks)
% p8 (11)
     expectation = sum(N_probs.*(1./(1+ks.*b^2)));
end

function expectation = exptN(N_probs, ks)
% p8 (12) checked
    expectation = sum(N_probs.*ks);
end

function probs = pdfNs(Rn, mu, sigma, nu, tau, lambda, ks)
% p8 (12) checked
    probs = ks;
    for k = ks
        probs(k+1) = pdf('Normal', Rn, mu+k*nu, sqrt(sigma^2+k*tau^2))*lambda^k/factorial(k); 
    end
    probs = probs/sum(probs);
end 

function data = read(filename, year_beg, year_end)
    table = table2timetable(readtable(filename));
    date_beg = datestr(datenum(year_beg,1,1),'00yy-mm-dd');
    date_end = datestr(datenum(year_end,12,31),'00yy-mm-dd');
    period = timerange(date_beg, date_end);
    data = flip(table(period,:),1);
end
