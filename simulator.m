%% Data
clear all
close all
data = read("^GSPC 00-19.csv",17,17);
y0 = data.Close(1);
dates = data.Time;
N = length(dates);
numSims = 4;
hist = read("^GSPC 00-19.csv",16,16);
h0 = hist.Close(1);
hdates = hist.Time;
hN = length(hdates);

%% MJD
clf
mu = 0.001124725615659;
sigma = 0.005735804232002;
nu = -8.552016486024466e-04;
tau = 0.005735804232002;
lambda = 0.819895882519280;
figure(1)
title("MJD")
hold on
plot(dates, data.Close, 'k')
hold on 
for j = 1:numSims
    y = MJD_Ito(y0, mu, sigma, nu, tau, lambda, N);
    figure(1)
    plot(dates, y,':')
    hold on    
    figure(2)
    subplot(2,2,j);
    qqplot(data.Close,y)
    hold on
end

%% GBM
clf
mu = 0.106496177725054;
sigma = 0.104283431591000;

figure(1)
title(["GBM mu = ",num2str(mu),"sigma = ",num2str(sigma)])
hold on 
plot(dates, data.Close, 'k')
hold on
figure(2)
title("qqplot: GBM")
hold on 
for j = 1:numSims
    y = GBM_Euler(y0, sigma, mu, N);
    figure(1)
    plot(dates, y,'-.')
    hold on 
    figure(2)
    subplot(2,2,j);
    qqplot(data.Close,y)   
end 

%% Functions
function price = GBM_Euler(price, sigma, mu, N)
    dt = 1/N;
    pd = makedist('Normal');
    for i = 1:N-1
        price(i+1) = (1 + mu*dt + sigma*sqrt(dt)*random(pd))*price(i);
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

function period = yearrange(year_beg, year_end)
    date_beg = datestr(datenum(year_beg,1,1),'00yy-mm-dd');
    date_end = datestr(datenum(year_end,12,31),'00yy-mm-dd');
    period = timerange(date_beg, date_end);
end 

function data = read(filename, year_beg, year_end)
    table = table2timetable(readtable(filename));
    date_beg = datestr(datenum(year_beg,1,1),'00yy-mm-dd');
    date_end = datestr(datenum(year_end,12,31),'00yy-mm-dd');
    period = timerange(date_beg, date_end);
    data = flip(table(period,:),1);
end