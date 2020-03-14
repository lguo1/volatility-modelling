clear
data = read("^GSPC 00-19.csv",16,16);
sigmas = GK(data, 252);
sigma = sigmas(end);
mu = log(data.Close(end))-log(data.Close(1));

%%
function vol = GK(data, bw)
    t = size(data,1);
    if (bw > t)
        error('The number of observations must be greater than the bandwidth.');
    end
    co = log(data.Close ./ data.Open);
    hl = log(data.High ./ data.Low);
    res = 0.5 * (hl .^ 2) - ((2 * log(2)) - 1) * (co .^ 2);
    param = sqrt(252 / bw);
    fun = @(x) param * sqrt(sum(x));
    win = get_rolling_windows(res,bw);
    win_len = length(win);
    win_dif = t - win_len;
    vol = NaN(t,1);
    for i = 1:win_len
        vol(i+win_dif) = fun(win{i});
    end
    vol(isnan(vol)) = [];
end

function data = fetch_data(varargin)

    persistent ip;

    if (isempty(ip))
        ip = inputParser();
        ip.addRequired('tkrs',@(x)validateattributes(x,{'cell','char'},{'vector','nonempty'}));
        ip.addRequired('date_beg',@(x)validateattributes(x,{'char'},{'nonempty','size',[1,NaN]}));
        ip.addRequired('date_end',@(x)validateattributes(x,{'char'},{'nonempty','size',[1,NaN]}));        
    end

    ip.parse(varargin{:});
    
    ip_res = ip.Results;
    date_beg = ip_res.date_beg;
    date_end = ip_res.date_end;
    
    try
        num_beg = datenum(date_beg,'yyyy-mm-dd');
        check = datestr(num_beg,'yyyy-mm-dd');

        if (~isequal(check,date_beg))
            error('Invalid end date specified.');
        end
    catch e
        rethrow(e);
    end

    try
        num_end = datenum(date_end,'yyyy-mm-dd');
        check = datestr(num_end,'yyyy-mm-dd');

        if (~isequal(check,date_end))
            error('Invalid end date specified.');
        end
    catch e
        rethrow(e);
    end

    if ((num_end - num_beg) < 30)
        error('The start date must be anterior to the end date by at least 30 days.');
    end

    data = cache_result(ip_res.tkrs,date_beg,date_end);

end

function data = fetch_data_internal(tkrs,date_beg,date_end)

    if (ischar(tkrs))
        tkrs = {tkrs};
    end
     
    tkrs_len = length(tkrs);
    data = cell(tkrs_len,1);
    
    date_orig = datenum('01-Jan-1970 00:00:00','dd-mmm-yyyy HH:MM:SS');
    date_beg = (datenum(date_beg,'yyyy-mm-dd') - date_orig) * 86400;
    date_end = (datenum(date_end,'yyyy-mm-dd') - date_orig) * 86400;

    url = ['https://query1.finance.yahoo.com/v8/finance/chart/%TICKER%?symbol=%TICKER%&period1=' num2str(date_beg) '&period2=' num2str(date_end) '&interval=1d'];

    bar = waitbar(0,'Fetching data from Yahoo! Finance...');
    
    try
        for i = 1:tkrs_len
            tkr = tkrs{i};
            
            tkr_data = webread(strrep(url,'%TICKER%',tkr));
            tkr_d = tkr_data.chart.result.timestamp;
            tkr_o = tkr_data.chart.result.indicators.quote.open;
            tkr_h = tkr_data.chart.result.indicators.quote.high;
            tkr_l = tkr_data.chart.result.indicators.quote.low;
            tkr_c = tkr_data.chart.result.indicators.quote.close;
            tkr_ac = tkr_data.chart.result.indicators.adjclose.adjclose;

            scl = tkr_ac ./ tkr_c;
            date = (tkr_d ./ 86400) + datenum(1970,1,1);
            open = tkr_o .* scl;
            high = tkr_h .* scl;              
            low = tkr_l .* scl;
            cls = tkr_c .* scl;
            ret = [NaN; diff(log(tkr_c))];

            data{i} = table(date,open,high,low,cls,ret,'VariableNames',{'Date' 'Open' 'High' 'Low' 'Close' 'Return'});

            waitbar((i / tkrs_len),bar);
        end
        
        close(bar);
        
        if (tkrs_len == 1)
            data = data{1};
        end
    catch e
        close(bar);
        rethrow(e);
    end

end

function varargout = cache_result(varargin)

    persistent cache;

    args = varargin;
    fun = @fetch_data_internal;
    now = cputime;
    
    key = [args {@fetch_data_internal,nargout}];
    key_inf = whos('key');
    key_sid = sprintf('s%.0f',key_inf.bytes);

    try
        pool = cache.(key_sid);

        for i = 1:length(pool.Inps)
            if (isequaln(key,pool.Inps{i}))
                varargout = pool.Outs{i};

                pool.Frqs(i) = pool.Frqs(i) + 1;
                pool.Lcnt = pool.Lcnt + 1;
                pool.Luse(i) = now;

                if (pool.Lcnt > cache.Cfg.ResFre)
                    [pool.Frqs,inds] = sort(pool.Frqs,'descend');

                    pool.Inps = pool.Inps(inds);
                    pool.Lcnt = 0;
                    pool.Luse = pool.Luse(inds);
                    pool.Outs = pool.Outs(inds);
                end

                cache.(key_sid) = pool;

                return;
            end
        end
    catch
        pool = struct('Frqs',{[]},'Inps',{{}},'Lcnt',{0},'Luse',{[]},'Outs',{{}});
    end

    if (~exist('varargout','var'))
        if (~isfield(cache,'Cfg'))
            cache.Cfg = struct();
            cache.Cfg.GrpSiz = 100;
            cache.Cfg.MaxSizKey = 100000000;
            cache.Cfg.MaxSizRes = 100000000;
            cache.Cfg.ResFre = 10;
        end

        [varargout{1:nargout}] = fun(varargin{:});
        var_inf = whos('varargout');

        if ((var_inf.bytes <= cache.Cfg.MaxSizRes) && (key_inf.bytes <= cache.Cfg.MaxSizKey))
            pool.Frqs(end+1) = 1;
            pool.Inps{end+1} = key;
            pool.Lcnt = 0;
            pool.Luse(end+1) = now;
            pool.Outs{end+1} = varargout;

            while (length(pool.Inps) > cache.Cfg.GrpSiz)
                [~,idx] = min(pool.Luse);

                pool.Luse(idx) = [];
                pool.Frqs(idx) = [];
                pool.Inps(idx) = [];
                pool.Outs(idx) = [];
            end

            cache.(key_sid) = pool;
        end
    end

end

% [INPUT]
% data = A numeric t-by-n matrix containing the time series.
% bw   = An integer representing the bandwidth (dimension) of each rolling window.
%
% [OUTPUT]
% win  = A vector of numeric bw-by-n matrices representing the rolling windows.
%
% [NOTES]
% If the number of observations is less than or equal to the specified bandwidth, a single rolling window containing all the observations is returned.

function win = get_rolling_windows(varargin)

    persistent ip;

    if (isempty(ip))
        ip = inputParser();
        ip.addRequired('data',@(x)validateattributes(x,{'numeric'},{'2d','nonempty'}));
        ip.addRequired('bw',@(x)validateattributes(x,{'numeric'},{'scalar','integer','real','finite','>=',2}));
    end

    ip.parse(varargin{:});
    ip_res = ip.Results;
    
    win = get_rolling_windows_internal(ip_res.data,ip_res.bw);

end

function win = get_rolling_windows_internal(data,bw)

    t = size(data,1);
    
    if (bw >= t)
        win = cell(1,1);
        win{1} = data;
        return;
    end

    lim = t - bw + 1;
    win = cell(lim,1);

    for i = 1:lim
        win{i} = data(i:bw+i-1,:);
    end

end

function data = read(filename, year_beg, year_end)
    table = table2timetable(readtable(filename));
    date_beg = datestr(datenum(year_beg,1,1),'00yy-mm-dd');
    date_end = datestr(datenum(year_end,12,31),'00yy-mm-dd');
    period = timerange(date_beg, date_end);
    data = flip(table(period,:),1);
end