% Load the data
load("X.mat");  

% a) Plot histograms of the samples with different numbers of bins
figure
subplot(2,2,1)
histogram(X, 10)
title('10 bins')

subplot(2,2,2)
histogram(X, 30)
title('30 bins')

subplot(2,2,3)
histogram(X, 100)
title('100 bins')

subplot(2,2,4)
histogram(X, 200)
title('200 bins')

sgtitle('Histograms of X with different number of bins')

% b) Approximate E[sin(X)] using Monte Carlo estimation
E = mean(sin(X));  % Monte Carlo approximation of the expectation

disp(['E[sin(X)] ≈ ', num2str(E)])

%%

Nsim = 10^3;        % number of simulated days
Tday = 720;        % minutes per day

% Distributions
arrivals = [1 2 3 4];  p_arr = [0.25 0.40 0.20 0.15];
ableT    = [2 3 4 5];  p_able = [0.30 0.28 0.25 0.17];
bakerT   = [3 4 5 6];  p_baker = [0.35 0.25 0.20 0.20];

% CDFs for sampling
arr_cdf = cumsum(p_arr);
abl_cdf = cumsum(p_able);
bak_cdf = cumsum(p_baker);

cars_Able = zeros(Nsim,1);
cars_Baker = zeros(Nsim,1);
wait_mean = zeros(Nsim,1);

for d = 1:Nsim
    tAble = 0; tBaker = 0; t = 0;  % time trackers
    total_wait = 0; cars = 0;
    a = 0; b = 0;
    
    while t < Tday
        % Next arrival
        u = rand; t = t + arrivals(find(u <= arr_cdf,1));
        if t >= Tday, break; end
        
        % Check which carhop is free first
        if tAble <= tBaker
            server = 1;   % Able
        else
            server = 2;   % Baker
        end
        start = max(t, min(tAble, tBaker));            
        total_wait = total_wait + (start - t);
        cars = cars + 1;
        
        % Service
        if server == 1
            u = rand; s = ableT(find(u <= abl_cdf,1));
            tAble = start + s; 
            a = a + 1;
        else
            u = rand; s = bakerT(find(u <= bak_cdf,1));
            tBaker = start + s; 
            b = b + 1;
        end
    end
    
    cars_Able(d) = a;
    cars_Baker(d) = b;
    wait_mean(d) = total_wait / cars;
end

% Averages over all simulated days
cars_Able  = mean(cars_Able);
cars_Baker = mean(cars_Baker);
mean_twait = mean(wait_mean);

disp(['cars_Able  = ', num2str(cars_Able)])
disp(['cars_Baker = ', num2str(cars_Baker)])
disp(['mean_twait = ', num2str(mean_twait), ' min'])



%%

clear; close all;

% Simulation of a service line with three sequential desks
% Service times: U[20,50], U[10,30], U[10,40]
% Inter-arrival times: N(35,20), rounded to positive integer
% (a) Count served customers during the simulation
% (b) Repeat and store average in avg_served

timestep  = 1;
finaltime = 500;
N_iter    = 100;       % number of simulations for part (b)

served_total = zeros(N_iter,1);

for iter = 1:N_iter

    % Initialization
    inter_arrival = max(1, round(normrnd(35,20)));   % first interarrival time
    job_time  = [inter_arrival 0 0 0 0 0];          % arrival + queue + 3 desks + output

    % Random times from uniform distributions with limits:
    % lower/upper for [arrival, queue, desk1, desk2, desk3, output]
    lower = [0 0 20 10 10 0];
    upper = [0 0 50 30 40 0];

    % Generate first service times for desks
    job_time(3) = irand(lower(3), upper(3), 1);
    job_time(4) = irand(lower(4), upper(4), 1);
    job_time(5) = irand(lower(5), upper(5), 1);

    ns  = length(job_time);
    s   = zeros(1, ns); % 0 idle, 1 occupied
    age = s;

    s_in = 1;
    s(1) = s_in;        % input always occupied
    s_out = 0;
    s(ns) = s_out;      % output always free
    Output = [];

    S = []; A = [];
    S(1,:) = s; A(1,:) = age;
    i = 1:(ns-1);

    for t = 1:timestep:finaltime

        occupied       = s(i) > 0;
        age(occupied)  = age(occupied) + timestep;
        job_done       = age(i) >= job_time(i);
        next_vacant    = s(i+1) == 0;
        next_vacant(1) = 1; % always accept arrivals
        move = find(occupied & job_done & next_vacant);

        if ~isempty(move)
            s(move+1)   = s(move+1) + 1;
            s(move)     = s(move) - 1;
            age(move+1) = 0;
            age(move)   = 0;

            % new service or interarrival times
            for k = move
                dst = k + 1;
                if dst >= 3 && dst <= 5
                    job_time(dst) = irand(lower(dst), upper(dst), 1);
                end
            end

            % new inter-arrival time when new job arrives
            if ismember(1, move)
                job_time(1) = max(1, round(normrnd(35,20)));
            end
        end

        s(1)  = s_in;   % keep occupied
        Output = [Output; s(ns)];
        s(ns) = s_out;  % always free

        S = [S; s];
        A = [A; age];
    end

    served_total(iter) = sum(Output);

    % Plot first iteration
    if iter == 1
        fprintf('Served customers in 1st simulation: %d\n', served_total(iter))
        io = find(Output==1);
        figure
        subplot(1,2,1)
        histogram(diff(io))
        title('Output time intervals')
        subplot(1,2,2)
        spy(S(1:400,:))
        title('States (first 40 steps)')
    end
end

avg_output = mean(served_total);
fprintf('Average served over %d simulations: %.2f\n', N_iter, avg_output);


%%
clc; close all; clearvars

N = 15;                  % max number of desks
topen = [8 16];          % open hours (8h)
n = 1000;                % simulations per desk count
L = 15;                  % max allowed expected waiters

% Arrival rates (per hour)
lam_a = 20; lam_b = 170;       % quiet, busy
meanIA_a = 1/lam_a; meanIA_b = 1/lam_b;

% Service times (h)
mu = 1.5/60; sigma = 1.3/60; lb = 0.5/60;

Nw_a = zeros(N,1); Nw_b = zeros(N,1);
meanIA_vec = [meanIA_a, meanIA_b]; 
for desks = 1:N
    Nwait_a = zeros(n,1); Nwait_b = zeros(n,1);

    for i = 1:n
        for mode = 1:2   % 1=quiet, 2=busy 
            meanIA = meanIA_vec(mode);
            tcur = topen(1); tarr = [];
            while tcur < topen(2)
                tint = exprnd(meanIA);
                tarr = [tarr tcur+tint]; tcur = tcur+tint;
            end
            narr = length(tarr)-1;
            tserv = max(lb, mu + sigma*randn(1,narr));
            twait = zeros(1,narr);
            free = repmat(topen(1),1,desks);
            for j = 1:narr
                [mint,mini] = min(free);
                twait(j) = max(tarr(j),mint)-tarr(j);
                free(mini) = tarr(j)+twait(j)+tserv(j);
            end
            if mode==1, Nwait_a(i)=sum(twait>0);
            else,        Nwait_b(i)=sum(twait>0); end
        end
    end

    Nw_a(desks)=mean(Nwait_a);
    Nw_b(desks)=mean(Nwait_b);
end

No_of_needed_desks_a = find(Nw_a<=L,1);
No_of_needed_desks_b = find(Nw_b<=L,1); 

fprintf('Quiet time desks: %d\nBusy time desks: %d\n',No_of_needed_desks_a, No_of_needed_desks_b)

%%
clear; close all;

% Parametets
N = 15;                 % max desks
T = 8;                  % hours (8:00-16:00) -> simulate on [0, T]
n = 1000;               % sims per desk count
L = 15;                 % expected max waiters

lam_a = 20; lam_b = 170;    % arrivals per hour (quiet/busy)

mu = 1.5/60;                 % service mean (h)
sigma = 1.3/60;              % service std  (h)
lb = 0.5/60;                 % service lower bound (h)

% (optional) common random numbers for variance reduction:
narr_a = poissrnd(lam_a*T, n, 1);
narr_b = poissrnd(lam_b*T, n, 1);

Nw_a = zeros(N,1); Nw_b = zeros(N,1);

for desks = 1:N
    % quiet
    wait_cnt = zeros(n,1);
    for i = 1:n
        narr = narr_a(i);
        if narr==0, continue; end
        tarr  = sort(rand(1,narr))*T;                  % arrival times in [0,T]
        tserv = max(lb, mu + sigma*randn(1, narr));    % truncated normal

        free = zeros(1, desks);                         % server next-free times
        w = 0;
        for j = 1:narr
            [mint, idx] = min(free);
            if tarr(j) < mint, w = w + 1; end
            start = max(tarr(j), mint);
            free(idx) = start + tserv(j);
        end
        wait_cnt(i) = w;
    end
    Nw_a(desks) = mean(wait_cnt);

    % busy 
    wait_cnt = zeros(n,1);
    for i = 1:n
        narr = narr_b(i);
        if narr==0, continue; end
        tarr  = sort(rand(1,narr))*T;
        tserv = max(lb, mu + sigma*randn(1, narr));

        free = zeros(1, desks);
        w = 0;
        for j = 1:narr
            [mint, idx] = min(free);
            if tarr(j) < mint, w = w + 1; end
            start = max(tarr(j), mint);
            free(idx) = start + tserv(j);
        end
        wait_cnt(i) = w;
    end
    Nw_b(desks) = mean(wait_cnt);
end

No_of_needed_desks_a = find(Nw_a <= L, 1, 'first'); if isempty(No_of_needed_desks_a), No_of_needed_desks_a = N; end
No_of_needed_desks_b = find(Nw_b <= L, 1, 'first'); if isempty(No_of_needed_desks_b), No_of_needed_desks_b = N; end

fprintf('Quiet time desks: %d\nBusy time desks: %d\n', No_of_needed_desks_a, No_of_needed_desks_b);

figure; hold on; grid on
plot(1:N, Nw_a, 'o-', 'DisplayName','Quiet (20/h)')
plot(1:N, Nw_b, 's-', 'DisplayName','Busy (170/h)')
yline(L, '--', 'DisplayName','Limit 15 waiting')
xlabel('Desks open'); ylabel('Mean passengers waiting')
legend('Location','northeast'); title(sprintf('n = %d sims/point', n))
hold off



%%
clear; close all;

N = 15; T = 8; n = 1000; L = 15;                   % desks, hours, sims, limit
lam_a = 20; lam_b = 170;                           % arrivals per hour
mu = 1.5/60; sigma = 1.3/60; lb = 0.5/60;          % service times (h)

Nw_a = zeros(N,1); Nw_b = zeros(N,1);
lam_vec = [lam_a lam_b];

for desks = 1:N
    w_a = zeros(n,1); w_b = zeros(n,1);
    for i = 1:n
        for mode = 1:2
            lam = lam_vec(mode);
            narr = poissrnd(lam*T);
            if narr==0, continue; end
            tarr = sort(rand(1,narr))*T;
            tserv = max(lb, mu + sigma*randn(1,narr));
            free = zeros(1,desks); w = 0;
            for j = 1:narr
                [mint,idx] = min(free);
                if tarr(j) < mint, w = w + 1; end
                free(idx) = max(tarr(j),mint) + tserv(j);
            end
            if mode==1, w_a(i)=w; else, w_b(i)=w; end
        end
    end
    Nw_a(desks)=mean(w_a);
    Nw_b(desks)=mean(w_b);
end

No_of_needed_desks_a = find(Nw_a<=L,1); if isempty(No_of_needed_desks_a),No_of_needed_desks_a=N;end
No_of_needed_desks_b = find(Nw_b<=L,1); if isempty(No_of_needed_desks_b),No_of_needed_desks_b=N;end

fprintf('Quiet time desks: %d\nBusy time desks: %d\n', No_of_needed_desks_a, No_of_needed_desks_b);

figure; hold on; grid on
plot(1:N,Nw_a,'o-','DisplayName','Quiet (20/h)')
plot(1:N,Nw_b,'s-','DisplayName','Busy (170/h)')
yline(L,'--','DisplayName','Limit 15 waiting')
xlabel('Desks open'); ylabel('Mean passengers waiting')
legend('Location','northeast')
title(sprintf('n = %d simulations per point', n))
hold off


%%
clear; close all; clearvars

% boundaries for sampling
xmin = -4; xmax = 4;
ymin = -4; ymax = 4;

% number of random points
N = 1e5;

% random samples
x1 = xmin + (xmax - xmin) * rand(N,1);
x2 = ymin + (ymax - ymin) * rand(N,1);

% accept if inside region (Z <= 0)
z = Z(x1, x2);
inside = z <= 0;

% area estimate using accept-reject
Ahat = mean(inside) * (xmax - xmin) * (ymax - ymin);

% scatter plot of accepted and rejected points
figure; hold on; axis equal; box on;
scatter(x1(~inside), x2(~inside), 6, [0.8 0.8 0.8], 'filled')
scatter(x1(inside), x2(inside), 6, [0 0.45 0.74], 'filled')
xlabel('x_1'); ylabel('x_2');
title(sprintf('Estimated area = %.4f', Ahat))
legend('Rejected','Accepted','Location','best')
hold off


%%

% Simulation of a service line with three sequential desks
% Service times: U[20,50], U[10,30], U[10,40]
% Inter-arrival times: N(35,20), rounded to positive integer
% (a) Count served customers during the simulation
% (b) Repeat and store average in avg_output

rng('default');

timestep  = 1;
finaltime = 500;
N_iter    = 1e3;          % tehtävän b) mukaan

served_total = zeros(N_iter,1);

for iter = 1:N_iter

    % --- init (pidetään samat rakenteet kuin sinulla) ---
    lower = [0 0 20 10 10 0];
    upper = [0 0 50 30 40 0];

    ns  = 6;              % [arrival, queue, d1, d2, d3, output]
    s   = zeros(1, ns);   % s(2) = jonon pituus, s(3:5) = 0/1 varattu
    age = zeros(1, ns);   % ikä: saapumiskello + asemien kulunut aika
    s(1) = 1;             % saapumiskello "käynnissä"

    % inter-arrival-timer ja palveluajat
    tIA    = max(1, round(normrnd(35,20)));
    tleft1 = 0; tleft2 = 0; tleft3 = 0;   % jäljellä olevat palveluajat

    % logit kuten sinulla
    Output = [];
    S = []; A = [];
    S(1,:) = s; A(1,:) = age;

    served = 0;

    for t = 1:timestep:finaltime

        % --- päivitys: kuluneet ajat ---
        age(1) = age(1) + timestep;       % saapumiskellon ikä
        if s(3)==1, age(3) = age(3)+1; end
        if s(4)==1, age(4) = age(4)+1; end
        if s(5)==1, age(5) = age(5)+1; end

        % --- mahdolliset valmistumiset (oikealta vasemmalle) ---
        % desk3 -> ulos
        s_out = 0;
        if s(5)==1
            tleft3 = tleft3 - 1;
            if tleft3 <= 0
                s(5) = 0;
                age(5) = 0;
                s_out = 1;          % yksi valmis tällä stepillä
                served = served + 1;
            end
        end

        % desk2 -> desk3
        if s(4)==1
            tleft2 = tleft2 - 1;
            if tleft2 <= 0 && s(5)==0
                s(4) = 0; age(4)=0;
                s(5) = 1; age(5)=0;
                tleft3 = irand(lower(5), upper(5), 1);   % U[10,40]
            end
        end

        % desk1 -> desk2
        if s(3)==1
            tleft1 = tleft1 - 1;
            if tleft1 <= 0 && s(4)==0
                s(3) = 0; age(3)=0;
                s(4) = 1; age(4)=0;
                tleft2 = irand(lower(4), upper(4), 1);   % U[10,30]
            end
        end

        % --- saapumiset (voi tulla useita samalla stepillä) ---
        tIA = tIA - 1;
        while tIA <= 0
            s(2) = s(2) + 1;         % jonon pituus kasvaa
            age(1) = 0;              % uusi saapuminen
            tIA = tIA + max(1, round(normrnd(35,20)));
        end

        % --- jonosta desk1:een ---
        if s(2) > 0 && s(3)==0
            s(2) = s(2) - 1;         % jonosta sisään
            s(3) = 1; age(3)=0;
            tleft1 = irand(lower(3), upper(3), 1);       % U[20,50]
        end

        % --- kirjanpito & ulostulo ---
        s(1) = 1;                % saapumiskello päällä
        Output = [Output; s_out];
        S = [S; s];
        A = [A; age];
    end

    served_total(iter) = sum(Output);

    if iter == 1
        fprintf('Served customers in 1st simulation: %d\n', served_total(iter))
        io = find(Output==1);
        figure
        subplot(1,2,1)
        if numel(io) > 1
            histogram(diff(io))
        else
            histogram(0)
        end
        title('Output time intervals')
        subplot(1,2,2)
        spy(S(1:min(40,size(S,1)),:))
        title('States (first 40 steps)')
    end
end

avg_output = mean(served_total);
fprintf('Average served over %d simulations: %.2f\n', N_iter, avg_output);

%%
% Simulation of a service line with three sequential desks
% Service times: U[20,50], U[10,30], U[10,40]
% Inter-arrival times: N(35,20), rounded to positive integer
% (a) Count served customers during the simulation
% (b) Repeat and store average in avg_output

rng('default');

timestep  = 1;
finaltime = 500;
N_iter    = 1000;

served_total = zeros(N_iter,1);

% helper: uniform integer in [a,b], inclusive
uuni = @(a,b) randi([a,b],1,1);

for iter = 1:N_iter

    % ---------- init ----------
    % s = [arrival, queue_len, d1_busy, d2_busy, d3_busy, output_flag]
    s   = zeros(1,6);
    age = zeros(1,6);
    s(1) = 1;                  % arrival clock "on"

    % remaining service times (minutes) for desks
    rem1 = 0; rem2 = 0; rem3 = 0;

    % queue length
    q = 0;

    % inter-arrival countdown (minutes to next arrival)
    ia = max(1, round(normrnd(35,20)));

    % logs (match your structure)
    Output = [];
    S = []; A = [];
    S(1,:) = s; A(1,:) = age;

    served = 0;

    % ---------- simulation ----------
    for t = 1:timestep:finaltime

        % advance ages
        age(1) = age(1) + 1;                % arrival clock age
        if s(3)==1, age(3) = age(3) + 1; end
        if s(4)==1, age(4) = age(4) + 1; end
        if s(5)==1, age(5) = age(5) + 1; end

        % decrement remaining service times if busy
        if s(3)==1, rem1 = rem1 - 1; end
        if s(4)==1, rem2 = rem2 - 1; end
        if s(5)==1, rem3 = rem3 - 1; end

        % arrivals: may be multiple in the same minute
        ia = ia - 1;
        while ia <= 0
            q = q + 1;
            age(1) = 0;                                % mark an arrival
            ia = ia + max(1, round(normrnd(35,20)));
        end

        % process completions & moves, right-to-left within the same minute
        out_flag = 0;

        % desk3 -> out
        if s(5)==1 && rem3 <= 0
            s(5) = 0; age(5) = 0;
            out_flag = 1;
            served = served + 1;
        end

        % desk2 -> desk3
        if s(4)==1 && rem2 <= 0 && s(5)==0
            s(4) = 0; age(4) = 0;
            s(5) = 1; age(5) = 0;
            rem3 = uuni(10,40);                           % U[10,40]
        end

        % desk1 -> desk2
        if s(3)==1 && rem1 <= 0 && s(4)==0
            s(3) = 0; age(3) = 0;
            s(4) = 1; age(4) = 0;
            rem2 = uuni(10,30);                           % U[10,30]
        end

        % queue -> desk1
        if q > 0 && s(3)==0
            q = q - 1;
            s(3) = 1; age(3) = 0;
            rem1 = uuni(20,50);                           % U[20,50]
        end

        % write state/output
        s(2) = q;                 % store visible queue length just like your vector
        s(1) = 1;                 % arrival clock stays "on"
        Output = [Output; out_flag];
        S = [S; s];
        A = [A; age];
    end

    served_total(iter) = sum(Output);

    % plots for the first run (match your style)
    if iter == 1
        fprintf('Served customers in 1st simulation: %d\n', served_total(iter))
        io = find(Output==1);
        figure
        subplot(1,2,1)
        if numel(io) > 1
            histogram(diff(io))
        else
            histogram(0)
        end
        title('Output time intervals')
        subplot(1,2,2)
        spy(S(1:min(40,size(S,1)),:))
        title('States (first 40 steps)')
    end
end

avg_output = mean(served_total);
fprintf('Average served over %d simulations: %.2f\n', N_iter, avg_output);
