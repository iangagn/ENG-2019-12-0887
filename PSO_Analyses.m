%% ANALYSIS CONFIG
clear; clc;

nb_runs = 1000;
fun = @sphere;

% Boundaries
lb = -5.12;
ub = abs(lb);

% Minimum value
minv = 0;

% Initial positions
initial_pos = [];

% Chaotic maps
names = {'Chebyshev', 'Circle', 'Gauss', 'Logistic', 'Sine', 'Tent'};
cindex = [1, 2, 3, 5, 7, 10];

nb_maps = numel(names);

%% RUN

for c_index = 1 : nb_maps
    
    chaos_index = cindex(c_index);
    function_name = func2str(fun);
    
    chaotic_pso_config = struct('fun',fun,          ...
        'nb_dim',              2,                   ...
        'initial_positions',   initial_pos,         ...
        'chaos_index',         chaos_index,         ...
        'lb',                  lb,                  ...
        'ub',                  ub,                  ...
        'w',                   1.5,                 ...
        'c1',                  0.7,                 ...
        'c2',                  0.9,                 ...
        'nb_particles',        40,                  ...
        'max_iter',            4000,                ...
        'known_best_fitness',  minv,                ...
        'tol',                 1e-2,                ...
        'positions_hist_flag', false,                ...
        'diversity_hist_flag', true);
    
    pso_config = struct('fun',fun,          ...
        'nb_dim',              2,           ...
        'initial_positions',   initial_pos, ...
        'lb',                  lb,          ...
        'ub',                  ub,          ...
        'w',                   1.5,         ...
        'c1',                  0.7,         ...
        'c2',                  0.9,         ...
        'nb_particles',        40,          ...
        'max_iter',            4000,        ...
        'known_best_fitness',  minv,        ...
        'tol',                 1e-2,        ...
        'positions_hist_flag', false,        ...
        'diversity_hist_flag', true);
    
    chaotic_results.(names{c_index}).(function_name) = cell(nb_runs, 1);
    results.(function_name) = cell(nb_runs, 1);
    for i = 1 : nb_runs
        chaotic_results.(names{c_index}).(function_name){i} = chaotic_pso(chaotic_pso_config);
        if c_index == nb_maps
            results.(function_name){i} = pso(pso_config);
        end
        disp(['Running ' function_name ' ' names{c_index} ' ' num2str(i/nb_runs*100) '%'])
    end
    
end

%% CONVERGENCE

nb_fun = numel(fun);
for map_index = 1 : nb_maps
    
    chaotic_map_name = names{map_index};
    
    for i = 1 : nb_fun
        chaotic_fht.(chaotic_map_name).(function_name) = nan(nb_runs, 1);
        fht.(function_name) = nan(nb_runs, 1);
    end
    
    for i = 1 : nb_runs
        chaotic_fht.(chaotic_map_name).(function_name)(i) = chaotic_results.(chaotic_map_name).(function_name){i}.first_hitting_time;
        fht.(function_name)(i) = results.(function_name){i};
    end
    
    converged_percentage = sum(~isnan(chaotic_fht.(chaotic_map_name).(function_name)))/nb_runs*100;
    disp([' CHAOTIC ALGORITHM - Function name : ' upper(function_name) '     Convergence rate : ' num2str(converged_percentage) '%'])
    converged_percentage = sum(~isnan(fht.(function_name)))/nb_runs*100;
    disp([' STANDARD ALGORITHM - Function name : ' upper(function_name) '     Convergence rate : ' num2str(converged_percentage) '%'])
    
    % Remove NaNs
    temp = chaotic_fht.(chaotic_map_name).(function_name);
    chaotic_fht.(chaotic_map_name).(function_name) = temp(~isnan(temp));
    temp = fht.(function_name);
    fht.(function_name) = temp(~isnan(temp));
    
end

% NHST for equal means (2-sample t-test)

map = zeros(1,nb_maps);
pso = zeros(1,nb_maps);
cpso = zeros(1,nb_maps);
pval = zeros(1,nb_maps);
r = zeros(1,nb_maps);

for i = 1 : nb_maps
    
    chaotic_map_name = names{i};
    m_original = median(fht.(function_name));
    m_chaotic = median(chaotic_fht.(chaotic_map_name).(function_name));
    [p, ~, stats] = ranksum(fht.(function_name), chaotic_fht.(chaotic_map_name).(function_name));
    
    pso(i) = m_original;
    cpso(i) = m_chaotic;
    pval(i) = p;
    r(i) = stats.zval/sqrt(nb_runs);
   
end

% Confidence intervals for median differences

difference = zeros(1, nb_maps);
difference_min = zeros(1, nb_maps);
difference_max = zeros(1, nb_maps);

ci = struct();
nb_resamples = 10000;
for i = 1 : nb_maps
    temp = zeros(nb_resamples, 1);
    map_name = names{i};
    difference(i) = median(chaotic_fht.(map_name).(function_name)) - median(fht.(function_name));
    for j = 1 : nb_resamples
        sample1 = randsample(fht.(function_name), nb_runs, true);
        sample2 = randsample(chaotic_fht.(map_name).(function_name), nb_runs, true);
        temp(j) = median(sample2) - median(sample1);
    end
    difference_min(i) = prctile(temp, 5);
    difference_max(i) = prctile(temp, 95);
end

data = [pso; cpso; pval; r; difference; difference_min; difference_max];
results = table(data(:,1), data(:,2), data(:,3), data(:,4), data(:,5), data(:,6), 'VariableNames', names, 'RowNames', {'PSO', 'CPSO', 'p', 'r', 'diff', 'diff_min', 'diff_max'});

disp(results)

%% POSTPROCESS

h = zeros(nb_maps + 1,1);

% Plot chaotic CDF
for i = 1 : nb_maps
    [f1,x1] = ksdensity(chaotic_fht.(names{i}).(function_name),'support', 'positive','Function','cdf','NumPoints',5000);
    hold on
    h(i+1) = plot(x1, f1, 'linewidth', 1.15);
    set(gca, 'color', [253,245,230]/255)
    m1 = median(chaotic_fht.(names{i}).(function_name));
end

% Plot standard CDF
[f2,x2] = ksdensity(fht.(function_name),'support', 'positive','Function','cdf','NumPoints',5000);
h(1) = plot(x2, f2, 'linewidth', 3, 'color', 'k');
m2 = median(fht.(function_name));

% Axes
xlabel('First hitting time [nb. evaluations]')
set(gca, 'fontname', 'times', 'fontsize', 14)
xlim([0,15e3])

% Other
grid on 
box on

% Kolmogorov-Smirnov Test
ks_p_values = zeros(nb_maps, 1);
for i = 1 : nb_maps
    [~, ks_p_values(i)] = kstest2(fht.(function_name), chaotic_fht.(names{i}).(function_name));
end

% Legend
names = {'Chebyshev', 'Circle', 'Gauss', 'Logistic', 'Sine', 'Tent'};
lgnd = legend(h,'Original','Chebyshev', 'Circle', 'Gauss', 'Logistic', 'Sine', 'Tent', 'location', 'southeast');
set(lgnd, 'color', [1 1 1])


for i = 1 : nb_maps
    disp(median(chaotic_fht.(names{i}).(function_name)))
end
disp(median(fht.(function_name)))


%% CONTOUR PLOT WITH CENTROIDS

index = 12;

contour_plot_results = chaotic_results.Circle.rosenbrock{index};

nb_epochs = contour_plot_results.nb_epochs;

% Contour plot
resolution = 1000;
xx = linspace(-5.12, 5.12, resolution);
yy = linspace(-5.12, 5.12, resolution);
zz = zeros(resolution);
for i = 1 : resolution
    for j = 1 : resolution
        zz(i,j) = fun([xx(i), yy(j)]);
    end
end
contourf(xx, yy, zz');
set(gca, 'fontname', 'times', 'fontsize', 14)
colormap(jet)
colorbar
hold on

% Centroids
cx = nan(nb_epochs + 1, 1);
cy = nan(nb_epochs + 1, 1);
for i = 1 : nb_epochs + 1
    positions = contour_plot_results.positions_history{i};
    cx(i) = mean(positions(:,1));
    cy(i) = mean(positions(:,2));
end

% Plot centroids
for i = 1 : nb_epochs + 1
    plot(cx, cy, '-+r')
    hold on
end

% Plot global minima
plot(1, 1, 'yp', 'markersize', 20, 'markerfacecolor', 'y', 'markeredgecolor', 'k')

%% CENTROID DISTANCE FROM GLOBAL MIN

centroid_dist_hist = sum(sqrt(([ccx, ccy] - repmat([1,1], size(ccx,1), 1)).^2), 2);
plot(centroid_dist_hist);

%% PLOT CHAOTIC MAPS

n = 50;
for i = 1 : 6
    subplot(3, 2, i)
    plot(chaos(cindex(i), n), 'k', 'linewidth', 0.75)
    hold on
    xlim([0, n])
    title(names{i}, 'fontweight', 'normal')
    set(gca, 'fontname', 'times', 'fontsize', 16, 'xtick', [], 'ytick', [])
end

disp(nanmedian(pso_fht.(function_name)))
for i = 1 : nb_maps
    disp(nanmedian(chaotic_pso_fht.(names{i}).(function_name))) 
end

%% CHAOTIC MAP EXAMPLES

n = 40;
time = 0:n-1;
x = zeros(n, 1);
y = zeros(n, 1);
random_number = 0.4;
x(1) = random_number;
y(1) = x(1) + 1e-6;
for i=2:n
    x(i)=mod(4*x(i-1)*(1-x(i-1)),1);
    y(i)=mod(4*y(i-1)*(1-y(i-1)),1);
end
h1 = plot(time, y, 'r', 'linewidth', 1);
hold on
h2 = plot(time, x, 'k', 'linewidth', 1);
xlim([0, n + 1])
ylim([-1, 2])
xlabel('Time step', 'fontname', 'times', 'fontsize', 14)
ylabel(' ', 'fontname', 'times', 'fontsize', 14)
set(gca, 'fontname', 'times', 'fontsize', 14)
lgnd = legend([h1,h2],'Original + 1e-8','Original');
set(lgnd, 'color', [1 1 1], 'fontsize', 16)
grid

h3 = plot(abs(x-y), 'k', 'linewidth', 1);
ylabel('Difference', 'fontname', 'times', 'fontsize', 14)
xlabel('Time step', 'fontname', 'times', 'fontsize', 14)
set(gca, 'fontname', 'times', 'fontsize', 14)
grid

%% EMPIRICAL PMF

histogram(chaos(1,100000), 'normalization', 'probability')
set(gca, 'fontname', 'times', 'fontsize', 14, 'color', [253,245,230]/255)

%% DIVERSITY

index = 20;
h1 = plot(chaotic_results.Chebyshev.(function_name){index}.diversity, 'k');
hold on
h2 = plot(results.(function_name){index}.diversity);
legend([h1, h2], 'Chaotic', 'Standard')

temp1 = zeros(nb_runs, 1);
temp2 = zeros(nb_runs, 1);
for i  = 1 : nb_runs
   temp1(i) = median(chaotic_results.Circle.(function_name){i}.diversity);
   temp2(i) = median(results.(function_name){i}.diversity);
end
disp(median(temp1))
disp(median(temp2))

%% ANIMATE

index = 15;
animate_results = results.(function_name){index}.diversity;
nb_results = numel(results);
for i = 30 : 34
    cla
    h = plot(results{i}(:,1), results{i}(:,2), 'ko', 'markerfacecolor', 'k');
    pause(0.25)
    drawnow
    hold on
end