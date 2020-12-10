function output = simulated_annealing(input)

% --- Input checks

% CHECK input field names
expected_fieldnames = {'fun';                 ...
                       'nb_dim';              ...
                       'initial_position';   ...
                       'lb';                  ...
                       'ub';                  ...
                       'max_iter';            ...
                       'max_temp';            ...
                       'known_best_fitness';  ...
                       'tol';                 ...
                       'positions_hist_flag'; ...
                       'diversity_hist_flag'};
                  
for i = 1 : numel(expected_fieldnames)
    input_fieldnames = fieldnames(input);
    if ~all(strcmpi(expected_fieldnames, input_fieldnames))
        error('SAAlgorithm:InputStruc', 'Missing fields in input structure.')
    end
    evalc([input_fieldnames{i} ' = input.' input_fieldnames{i}]);
end

% CHECK objective function
if ~isa(fun, 'function_handle')
    error('SimulatedAnnealing:ObjectiveFunction','Invalid objective function.')
end

% CHECK input classes and attributes
inputs = {nb_dim, lb, ub, max_iter, max_temp, tol};
for i = 1 : 5
   validateattributes(inputs{i}, {'numeric'}, {'finite', 'nonempty', 'nonnan'}) 
end

% CHECK number of dimensions
if nb_dim <= 0
    error('SimulatedAnnealing:Dimensions','Invalid number of dimensions.')
end

% CHECK boundaries
if ~(lb < ub)
   error('SimulatedAnnealing:Boundaries','Invalid boundaries.') 
end

% Initial position
if isempty(initial_position)
    
    % Random initial position
    position = lb * ones(1,nb_dim) + rand(1,nb_dim) * (ub-lb);

else
    
    % CHECK if dimensions match
    if len(initial_position) ~= nb_dim
        error('SimulatedAnnealing:NbDimensions', 'Initial position and number of dimensions do not match.')
    end
    
    % Input initial position
    position = initial_position;
    
end

% Initialize positions history
if positions_hist_flag
    positions_history = position;
else
    positions_history = [];
end

% CHECK max_temp
if max_temp <= 0 
    error('SimulatedAnnealing:MaximumTemperature','Invalid maximum temperature.') 
end

% --- MAIN

% Initial best position
best_position = position;

% Initial fitness value
fitness = fun(position);

% Initial number of objective function evaluations
nb_fun_eval = 1;

% Initial best fitness value
best_fitness = fitness;

% Initial first hitting time
first_hitting_time = nan;

% Calculate delta temperature for linearly decreasing schedule
delta_temp = max_temp / (max_iter - 1);

% Initial acceptance history
acceptance_history = [];

% Main loop
for iteration = 1 : max_iter 
    
    % Calculate temperature
    temperature = max_temp - (iteration-1) * delta_temp;

    % New position
    rand_integer = randi(2)-1;
    if rand_integer == 0
        step = ub - position;
    else
        step = position - lb;
    end
    new_position = position + (-1)^rand_integer .* rand(1, nb_dim) .* step;
    
    % New fitness
    new_fitness = fun(new_position);
    
    % Increment number of objetive function evaluations
    nb_fun_eval = nb_fun_eval + 1;
    
    p = rand() < exp(-(new_fitness-fitness)/temperature);
    if new_fitness <= fitness
        acceptance_history = [acceptance_history; iteration,1];
    else
        if p == 0
            p = -1;
        end
        acceptance_history = [acceptance_history; iteration,p];
    end
    
    if new_fitness < fitness || p
        
        % Update position and fitness
        position = new_position;
        fitness = new_fitness;
        
        % Update best position and fitness
        if new_fitness < best_fitness
            best_position = position;
            best_fitness = fitness;
        end
        
        % Check for convergence
        if abs(best_fitness - known_best_fitness) <= tol
            first_hitting_time = nb_fun_eval;
            break
        end
        
    end
    
    % Update positions history
    if positions_hist_flag
        positions_history = [positions_history; position];
    end
    
    % Check for convergence
    if abs(best_fitness - known_best_fitness) <= tol
        first_hitting_time = nb_fun_eval;
        break
    end
    
end

%% OUTPUT

output.best_fitness = best_fitness;
output.best_position = best_position;
output.first_hitting_time = first_hitting_time;
output.acceptance_history = acceptance_history;
output.positions_history = positions_history;

end