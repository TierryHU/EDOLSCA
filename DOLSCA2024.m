function [ge, pop] = DOLSCA(pop, Gm_o, D, Np, lb, ub, fobj)
    % Dynamic Opposite Learning Sine Cosine Algorithm (DOLSCA) for optimization
    % Inputs:
    %   pop - Initial population matrix
    %   Gm_o - Maximum number of iterations
    %   D - Dimensionality of the problem
    %   Np - Number of particles in the population
    %   lb - Lower bound for each dimension
    %   ub - Upper bound for each dimension
    %   fobj - Objective function handle

    disp('DOLSCA Algorithm Started');
    Jr = 0.5;        % Probability threshold for dynamic opposite phase
    w = 12;          % Weight for opposite learning term
    Gm = Gm_o;       % Total generations
    G = 1;           % Initialize generation counter
    
    % Initialize destination (best) position and fitness
    Destination_position = zeros(1, D);
    Destination_fitness = inf;
    Convergence_curve = zeros(1, Gm);
    Objective_values = zeros(1, Np);

    % Evaluate initial population fitness
    for i = 1:Np
        Objective_values(i) = fobj(pop(i, :), D);
        if Objective_values(i) < Destination_fitness
            Destination_position = pop(i, :);
            Destination_fitness = Objective_values(i);
        end
    end
    
    % Sine-Cosine Algorithm Main Loop
    while G <= Gm
        a = 2;
        r1 = a - G * (a / Gm); % Decrease r1 linearly with iterations

        % Generate new solutions using Sine and Cosine transformations
        for i = 1:Np
            for j = 1:D
                r2 = 2 * pi * rand();
                r3 = 2 * rand();
                r4 = rand();

                if r4 < 0.5
                    pop_new(i, j) = pop(i, j) + r1 * sin(r2) * abs(r3 * Destination_position(j) - pop(i, j));
                else
                    pop_new(i, j) = pop(i, j) + r1 * cos(r2) * abs(r3 * Destination_position(j) - pop(i, j));
                end
            end
        end

        % Boundary correction
        for i = 1:Np
            for j = 1:D
                if pop_new(i, j) < lb(j)
                    pop_new(i, j) = lb(j) + rand * (ub(j) - lb(j));
                elseif pop_new(i, j) > ub(j)
                    pop_new(i, j) = ub(j) - rand * (ub(j) - lb(j));
                end
            end
        end

        % Update particles based on fitness
        for i = 1:Np
            fitnew = fobj(pop_new(i, :), D);
            if fitnew < Objective_values(i)
                pop(i, :) = pop_new(i, :);
                Objective_values(i) = fitnew;
                if fitnew < Destination_fitness
                    Destination_fitness = fitnew;
                    Destination_position = pop_new(i, :);
                end
            end
        end

        % Record the best fitness at each generation
        ge(G) = Destination_fitness;
        G = G + 1;

        % Dynamic Opposite Phase
        if rand < Jr
            for i = 1:Np
                for j = 1:D
                    Upperbound1 = max(pop(:, j));
                    Lowerbound1 = min(pop(:, j));
                    op(i, j) = Upperbound1 + Lowerbound1 - pop(i, j);
                end
                popO(i, :) = pop(i, :) + w * rand * (rand * op(i, :) - pop(i, :));
                popO_new(i, :) = pop(i, :);
                popO_new(Np + i, :) = popO(i, :);
            end

            % Boundary correction for opposite solutions
            for i = 1:Np
                for j = 1:D
                    if popO_new(i, j) < lb(j)
                        popO_new(i, j) = lb(j) + rand * (ub(j) - lb(j));
                    elseif popO_new(i, j) > ub(j)
                        popO_new(i, j) = ub(j) - rand * (ub(j) - lb(j));
                    end
                end
            end

            % Evaluate and update opposite solutions
            for i = 1:Np
                fitnew = fobj(popO_new(i, :), D);
                if fitnew < Objective_values(i)
                    pop(i, :) = popO_new(i, :);
                    Objective_values(i) = fitnew;
                    if fitnew < Destination_fitness
                        Destination_fitness = fitnew;
                        Destination_position = popO_new(i, :);
                    end
                end
            end

            % Update the convergence curve for the current generation
            if G == 1
                ge(G) = Destination_fitness;
            else
                if Destination_fitness < ge(G - 1)
                    ge(G) = Destination_fitness;
                else
                    ge(G) = ge(G - 1);
                end
            end
            G = G + 1;
        end
    end

    % Clean up the convergence curve if excess iterations are present
    if size(ge, 2) == 30001
        ge(:, 1) = [];
    elseif size(ge, 2) == 30002
        ge(:, 1:2) = [];
    end
end
