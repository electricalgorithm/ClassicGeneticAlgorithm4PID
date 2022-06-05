function new_chromosome = CGA_Crossover(choromosome_ma, choromosome_pa, ...
                                        crossover_chance)
    
    % Check if the chromosomes have same gen amount.
    if (size(choromosome_ma, 2) ~= size(choromosome_pa, 2))
        error("Chromosomes must be in same length.");
    end
    
    % If crossover chance reached, do the process.
    active_chance = randi(100) / 100;
    if active_chance >= crossover_chance
        % Pre-allocate the array.
        new_chromosome = repmat('0', 1, size(choromosome_pa, 2));
    
        % Decide the slice edge for chromosome_ma.
        edge_percentage = randi(100) / 100;
        selection_edge = ceil(size(choromosome_ma, 2) * edge_percentage);
    
        % Assign the values.
        new_chromosome(1, 1:selection_edge) = choromosome_ma(1, 1:selection_edge);
        new_chromosome(1, (selection_edge+1):end) = ...
                  choromosome_pa(1, (selection_edge+1):size(choromosome_pa, 2));
    
        % If crossover chance couldn't achieved, assign the mother chromosome
    % into the child.
    else
        new_chromosome = choromosome_ma;
    end
end