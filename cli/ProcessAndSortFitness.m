function SortedFitness = ProcessAndSortFitness(PID_Population, PlantObject, RefValue)
% PROCESSANDSORTFITNESS The function returns a matrix, which is sorted
% and has four row. The first three row is P-I-D gains, and the last
% one is the FitnessFunction value.
    
    [~, population_size] = size(PID_Population);
    SortedFitness = zeros([4, population_size]);
    
    % Assign fourth row as Fitness Values.
    for index = 1:1:population_size
        FitnesValue = FitnessFunction(PlantObject, RefValue, PID_Population(:, index));
        SortedFitness(:, index) = [PID_Population(:, index); FitnesValue];
    end

    % Sort the PID Gains according to their Fitness Values.
    SortedFitness = sortrows(SortedFitness', 4, "descend")';
end