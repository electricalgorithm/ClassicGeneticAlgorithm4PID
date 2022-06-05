function [group_1, group_2] = SelectionMethod(Individiuals, PairCount)
% SELECTIONBYELITMETHOD This function returns the pairs of individuals
% selected with Elit Method of Selection in Classic Genetic Algorithm.
    
    % Get the sizes of individuals.
    IndividiualsSize = size(Individiuals, 2);
    
    % Select the groups which will be fathers and mothers.
    group_1 = [Individiuals(:, randi(IndividiualsSize, ...
                                  [1, PairCount]))];
    group_2 = [Individiuals(:, randi(IndividiualsSize, ...
                                  [1, PairCount]))];

    % Shuffle the groups.
    group_1 = group_1(:, randi(size(group_1, 2), [1, PairCount]));
    group_2 = group_2(:, randi(size(group_2, 2), [1, PairCount]));
end