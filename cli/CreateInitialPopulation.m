function initial_population = CreateInitialPopulation(PopSize)
% CREATEINITIALPOPULATION This function creates a 3xPopSize matrix,
% which is PopSize different P-I-D gains.
    initial_population = 100 * rand([3 PopSize]);
end