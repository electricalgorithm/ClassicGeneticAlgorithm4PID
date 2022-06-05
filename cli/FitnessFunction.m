function fitness_value = FitnessFunction(PlantObject, RefValue, PID_coeff)
    [amps, timescale] = SimulatePIDSystem(RefValue, PID_coeff, PlantObject);
    [overshoot, stability_time, err_given_t] = AnalyseSystemResult(timescale, amps, RefValue);
   
    % Definitations for Fitness Function
    weight_t = 0.5;
    weight_o = 0.4;
    weight_e = 0.1;

    % Get the last time of time scale.
    stability_time_ratio = stability_time / timescale(end);
    
    % Create the fitness value. Bigger value means the most fittest.
    fitness_value = -(weight_t * stability_time_ratio + ...
                      weight_o * overshoot            + ...
                      weight_e * err_given_t            ...
                      );
                
end