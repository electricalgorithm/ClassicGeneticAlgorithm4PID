classdef ClassicGeneticAlgorithmForPID_asMFile < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        ClassicGeneticAlgorithmforPIDTunningUIFigure  matlab.ui.Figure
        GridLayout                    matlab.ui.container.GridLayout
        ImportaSimulationButton       matlab.ui.control.Button
        StoptheSimulationButton       matlab.ui.control.Button
        SimulationResultsPanel        matlab.ui.container.Panel
        SavetheSimulationButton       matlab.ui.control.Button
        SaveResultsasTXTButton        matlab.ui.control.Button
        KdTextSR                      matlab.ui.control.Label
        kd_output                     matlab.ui.control.NumericEditField
        KiTextSR                      matlab.ui.control.Label
        ki_output                     matlab.ui.control.NumericEditField
        fitnessvalue_output           matlab.ui.control.NumericEditField
        FitnessValueTextSR            matlab.ui.control.Label
        KpTextSR                      matlab.ui.control.Label
        kp_output                     matlab.ui.control.NumericEditField
        generation_output             matlab.ui.control.NumericEditField
        GenerationTextSR              matlab.ui.control.Label
        BestGuessTextSR               matlab.ui.control.Label
        SimulationStepsTextSR         matlab.ui.control.Label
        OutputGraph                   matlab.ui.control.UIAxes
        StarttheSimulationButton      matlab.ui.control.Button
        TabGroup                      matlab.ui.container.TabGroup
        PlantModelSettingsTab         matlab.ui.container.Tab
        SaveButton                    matlab.ui.control.Button
        SavedStatusOutput             matlab.ui.control.Label
        DescriptionTextPlantModel_2   matlab.ui.control.Label
        SPlusOneTextPlantModel        matlab.ui.control.Label
        firstorder_input              matlab.ui.control.NumericEditField
        sTextPlantModel               matlab.ui.control.Label
        EulerTextPlantModel           matlab.ui.control.Label
        NegativeTextPlantModel        matlab.ui.control.Label
        delay_input                   matlab.ui.control.NumericEditField
        DotTextPlantModel             matlab.ui.control.Label
        gain_input                    matlab.ui.control.NumericEditField
        UnderlineFractionTextPlantModel  matlab.ui.control.Label
        DescriptionTextPlantModel     matlab.ui.control.Label
        PlantModelText                matlab.ui.control.Label
        UnderlineTextPlantModel       matlab.ui.control.Label
        AlgorithmSettingsTab          matlab.ui.container.Tab
        inversionchance_input         matlab.ui.control.NumericEditField
        mutationchance_input          matlab.ui.control.NumericEditField
        crossoverchance_input         matlab.ui.control.NumericEditField
        InversionTextCGAS             matlab.ui.control.Label
        MutationTextCGAS              matlab.ui.control.Label
        CrossoverTextCGAS             matlab.ui.control.Label
        RequiredChancesText           matlab.ui.control.Label
        MaxGenerationTextCGAS         matlab.ui.control.Label
        FitnessFunctionTextCGAS       matlab.ui.control.Label
        maxgeneration_input           matlab.ui.control.NumericEditField
        fitnessbigger_input           matlab.ui.control.NumericEditField
        EndingConditionText           matlab.ui.control.Label
        referencevalue_input          matlab.ui.control.NumericEditField
        ReferenceValueEditFieldLabel  matlab.ui.control.Label
        populationsize_input          matlab.ui.control.NumericEditField
        PopulationSizeEditFieldLabel  matlab.ui.control.Label
        EnvironmentSettingsText       matlab.ui.control.Label
        UnderlineTextCGAS             matlab.ui.control.Label
        ClassicGeneticAlgorithmSettingsLabel  matlab.ui.control.Label
    end

    
    properties (Access = private)
        % Application Info
        VersionInfo = "0.2 Beta";
        GitHubLink = "https://github.com/electricalgorithm/ClassicGeneticAlgorithm4PID";
        
        % Algorithm Details
        LIMIT_NUMBER = 128;
        PlantModel          % Which stores the TF object of the Plant Model.
        PlantModelInputs = [0, 0, 1]
        Chances = [0.5, 0.5, 0.5]; % Stores the chances for Mut, CO and Inv.
        MaxGeneration = 20; % Stores the maximum generation of the same best.
        ReferenceValue = 1; % Stores the reference value for system.
        FitnessValueToStop = -0.4; % Stores the ending condition of FF.
        PopulationSize = 20; % Stores the population size genetered.
        
        % Execution
        StopFlag = false;
        
        % Results
        BestPIDResult;
        BestFitnessValue;
        GenerationCreated;
        LastPopulation;
        IsImporting = false;
    end
    
    methods (Access = private)
        
        function ClassicGeneticAlgorithm(app)
            % Create the initial population, and start the algorithm.
            population = app.CreateInitialPopulation();
            sorted_population = app.ProcessAndSortFitness(population);
            
            % If the generation will be same for EndingCondGenCount 
            % generation, exit the program.
            being_same_generation = 0;
            previous_best = sorted_population(:, 1);
            
            % Count the generation for investigation purposes.
            app.GenerationCreated = 1;
        
            % If the ending condition is not met, continue the process.
            while sorted_population(4, 1) <= app.FitnessValueToStop && ...
                    being_same_generation <= app.MaxGeneration

                % Exit the loop.
                if (app.StopFlag == true)
                    break;
                end
                
                pause(0.0001);
                app.generation_output.Value = app.GenerationCreated;
        
                % Select best two individuals, put them into the population.
                best_two = [sorted_population(:, 1), sorted_population(:, 2)];
                population(:, [1, 2]) = best_two(1:3, :);
        
                % Select individuals for crossover, mutation and inversion.
                [group_mas, group_pas] = app.SelectionMethod(sorted_population);
        
                % Travel throught to apply mutation, inversion, and crossover.
                for index = 1:size(group_pas, 2)
                    % Convert PID values into gens.
                    MA_PID_chromosomes = app.PIDtoBinaryGens(group_mas(1:3, index));
                    PA_PID_chromosomes = app.PIDtoBinaryGens(group_pas(1:3, index));
        
                    % Apply Crossover
                    crossovered_chromosome(1, :) = ...
                        app.CGA_Crossover(MA_PID_chromosomes(1, :),...
                                          PA_PID_chromosomes(1, :));
                    crossovered_chromosome(2, :) = ...
                        app.CGA_Crossover(MA_PID_chromosomes(2, :),...
                                          PA_PID_chromosomes(2, :));
                    crossovered_chromosome(3, :) = ...
                        app.CGA_Crossover(MA_PID_chromosomes(3, :),...
                                          PA_PID_chromosomes(3, :));
                    
                    % Apply Mutation
                    mutated_chromosome(1, :) = ...
                        app.CGA_Mutation(crossovered_chromosome(1, :));
                    mutated_chromosome(2, :) = ...
                        app.CGA_Mutation(crossovered_chromosome(2, :));
                    mutated_chromosome(3, :) = ...
                        app.CGA_Mutation(crossovered_chromosome(3, :));
        
                    % Apply Inversion
                    inverted_chromosome(1, :) = ...
                        app.CGA_Inversion(mutated_chromosome(1, :));
                    inverted_chromosome(2, :) = ...
                        app.CGA_Inversion(mutated_chromosome(2, :));
                    inverted_chromosome(3, :) = ...
                        app.CGA_Inversion(mutated_chromosome(3, :));
        
                    % Convert binary gens into decimal PID format.
                    PID_format = app.BinaryGenstoPID(inverted_chromosome);
                    
                    % Check for PID if 0 or LIMIT_NUMBER.
                    if PID_format(1) == 0, PID_format(1) = 1; end
                    if PID_format(1) == app.LIMIT_NUMBER, PID_format(1) = app.LIMIT_NUMBER - 1; end
                    if PID_format(2) == 0, PID_format(2) = 1; end
                    if PID_format(2) == app.LIMIT_NUMBER, PID_format(2) = app.LIMIT_NUMBER - 1; end
                    if PID_format(3) == 0, PID_format(3) = 1; end
                    if PID_format(3) == app.LIMIT_NUMBER, PID_format(3) = app.LIMIT_NUMBER - 1; end
        
                    % Assert into the population.
                    population(1:3, index + 2) = PID_format;
                end
        
                % New generation created after previous for loop.
                app.GenerationCreated = app.GenerationCreated + 1;
        
                % Put them into fitness function, and sort the new generation.
                sorted_population = app.ProcessAndSortFitness(population);
                
                % If the previous_best is the same as this best, increase
                % the number of being_same_generation.
                if (previous_best(:) == sorted_population(:, 1))
                    being_same_generation = being_same_generation + 1;
                else
                    % Update the previous best.
                    previous_best(:) = sorted_population(:, 1);
                end
            end
        
            % If the ending condition met, stop the process, and return.
            app.BestPIDResult = sorted_population([1, 2, 3], 1);
            app.BestFitnessValue = sorted_population(4, 1);
            app.LastPopulation = population;
        
        end
        
        function [overshoot_error, stability_time, error_given_time] = AnalyseSystemResult(~, res_t, res_y, ref_val, error_at_time)

            % Get the size of timing.
            [size_of_time, ~] = size(res_t);
        
            % Check for the parameters if they are given.
            if ~exist("error_at_time", "var")
                error_at_time = res_t(size_of_time - 5);
            end
            
            % Find the overshoot:
            %   Overshoot is the maximum distance that our function
            %   rises comparing to our reference.
            overshoot = max(res_y);
            if (overshoot < ref_val)
                overshoot_error = 0;
            else
                overshoot_error = (overshoot - ref_val) / ref_val;
            end
        
            % Find the necessary time to achieve stability,
            % and if it does not until the max time, find
            % the error rate at the ending time.
            % According to needs, change the stability criteria.
            stability_criteria = 0.05;
            stability_time = 0;
        
            for index = 1 : (size_of_time - 5)
                checking_values = res_y(index : (index + 5));
                average_values = mean(checking_values);
        
                % Check if the average inside the stability criteria.
                current_error = abs(average_values - ref_val) / ref_val;
                if current_error <= stability_criteria
                    stability_time = res_t(index);
                end
        
                % Check if it is in error_at_time state.
                if res_t(index) == error_at_time
                    error_given_time = current_error;
                end
        
                % If loop reaches the end of the time, stability
                % doesn't exist, and return -1.
                if index == (size_of_time - 5) && stability_time == 0
                    stability_time = -1;
                end
            end
        end    
        
        function PID_Coefficents = BinaryGenstoPID(~, CombinedGen)
            [~, NucleoditSize] = size(CombinedGen);
        
           % Check if nucleodit size is even number.
            if (mod(NucleoditSize, 2) ~= 0)
                error("Nucleodit size has to be even number.")
            end
        
            PID_Coefficents = [
                bin2dec(CombinedGen(1, 1:(NucleoditSize/2))) + 0.01 * bin2dec(CombinedGen(1, (NucleoditSize/2 + 1):end));
                bin2dec(CombinedGen(2, 1:(NucleoditSize/2))) + 0.01 * bin2dec(CombinedGen(2, (NucleoditSize/2 + 1):end));
                bin2dec(CombinedGen(3, 1:(NucleoditSize/2))) + 0.01 * bin2dec(CombinedGen(3, (NucleoditSize/2 + 1):end))
           ];
        end
        
        function new_chromosome = CGA_Crossover(app, choromosome_ma, choromosome_pa)
    
            % Check if the chromosomes have same gen amount.
            if (size(choromosome_ma, 2) ~= size(choromosome_pa, 2))
                error("Chromosomes must be in same length.");
            end
            
            % If crossover chance reached, do the process.
            active_chance = randi(100) / 100;
            if active_chance >= app.Chances(1)
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
        
        function chromosome = CGA_Inversion(app, chromosome)
            % If the chance achieved, then do the inversion.
            active_chance = randi(100) / 100;
            if active_chance >= app.Chances(3)
        
                % Init two random value, and assign them into start and end indices.
                rand_1 = randi(size(chromosome, 2));
                rand_2 = randi(size(chromosome, 2));
                
                if (rand_1 < rand_2)
                    starting_index = rand_1;
                    ending_index = rand_2;
                else
                    starting_index = rand_2;
                    ending_index = rand_1;
                end
            
                % Flip the bits.
                for gen_index = starting_index:ending_index
                    if (chromosome(1, gen_index) == '0')
                        chromosome(1, gen_index) = '1';
                    else
                        chromosome(1, gen_index) = '0';
                    end
                end
        
        
            end
        end
        
        function mutated_chromosome = CGA_Mutation(app, original_chromosome)
            mutated_chromosome = original_chromosome;
        
            for index = 1:size(original_chromosome, 2)
                random_chance = randi(100) / 100;
        
                % If the randomized chance bigger then chance of change, flip it.
                if random_chance >= app.Chances(3)
                    char_chance = randi(100);
                    
                    % Change char_chance into the real gen value.
                    if char_chance >= 50
                        char_chance = '1';
                    else
                        char_chance = '0';
                    end
        
                    mutated_chromosome(1, index) = char_chance;
                end
        
            end
        end
        
        function initial_population = CreateInitialPopulation(app)
            initial_population = 100 * rand([3 app.PopulationSize]);
        end

        function plant = CreatePlantObject(~, k, T, tau)
            plant = tf(k, [T 1], 'InputDelay', tau);
        end
        
        function fitness_value = FitnessFunction(app, PID_coeff)
            [amps, timescale] = app.SimulatePIDSystem(app.ReferenceValue, PID_coeff, app.PlantModel);
            [overshoot, stability_time, err_given_t] = app.AnalyseSystemResult(timescale, amps, app.ReferenceValue);
           
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
        
        function [EachCombined, BinaryGens] = PIDtoBinaryGens(app, PIDCoefficents)
            % Check for argument's size.
            if size(PIDCoefficents) ~= [3, 1]
                error("The argument provided is incorrect size.");
            end
        
            Kp = PIDCoefficents(1);
            Ki = PIDCoefficents(2);
            Kd = PIDCoefficents(3);
        
            % Save 13.1257 as [13, 12], or 53.987 as [53, 98]
            Kp_arr = [floor(Kp), floor((Kp - floor(Kp))*100)];
            Ki_arr = [floor(Ki), floor((Ki - floor(Ki))*100)];
            Kd_arr = [floor(Kd), floor((Kd - floor(Kd))*100)];
        
            % Check the numbers if they are bigger than limit number.
            if Kp_arr(1) >= app.LIMIT_NUMBER || ...
               Kp_arr(2) >= app.LIMIT_NUMBER || ...
               Ki_arr(1) >= app.LIMIT_NUMBER || ...
               Ki_arr(2) >= app.LIMIT_NUMBER || ...
               Kd_arr(1) >= app.LIMIT_NUMBER || ...
               Kd_arr(2) >= app.LIMIT_NUMBER
                warning("You shouldn't put a number which has more than 7 digits" + ...
                      " in binary format, or put a number that is bigger than 127.")
        
            end
        
            % Correct each setting if they are 0 or 128.
            Kp_arr(Kp_arr == 0) = Kp_arr(Kp_arr == 0) + 1;
            Kp_arr(Kp_arr == app.LIMIT_NUMBER) = Kp_arr(Kp_arr == app.LIMIT_NUMBER) - 1;
            Ki_arr(Ki_arr == 0) = Ki_arr(Ki_arr == 0) + 1;
            Ki_arr(Ki_arr == app.LIMIT_NUMBER) = Ki_arr(Ki_arr == app.LIMIT_NUMBER) - 1;
            Kd_arr(Kd_arr == 0) = Kd_arr(Kd_arr == 0) + 1;
            Kd_arr(Kd_arr == app.LIMIT_NUMBER) = Kd_arr(Kd_arr == app.LIMIT_NUMBER) - 1;
        
            % Convert each decimal into binary.
            Kp_bin = dec2bin(Kp_arr, floor(log2(app.LIMIT_NUMBER)));
            Ki_bin = dec2bin(Ki_arr, floor(log2(app.LIMIT_NUMBER)));
            Kd_bin = dec2bin(Kd_arr, floor(log2(app.LIMIT_NUMBER)));
        
            % Return BinaryGens.
            % Format: [Kp_int, Kp_float; Ki_int, Ki_float, Kd_int, Kd_float];
            BinaryGens = [Kp_bin; Ki_bin; Kd_bin];
            
            % Return EachCombined
            Kp_combined(1:7) = Kp_bin(1, :);
            Kp_combined(8:14) = Kp_bin(2, :);
            Ki_combined(1:7) = Ki_bin(1, :);
            Ki_combined(8:14) = Ki_bin(2, :);
            Kd_combined(1:7) = Kd_bin(1, :);
            Kd_combined(8:14) = Kd_bin(2, :);
            EachCombined = [Kp_combined; Ki_combined; Kd_combined];
        end
        
        function SortedFitness = ProcessAndSortFitness(app, PID_Population)
            [~, population_size] = size(PID_Population);
            SortedFitness = zeros([4, population_size]);
            
            % Assign fourth row as Fitness Values.
            for index = 1:1:population_size
                FitnesValue = app.FitnessFunction(PID_Population(:, index));
                SortedFitness(:, index) = [PID_Population(:, index); FitnesValue];
            end
        
            % Sort the PID Gains according to their Fitness Values.
            SortedFitness = sortrows(SortedFitness', 4, "descend")';
            
            % Display in the GUI
            app.kp_output.Value = SortedFitness(1, 1);
            app.ki_output.Value = 1/SortedFitness(2, 1);
            app.kd_output.Value = 1/SortedFitness(3, 1); 
            app.fitnessvalue_output.Value = SortedFitness(4, 1);
        end
        
        function [group_1, group_2] = SelectionMethod(app, Individiuals)
            % Get the sizes of individuals.
            IndividiualsSize = size(Individiuals, 2);
            
            % Select the groups which will be fathers and mothers.
            group_1 = [Individiuals(:, randi(IndividiualsSize, ...
                                        [1, app.PopulationSize - 2]))];
            group_2 = [Individiuals(:, randi(IndividiualsSize, ...
                                          [1, app.PopulationSize - 2]))];
        
            % Shuffle the groups.
            group_1 = group_1(:, randi(size(group_1, 2), ...
                                [1, app.PopulationSize - 2]));
            group_2 = group_2(:, randi(size(group_2, 2), ...
                                [1, app.PopulationSize - 2]));
        end
        
        function [y_data, x_data] = SimulatePIDSystem(~, ref_val, PID_gains, plant_model)
            % Create the PID controller with given gains.
            pid_controller = pid(PID_gains(1), 1/PID_gains(2), 1/PID_gains(3));
            % Create a system to simulate.
            feedback_system = feedback(pid_controller * plant_model, 1);
            % Simulate the system with given reference step function.
            step_options = stepDataOptions('StepAmplitude', ref_val);
            [y_data, x_data] = step(feedback_system, step_options);
        end
        
        function SimulationEndedConditions(app)
            % Show the Results in the GUI
            app.generation_output.Value = app.GenerationCreated;
            app.kp_output.Value = app.BestPIDResult(1);
            app.ki_output.Value = 1/app.BestPIDResult(2);
            app.kd_output.Value = 1/app.BestPIDResult(3);
            app.fitnessvalue_output.Value = app.BestFitnessValue;

            % Print to figure.
            app.OutputGraph.Visible = "on";
            [y_data, x_data] = app.SimulatePIDSystem(app.ReferenceValue,...
                                                     app.BestPIDResult,...
                                                     app.PlantModel);
            plot(app.OutputGraph, x_data, y_data);

            % Make the Buttons Visible
            app.SaveResultsasTXTButton.Enable = "on";
            app.SavetheSimulationButton.Enable = "on";

            % Disable the Stop Simulation Button.
            app.StoptheSimulationButton.Enable = "off";
            app.StarttheSimulationButton.Enable = 'on';
        end
        
        function ContinueClassicGeneticAlgorithm(app)
            % Create the initial population, and start the algorithm.
            population = app.LastPopulation;
            sorted_population = app.ProcessAndSortFitness(population);
            
            % If the generation will be same for EndingCondGenCount 
            % generation, exit the program.
            being_same_generation = 0;
            previous_best = sorted_population(:, 1);
        
            % If the ending condition is not met, continue the process.
            while sorted_population(4, 1) <= app.FitnessValueToStop && ...
                    being_same_generation <= app.MaxGeneration

                pause(0.0001);
                
                if app.StopFlag == true
                    break;
                end
                
                app.generation_output.Value = app.GenerationCreated;
        
                % Select best two individuals, put them into the population.
                best_two = [sorted_population(:, 1), sorted_population(:, 2)];
                population(:, [1, 2]) = best_two(1:3, :);
        
                % Select individuals for crossover, mutation and inversion.
                [group_mas, group_pas] = app.SelectionMethod(sorted_population);
        
                % Travel throught to apply mutation, inversion, and crossover.
                for index = 1:size(group_pas, 2)
                    % Convert PID values into gens.
                    MA_PID_chromosomes = app.PIDtoBinaryGens(group_mas(1:3, index));
                    PA_PID_chromosomes = app.PIDtoBinaryGens(group_pas(1:3, index));
        
                    % Apply Crossover
                    crossovered_chromosome(1, :) = ...
                        app.CGA_Crossover(MA_PID_chromosomes(1, :),...
                                          PA_PID_chromosomes(1, :));
                    crossovered_chromosome(2, :) = ...
                        app.CGA_Crossover(MA_PID_chromosomes(2, :),...
                                          PA_PID_chromosomes(2, :));
                    crossovered_chromosome(3, :) = ...
                        app.CGA_Crossover(MA_PID_chromosomes(3, :),...
                                          PA_PID_chromosomes(3, :));
                    
                    % Apply Mutation
                    mutated_chromosome(1, :) = ...
                        app.CGA_Mutation(crossovered_chromosome(1, :));
                    mutated_chromosome(2, :) = ...
                        app.CGA_Mutation(crossovered_chromosome(2, :));
                    mutated_chromosome(3, :) = ...
                        app.CGA_Mutation(crossovered_chromosome(3, :));
        
                    % Apply Inversion
                    inverted_chromosome(1, :) = ...
                        app.CGA_Inversion(mutated_chromosome(1, :));
                    inverted_chromosome(2, :) = ...
                        app.CGA_Inversion(mutated_chromosome(2, :));
                    inverted_chromosome(3, :) = ...
                        app.CGA_Inversion(mutated_chromosome(3, :));
        
                    % Convert binary gens into decimal PID format.
                    PID_format = app.BinaryGenstoPID(inverted_chromosome);
                    
                    % Check for PID if 0 or LIMIT_NUMBER.
                    if PID_format(1) == 0, PID_format(1) = 1; end
                    if PID_format(1) == app.LIMIT_NUMBER, PID_format(1) = app.LIMIT_NUMBER - 1; end
                    if PID_format(2) == 0, PID_format(2) = 1; end
                    if PID_format(2) == app.LIMIT_NUMBER, PID_format(2) = app.LIMIT_NUMBER - 1; end
                    if PID_format(3) == 0, PID_format(3) = 1; end
                    if PID_format(3) == app.LIMIT_NUMBER, PID_format(3) = app.LIMIT_NUMBER - 1; end
        
                    % Assert into the population.
                    population(1:3, index + 2) = PID_format;
                end
        
                % New generation created after previous for loop.
                app.GenerationCreated = app.GenerationCreated + 1;
        
                % Put them into fitness function, and sort the new generation.
                sorted_population = app.ProcessAndSortFitness(population);
                
                % If the previous_best is the same as this best, increase
                % the number of being_same_generation.
                if (previous_best(:) == sorted_population(:, 1))
                    being_same_generation = being_same_generation + 1;
                else
                    % Update the previous best.
                    previous_best(:) = sorted_population(:, 1);
                end
            end
        
            % If the ending condition met, stop the process, and return.
            app.BestPIDResult = sorted_population([1, 2, 3], 1);
            app.BestFitnessValue = sorted_population(4, 1);
            app.LastPopulation = population;
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: SaveButton
        function SaveButtonPushed(app, event)
            app.IsImporting = false;

            app.PlantModelInputs = [app.gain_input.Value,...
                                    app.firstorder_input.Value,...
                                    app.delay_input.Value];
            
            app.PlantModel = app.CreatePlantObject(app.gain_input.Value,...
                                                   app.firstorder_input.Value,...
                                                   app.delay_input.Value);
            app.SavedStatusOutput.Text = "Saved";

            app.StarttheSimulationButton.Enable = "on";
            app.StoptheSimulationButton.Enable = "off";
        end

        % Button pushed function: StarttheSimulationButton
        function StarttheSimulationButtonPushed(app, event)
            app.StarttheSimulationButton.Enable = 'off';

            % If using the textboxes,
            if app.IsImporting == false
                % Get the data from textboxes
                app.ReferenceValue = app.referencevalue_input.Value;
                app.PopulationSize = app.populationsize_input.Value;
                app.MaxGeneration = app.maxgeneration_input.Value;
                app.FitnessValueToStop = app.fitnessbigger_input.Value;
                app.Chances(1) = 1 - app.crossoverchance_input.Value / 100;
                app.Chances(2) = 1 - app.mutationchance_input.Value / 100;
                app.Chances(3) = 1 - app.inversionchance_input.Value / 100;
    
                % Open the Stop Button
                app.StoptheSimulationButton.Enable = 'on';
    
                % Start the Algorithm
                app.StopFlag = false;
                app.ClassicGeneticAlgorithm();
    
                % Call the Ending Function
                app.SimulationEndedConditions();
            % If the details are imported
            else
                % Open the Stop Button
                app.StoptheSimulationButton.Enable = 'on';
    
                % Start the Algorithm
                app.StopFlag = false;
                app.ContinueClassicGeneticAlgorithm();

                % Call the Ending Function
                app.SimulationEndedConditions();
                
            end

            % Fix the Stop Button's text.
            app.StoptheSimulationButton.Text = "Stop the Simulation";
            
        end

        % Value changed function: gain_input
        function gain_inputValueChanged(app, event)
           if app.gain_input.Value ~= app.PlantModelInputs(1)
               app.SavedStatusOutput.Text = "Not Saved";
           end
        end

        % Value changed function: delay_input
        function delay_inputValueChanged(app, event)
            if app.delay_input.Value ~= app.PlantModelInputs(3)
               app.SavedStatusOutput.Text = "Not Saved";
           end
        end

        % Value changed function: firstorder_input
        function firstorder_inputValueChanged(app, event)
            if app.firstorder_input.Value ~= app.PlantModelInputs(2)
               app.SavedStatusOutput.Text = "Not Saved";
           end
        end

        % Button pushed function: StoptheSimulationButton
        function StoptheSimulationButtonPushed(app, event)
            app.StopFlag = true;
            app.StoptheSimulationButton.Enable = "off";
            
            app.SimulationEndedConditions();
        end

        % Button pushed function: SaveResultsasTXTButton
        function SaveResultsasTXTButtonPushed(app, event)
            filter = {'*.txt', 'Text Files (txt)'};
            [file, path] = uiputfile(filter, "Save the Results");
            full_address_file = fullfile(path, file);
            
            % Check if user clicked SAVE
            if ~isequal(file, 0) || ~isequal(path, 0)
               fileID = fopen(full_address_file, "w");
               fprintf(fileID, "==========================================\n");
               fprintf(fileID, " PID Tuner with Classic Genetic Algorithm \n");
               fprintf(fileID, "==========================================\n");
               fprintf(fileID, "Version: %s\n", app.VersionInfo);
               fprintf(fileID, "GitHub: %s\n", app.GitHubLink);
               fprintf(fileID, "Report Author: %s\n", char(java.lang.System.getProperty('user.name')));
               fprintf(fileID, "Report Created at: %s\n\n", datestr(datetime(now,'ConvertFrom','datenum')));
              
               fprintf(fileID, "######## Plant #######\n");
               fprintf(fileID, "Model: %.2f/(%.2f s + 1) * e^(-%.2f s)\n", app.PlantModelInputs(1), app.PlantModelInputs(2), app.PlantModelInputs(3));
               fprintf(fileID, "Model Inputs: %.2f, %.2f, %.2f\n", app.PlantModelInputs(1), app.PlantModelInputs(2), app.PlantModelInputs(3));
               fprintf(fileID, "######################\n\n");

               fprintf(fileID, "###### Algorithm ######\n");
               fprintf(fileID, "Population Size: %d\n", app.PopulationSize);
               fprintf(fileID, "Reference Value: %.2f\n", app.ReferenceValue);
               fprintf(fileID, "Seeking Fitness Value >=: %f\n", app.FitnessValueToStop);
               fprintf(fileID, "Maximum Same Generation: %d\n", app.MaxGeneration);
               fprintf(fileID, "Crossover Chance: %d\n", (1-app.Chances(1))*100);
               fprintf(fileID, "Mutation Chance: %d\n", (1-app.Chances(2))*100);
               fprintf(fileID, "Inversion Chance: %d\n", (1-app.Chances(3))*100);
               fprintf(fileID, "#####################\n\n");

               fprintf(fileID, "####### Results ######\n");
               fprintf(fileID, "Generation Created: %d\n", app.GenerationCreated);
               fprintf(fileID, "Best Fitness Value: %f\n", app.BestFitnessValue);
               fprintf(fileID, "Proportional Gain: %.2f\n", app.BestPIDResult(1));
               fprintf(fileID, "Integration Gain: %.4f\n", 1/app.BestPIDResult(2));
               fprintf(fileID, "Derivation Gain: %.4f\n", 1/app.BestPIDResult(3));
               fprintf(fileID, "######################\n");
               fclose(fileID);
            end
        end

        % Button pushed function: SavetheSimulationButton
        function SavetheSimulationButtonPushed(app, event)
            filter = {'*.xml', 'XML Files (xml)'};
            [file, path] = uiputfile(filter, "Save the Simulation");
            full_address_file = fullfile(path, file);
            
            % Check if user clicked SAVE
            if ~isequal(file, 0) || ~isequal(path, 0)
               fileID = fopen(full_address_file, "w");
               
               % Create the struct
               SaveStruct.PlantModelInputs = app.PlantModelInputs;
               SaveStruct.LIMIT_NUMBER = app.LIMIT_NUMBER;
               SaveStruct.Chances = app.Chances;
               SaveStruct.MaxGeneration = app.MaxGeneration;
               SaveStruct.ReferenceValue = app.ReferenceValue;
               SaveStruct.FitnessValueToStop = app.FitnessValueToStop;
               SaveStruct.PopulationSize = app.PopulationSize;
               SaveStruct.StopFlag = app.StopFlag;
               SaveStruct.BestPIDResult = app.BestPIDResult;
               SaveStruct.BestFitnessValue = app.BestFitnessValue;
               SaveStruct.GenerationCreated = app.GenerationCreated;
               SaveStruct.LastPopulationP = app.LastPopulation(1,:);
               SaveStruct.LastPopulationI = app.LastPopulation(2, :);
               SaveStruct.LastPopulationD = app.LastPopulation(3, :);
               
               % Write it to file
               writestruct(SaveStruct, full_address_file);

               fclose(fileID);
            end
        end

        % Button pushed function: ImportaSimulationButton
        function ImportaSimulationButtonPushed(app, event)
            filter = {'*.xml', 'XML Files (xml)'};
            [file, path] = uigetfile(filter, "Import the Simulation");
            full_address_file = fullfile(path, file);

            % Check if user clicked OPEN
            if ~isequal(file, 0) || ~isequal(path, 0)
               OpenStruct = readstruct(full_address_file);
               
               % Create the struct
               app.PlantModel = app.CreatePlantObject(OpenStruct.PlantModelInputs(1), OpenStruct.PlantModelInputs(2), OpenStruct.PlantModelInputs(3));
               app.PlantModelInputs = OpenStruct.PlantModelInputs;
               app.LIMIT_NUMBER = OpenStruct.LIMIT_NUMBER;
               app.Chances = OpenStruct.Chances;
               app.MaxGeneration = OpenStruct.MaxGeneration;
               app.ReferenceValue = OpenStruct.ReferenceValue;
               app.FitnessValueToStop = OpenStruct.FitnessValueToStop;
               app.PopulationSize = OpenStruct.PopulationSize;
               app.StopFlag = OpenStruct.StopFlag;
               app.BestPIDResult = OpenStruct.BestPIDResult;
               app.BestFitnessValue = OpenStruct.BestFitnessValue;
               app.GenerationCreated = OpenStruct.GenerationCreated;
               app.LastPopulation = [OpenStruct.LastPopulationP; OpenStruct.LastPopulationI; OpenStruct.LastPopulationD];
            end
            
            % Prepare the environment.
            app.gain_input.Value = app.PlantModelInputs(1);
            app.firstorder_input.Value = app.PlantModelInputs(2);
            app.delay_input.Value = app.PlantModelInputs(3);
            app.SavedStatusOutput.Text = "Saved";
            app.populationsize_input.Value = app.PopulationSize;
            app.referencevalue_input.Value = app.ReferenceValue;
            app.fitnessbigger_input.Value = app.FitnessValueToStop;
            app.maxgeneration_input.Value = app.MaxGeneration;
            app.crossoverchance_input.Value = (1 - app.Chances(1))*100;
            app.mutationchance_input.Value = (1 - app.Chances(2))*100;
            app.inversionchance_input.Value = (1 - app.Chances(3))*100;

            app.StarttheSimulationButton.Enable = "on";
            app.SimulationEndedConditions();
            app.IsImporting = true;

        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create ClassicGeneticAlgorithmforPIDTunningUIFigure and hide until all components are created
            app.ClassicGeneticAlgorithmforPIDTunningUIFigure = uifigure('Visible', 'off');
            app.ClassicGeneticAlgorithmforPIDTunningUIFigure.NumberTitle = 'on';
            app.ClassicGeneticAlgorithmforPIDTunningUIFigure.AutoResizeChildren = 'off';
            app.ClassicGeneticAlgorithmforPIDTunningUIFigure.Position = [100 100 640 480];
            app.ClassicGeneticAlgorithmforPIDTunningUIFigure.Name = 'Classic Genetic Algorithm for PID Tunning';
            app.ClassicGeneticAlgorithmforPIDTunningUIFigure.Resize = 'off';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.ClassicGeneticAlgorithmforPIDTunningUIFigure);
            app.GridLayout.ColumnWidth = {'50x', '25x', '25x'};
            app.GridLayout.RowHeight = {'5x', '5x', '100x'};

            % Create TabGroup
            app.TabGroup = uitabgroup(app.GridLayout);
            app.TabGroup.AutoResizeChildren = 'off';
            app.TabGroup.Layout.Row = [1 3];
            app.TabGroup.Layout.Column = 1;

            % Create PlantModelSettingsTab
            app.PlantModelSettingsTab = uitab(app.TabGroup);
            app.PlantModelSettingsTab.AutoResizeChildren = 'off';
            app.PlantModelSettingsTab.Title = 'Plant Model Settings';

            % Create UnderlineTextPlantModel
            app.UnderlineTextPlantModel = uilabel(app.PlantModelSettingsTab);
            app.UnderlineTextPlantModel.FontSize = 25;
            app.UnderlineTextPlantModel.Position = [20 352 259 32];
            app.UnderlineTextPlantModel.Text = '__________________';

            % Create PlantModelText
            app.PlantModelText = uilabel(app.PlantModelSettingsTab);
            app.PlantModelText.FontSize = 25;
            app.PlantModelText.Position = [20 366 138 32];
            app.PlantModelText.Text = 'Plant Model';

            % Create DescriptionTextPlantModel
            app.DescriptionTextPlantModel = uilabel(app.PlantModelSettingsTab);
            app.DescriptionTextPlantModel.Position = [20 296 274 42];
            app.DescriptionTextPlantModel.Text = {'You have to set the necesarry coefficients for your'; 'plant. This plant has to be a first-order with inertia'; 'model.'};

            % Create UnderlineFractionTextPlantModel
            app.UnderlineFractionTextPlantModel = uilabel(app.PlantModelSettingsTab);
            app.UnderlineFractionTextPlantModel.HorizontalAlignment = 'center';
            app.UnderlineFractionTextPlantModel.Position = [66 222 74 24];
            app.UnderlineFractionTextPlantModel.Text = '__________';

            % Create gain_input
            app.gain_input = uieditfield(app.PlantModelSettingsTab, 'numeric');
            app.gain_input.ValueChangedFcn = createCallbackFcn(app, @gain_inputValueChanged, true);
            app.gain_input.HorizontalAlignment = 'center';
            app.gain_input.Position = [87 230 31 22];

            % Create DotTextPlantModel
            app.DotTextPlantModel = uilabel(app.PlantModelSettingsTab);
            app.DotTextPlantModel.FontSize = 20;
            app.DotTextPlantModel.Position = [139 220 25 26];
            app.DotTextPlantModel.Text = '.';

            % Create delay_input
            app.delay_input = uieditfield(app.PlantModelSettingsTab, 'numeric');
            app.delay_input.Limits = [1 Inf];
            app.delay_input.ValueChangedFcn = createCallbackFcn(app, @delay_inputValueChanged, true);
            app.delay_input.Position = [169 222 31 22];
            app.delay_input.Value = 1;

            % Create NegativeTextPlantModel
            app.NegativeTextPlantModel = uilabel(app.PlantModelSettingsTab);
            app.NegativeTextPlantModel.Position = [162 222 25 22];
            app.NegativeTextPlantModel.Text = '-';

            % Create EulerTextPlantModel
            app.EulerTextPlantModel = uilabel(app.PlantModelSettingsTab);
            app.EulerTextPlantModel.Position = [153 216 11 22];
            app.EulerTextPlantModel.Text = 'e';

            % Create sTextPlantModel
            app.sTextPlantModel = uilabel(app.PlantModelSettingsTab);
            app.sTextPlantModel.Position = [201 222 25 22];
            app.sTextPlantModel.Text = 's';

            % Create firstorder_input
            app.firstorder_input = uieditfield(app.PlantModelSettingsTab, 'numeric');
            app.firstorder_input.ValueChangedFcn = createCallbackFcn(app, @firstorder_inputValueChanged, true);
            app.firstorder_input.Position = [66 201 34 22];

            % Create SPlusOneTextPlantModel
            app.SPlusOneTextPlantModel = uilabel(app.PlantModelSettingsTab);
            app.SPlusOneTextPlantModel.Position = [102 201 31 22];
            app.SPlusOneTextPlantModel.Text = 's + 1';

            % Create DescriptionTextPlantModel_2
            app.DescriptionTextPlantModel_2 = uilabel(app.PlantModelSettingsTab);
            app.DescriptionTextPlantModel_2.Position = [19 18 43 22];
            app.DescriptionTextPlantModel_2.Text = 'Status:';

            % Create SavedStatusOutput
            app.SavedStatusOutput = uilabel(app.PlantModelSettingsTab);
            app.SavedStatusOutput.FontWeight = 'bold';
            app.SavedStatusOutput.Position = [59 18 64 22];
            app.SavedStatusOutput.Text = 'Not Saved';

            % Create SaveButton
            app.SaveButton = uibutton(app.PlantModelSettingsTab, 'push');
            app.SaveButton.ButtonPushedFcn = createCallbackFcn(app, @SaveButtonPushed, true);
            app.SaveButton.Position = [186 18 100 22];
            app.SaveButton.Text = 'Save';

            % Create AlgorithmSettingsTab
            app.AlgorithmSettingsTab = uitab(app.TabGroup);
            app.AlgorithmSettingsTab.AutoResizeChildren = 'off';
            app.AlgorithmSettingsTab.Title = 'Algorithm Settings';

            % Create ClassicGeneticAlgorithmSettingsLabel
            app.ClassicGeneticAlgorithmSettingsLabel = uilabel(app.AlgorithmSettingsTab);
            app.ClassicGeneticAlgorithmSettingsLabel.FontSize = 25;
            app.ClassicGeneticAlgorithmSettingsLabel.Position = [20 352 209 60];
            app.ClassicGeneticAlgorithmSettingsLabel.Text = {'Classic Genetic'; 'Algorithm Settings'};

            % Create UnderlineTextCGAS
            app.UnderlineTextCGAS = uilabel(app.AlgorithmSettingsTab);
            app.UnderlineTextCGAS.FontSize = 25;
            app.UnderlineTextCGAS.Position = [20 337 259 32];
            app.UnderlineTextCGAS.Text = '__________________';

            % Create EnvironmentSettingsText
            app.EnvironmentSettingsText = uilabel(app.AlgorithmSettingsTab);
            app.EnvironmentSettingsText.FontSize = 14;
            app.EnvironmentSettingsText.FontWeight = 'bold';
            app.EnvironmentSettingsText.Position = [65 299 150 22];
            app.EnvironmentSettingsText.Text = 'Environment Settings';

            % Create PopulationSizeEditFieldLabel
            app.PopulationSizeEditFieldLabel = uilabel(app.AlgorithmSettingsTab);
            app.PopulationSizeEditFieldLabel.Position = [22 269 106 22];
            app.PopulationSizeEditFieldLabel.Text = 'Population Size';

            % Create populationsize_input
            app.populationsize_input = uieditfield(app.AlgorithmSettingsTab, 'numeric');
            app.populationsize_input.Limits = [10 1000];
            app.populationsize_input.RoundFractionalValues = 'on';
            app.populationsize_input.HorizontalAlignment = 'center';
            app.populationsize_input.Position = [143 269 126 22];
            app.populationsize_input.Value = 20;

            % Create ReferenceValueEditFieldLabel
            app.ReferenceValueEditFieldLabel = uilabel(app.AlgorithmSettingsTab);
            app.ReferenceValueEditFieldLabel.Position = [23 233 94 22];
            app.ReferenceValueEditFieldLabel.Text = 'Reference Value';

            % Create referencevalue_input
            app.referencevalue_input = uieditfield(app.AlgorithmSettingsTab, 'numeric');
            app.referencevalue_input.HorizontalAlignment = 'center';
            app.referencevalue_input.Position = [143 233 127 22];
            app.referencevalue_input.Value = 1;

            % Create EndingConditionText
            app.EndingConditionText = uilabel(app.AlgorithmSettingsTab);
            app.EndingConditionText.FontSize = 14;
            app.EndingConditionText.FontWeight = 'bold';
            app.EndingConditionText.Position = [79 182 122 22];
            app.EndingConditionText.Text = 'Ending Condition';

            % Create fitnessbigger_input
            app.fitnessbigger_input = uieditfield(app.AlgorithmSettingsTab, 'numeric');
            app.fitnessbigger_input.Position = [22 151 111 22];

            % Create maxgeneration_input
            app.maxgeneration_input = uieditfield(app.AlgorithmSettingsTab, 'numeric');
            app.maxgeneration_input.Limits = [1 Inf];
            app.maxgeneration_input.RoundFractionalValues = 'on';
            app.maxgeneration_input.Position = [142 151 125 22];
            app.maxgeneration_input.Value = 250;

            % Create FitnessFunctionTextCGAS
            app.FitnessFunctionTextCGAS = uilabel(app.AlgorithmSettingsTab);
            app.FitnessFunctionTextCGAS.Position = [22 130 111 22];
            app.FitnessFunctionTextCGAS.Text = 'Fitness Function >=';

            % Create MaxGenerationTextCGAS
            app.MaxGenerationTextCGAS = uilabel(app.AlgorithmSettingsTab);
            app.MaxGenerationTextCGAS.Position = [142 130 129 22];
            app.MaxGenerationTextCGAS.Text = 'Max Same Generation:';

            % Create RequiredChancesText
            app.RequiredChancesText = uilabel(app.AlgorithmSettingsTab);
            app.RequiredChancesText.FontSize = 14;
            app.RequiredChancesText.FontWeight = 'bold';
            app.RequiredChancesText.Position = [77 85 128 22];
            app.RequiredChancesText.Text = 'Required Chances';

            % Create CrossoverTextCGAS
            app.CrossoverTextCGAS = uilabel(app.AlgorithmSettingsTab);
            app.CrossoverTextCGAS.Position = [23 56 74 22];
            app.CrossoverTextCGAS.Text = 'Crossover %';

            % Create MutationTextCGAS
            app.MutationTextCGAS = uilabel(app.AlgorithmSettingsTab);
            app.MutationTextCGAS.Position = [117 56 66 22];
            app.MutationTextCGAS.Text = 'Mutation %';

            % Create InversionTextCGAS
            app.InversionTextCGAS = uilabel(app.AlgorithmSettingsTab);
            app.InversionTextCGAS.Position = [202 56 68 22];
            app.InversionTextCGAS.Text = 'Inversion %';

            % Create crossoverchance_input
            app.crossoverchance_input = uieditfield(app.AlgorithmSettingsTab, 'numeric');
            app.crossoverchance_input.Limits = [0 100];
            app.crossoverchance_input.RoundFractionalValues = 'on';
            app.crossoverchance_input.HorizontalAlignment = 'center';
            app.crossoverchance_input.Position = [23 35 67 22];
            app.crossoverchance_input.Value = 50;

            % Create mutationchance_input
            app.mutationchance_input = uieditfield(app.AlgorithmSettingsTab, 'numeric');
            app.mutationchance_input.Limits = [0 100];
            app.mutationchance_input.RoundFractionalValues = 'on';
            app.mutationchance_input.HorizontalAlignment = 'center';
            app.mutationchance_input.Position = [117 35 60 22];
            app.mutationchance_input.Value = 50;

            % Create inversionchance_input
            app.inversionchance_input = uieditfield(app.AlgorithmSettingsTab, 'numeric');
            app.inversionchance_input.Limits = [0 100];
            app.inversionchance_input.RoundFractionalValues = 'on';
            app.inversionchance_input.HorizontalAlignment = 'center';
            app.inversionchance_input.Position = [202 35 66 22];
            app.inversionchance_input.Value = 10;

            % Create StarttheSimulationButton
            app.StarttheSimulationButton = uibutton(app.GridLayout, 'push');
            app.StarttheSimulationButton.ButtonPushedFcn = createCallbackFcn(app, @StarttheSimulationButtonPushed, true);
            app.StarttheSimulationButton.Enable = 'off';
            app.StarttheSimulationButton.Layout.Row = 2;
            app.StarttheSimulationButton.Layout.Column = 2;
            app.StarttheSimulationButton.Text = 'Start the Simulation';

            % Create SimulationResultsPanel
            app.SimulationResultsPanel = uipanel(app.GridLayout);
            app.SimulationResultsPanel.AutoResizeChildren = 'off';
            app.SimulationResultsPanel.Title = 'Simulation Results';
            app.SimulationResultsPanel.Layout.Row = 3;
            app.SimulationResultsPanel.Layout.Column = [2 3];

            % Create OutputGraph
            app.OutputGraph = uiaxes(app.SimulationResultsPanel);
            title(app.OutputGraph, 'Output')
            xlabel(app.OutputGraph, 'Time')
            ylabel(app.OutputGraph, 'Amplitude')
            app.OutputGraph.XGrid = 'on';
            app.OutputGraph.YGrid = 'on';
            app.OutputGraph.Visible = 'off';
            app.OutputGraph.Position = [17 59 269 185];

            % Create SimulationStepsTextSR
            app.SimulationStepsTextSR = uilabel(app.SimulationResultsPanel);
            app.SimulationStepsTextSR.HorizontalAlignment = 'center';
            app.SimulationStepsTextSR.FontWeight = 'bold';
            app.SimulationStepsTextSR.Position = [31 337 103 22];
            app.SimulationStepsTextSR.Text = 'Simulation Steps';

            % Create BestGuessTextSR
            app.BestGuessTextSR = uilabel(app.SimulationResultsPanel);
            app.BestGuessTextSR.HorizontalAlignment = 'center';
            app.BestGuessTextSR.FontWeight = 'bold';
            app.BestGuessTextSR.Position = [193 337 72 22];
            app.BestGuessTextSR.Text = 'Best Guess';

            % Create GenerationTextSR
            app.GenerationTextSR = uilabel(app.SimulationResultsPanel);
            app.GenerationTextSR.HorizontalAlignment = 'center';
            app.GenerationTextSR.Position = [13 310 81 22];
            app.GenerationTextSR.Text = 'Generation:';

            % Create generation_output
            app.generation_output = uieditfield(app.SimulationResultsPanel, 'numeric');
            app.generation_output.Editable = 'off';
            app.generation_output.HorizontalAlignment = 'center';
            app.generation_output.Position = [95 310 58 22];

            % Create kp_output
            app.kp_output = uieditfield(app.SimulationResultsPanel, 'numeric');
            app.kp_output.ValueDisplayFormat = '%.2f';
            app.kp_output.Editable = 'off';
            app.kp_output.HorizontalAlignment = 'center';
            app.kp_output.Position = [219 310 55 22];

            % Create KpTextSR
            app.KpTextSR = uilabel(app.SimulationResultsPanel);
            app.KpTextSR.HorizontalAlignment = 'center';
            app.KpTextSR.Position = [190 310 30 22];
            app.KpTextSR.Text = 'Kp:';

            % Create FitnessValueTextSR
            app.FitnessValueTextSR = uilabel(app.SimulationResultsPanel);
            app.FitnessValueTextSR.HorizontalAlignment = 'center';
            app.FitnessValueTextSR.Position = [12 283 81 22];
            app.FitnessValueTextSR.Text = 'Fitness Value:';

            % Create fitnessvalue_output
            app.fitnessvalue_output = uieditfield(app.SimulationResultsPanel, 'numeric');
            app.fitnessvalue_output.ValueDisplayFormat = '%.4f';
            app.fitnessvalue_output.Editable = 'off';
            app.fitnessvalue_output.Position = [95 283 58 22];

            % Create ki_output
            app.ki_output = uieditfield(app.SimulationResultsPanel, 'numeric');
            app.ki_output.ValueDisplayFormat = '%.4f';
            app.ki_output.Editable = 'off';
            app.ki_output.HorizontalAlignment = 'center';
            app.ki_output.Position = [219 286 55 22];

            % Create KiTextSR
            app.KiTextSR = uilabel(app.SimulationResultsPanel);
            app.KiTextSR.HorizontalAlignment = 'center';
            app.KiTextSR.Position = [190 286 30 22];
            app.KiTextSR.Text = 'Ki:';

            % Create kd_output
            app.kd_output = uieditfield(app.SimulationResultsPanel, 'numeric');
            app.kd_output.ValueDisplayFormat = '%.4f';
            app.kd_output.Editable = 'off';
            app.kd_output.HorizontalAlignment = 'center';
            app.kd_output.Position = [219 262 55 22];

            % Create KdTextSR
            app.KdTextSR = uilabel(app.SimulationResultsPanel);
            app.KdTextSR.HorizontalAlignment = 'center';
            app.KdTextSR.Position = [190 262 30 22];
            app.KdTextSR.Text = 'Kd:';

            % Create SaveResultsasTXTButton
            app.SaveResultsasTXTButton = uibutton(app.SimulationResultsPanel, 'push');
            app.SaveResultsasTXTButton.ButtonPushedFcn = createCallbackFcn(app, @SaveResultsasTXTButtonPushed, true);
            app.SaveResultsasTXTButton.Enable = 'off';
            app.SaveResultsasTXTButton.Position = [16 14 134 22];
            app.SaveResultsasTXTButton.Text = 'Save Results as TXT';

            % Create SavetheSimulationButton
            app.SavetheSimulationButton = uibutton(app.SimulationResultsPanel, 'push');
            app.SavetheSimulationButton.ButtonPushedFcn = createCallbackFcn(app, @SavetheSimulationButtonPushed, true);
            app.SavetheSimulationButton.Enable = 'off';
            app.SavetheSimulationButton.Position = [160 14 126 22];
            app.SavetheSimulationButton.Text = 'Save the Simulation';

            % Create StoptheSimulationButton
            app.StoptheSimulationButton = uibutton(app.GridLayout, 'push');
            app.StoptheSimulationButton.ButtonPushedFcn = createCallbackFcn(app, @StoptheSimulationButtonPushed, true);
            app.StoptheSimulationButton.Enable = 'off';
            app.StoptheSimulationButton.Layout.Row = 2;
            app.StoptheSimulationButton.Layout.Column = 3;
            app.StoptheSimulationButton.Text = 'Stop the Simulation';

            % Create ImportaSimulationButton
            app.ImportaSimulationButton = uibutton(app.GridLayout, 'push');
            app.ImportaSimulationButton.ButtonPushedFcn = createCallbackFcn(app, @ImportaSimulationButtonPushed, true);
            app.ImportaSimulationButton.Layout.Row = 1;
            app.ImportaSimulationButton.Layout.Column = [2 3];
            app.ImportaSimulationButton.Text = 'Import a Simulation';

            % Show the figure after all components are created
            app.ClassicGeneticAlgorithmforPIDTunningUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = ClassicGeneticAlgorithmForPID_asMFile

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.ClassicGeneticAlgorithmforPIDTunningUIFigure)
            else

                % Focus the running singleton app
                figure(runningApp.ClassicGeneticAlgorithmforPIDTunningUIFigure)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.ClassicGeneticAlgorithmforPIDTunningUIFigure)
        end
    end
end
