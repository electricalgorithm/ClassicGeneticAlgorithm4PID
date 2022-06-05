function mutated_chromosome = CGA_Mutation(original_chromosome, change_chance)
% CGA_MUTATION This function changes the gens for given chance of change.
    
    mutated_chromosome = original_chromosome;

    for index = 1:size(original_chromosome, 2)
        random_chance = randi(100) / 100;

        % If the randomized chance bigger then chance of change, flip it.
        if random_chance >= change_chance
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