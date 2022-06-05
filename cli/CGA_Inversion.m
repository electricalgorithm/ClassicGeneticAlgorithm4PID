function chromosome = CGA_Inversion(chromosome, inversion_chance)
% CGA_INVERSION This function flips the gens starting from a random index
% and ending from another random index.
    
    % If the chance achieved, then do the inversion.
    active_chance = randi(100) / 100;
    if active_chance >= inversion_chance

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