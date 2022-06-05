function [overshoot_error, stability_time, error_given_time] ...
         = AnalyseSystemResult(res_t, res_y, ref_val, error_at_time)

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