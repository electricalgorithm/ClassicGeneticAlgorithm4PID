function [EachCombined, BinaryGens] = PIDtoBinaryGens(PIDCoefficents)
% PIDTOBINARYGENS The function converts three float number
% into its gene forms.
% Example:
%   Input P=45,789 -> [45, 78] -> ['0101101'; '1001110']
    
    % LIMIT_NUMBER
    LIMIT_NUMBER = 128;

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
    LIMIT_NUMBER = 128;
    if Kp_arr(1) >= LIMIT_NUMBER || ...
       Kp_arr(2) >= LIMIT_NUMBER || ...
       Ki_arr(1) >= LIMIT_NUMBER || ...
       Ki_arr(2) >= LIMIT_NUMBER || ...
       Kd_arr(1) >= LIMIT_NUMBER || ...
       Kd_arr(2) >= LIMIT_NUMBER
        warning("You shouldn't put a number which has more than 7 digits" + ...
              " in binary format, or put a number that is bigger than 127.")

    end

    % Correct each setting if they are 0 or 128.
    Kp_arr(Kp_arr == 0) = Kp_arr(Kp_arr == 0) + 1;
    Kp_arr(Kp_arr == LIMIT_NUMBER) = Kp_arr(Kp_arr == LIMIT_NUMBER) - 1;
    Ki_arr(Ki_arr == 0) = Ki_arr(Ki_arr == 0) + 1;
    Ki_arr(Ki_arr == LIMIT_NUMBER) = Ki_arr(Ki_arr == LIMIT_NUMBER) - 1;
    Kd_arr(Kd_arr == 0) = Kd_arr(Kd_arr == 0) + 1;
    Kd_arr(Kd_arr == LIMIT_NUMBER) = Kd_arr(Kd_arr == LIMIT_NUMBER) - 1;

    % Convert each decimal into binary.
    Kp_bin = dec2bin(Kp_arr, floor(log2(LIMIT_NUMBER)));
    Ki_bin = dec2bin(Ki_arr, floor(log2(LIMIT_NUMBER)));
    Kd_bin = dec2bin(Kd_arr, floor(log2(LIMIT_NUMBER)));

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