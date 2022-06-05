function PID_Coefficents = BinaryGenstoPID(CombinedGen)
% BINARYGENSTOPID This function converts a Binary Gen into float
% number format.

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