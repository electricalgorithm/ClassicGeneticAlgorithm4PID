function plant = CreatePlantObject(k, T, tau)
% CREATEPLANTOBJECT The function will create
% a first order and a time delay.
    plant = tf(k, [T 1], 'InputDelay', tau);
end