# Classic Genetic Algorithm for PID Tunning
This project aims to create an GUI application using Classic Genetic Algorithm with Elite Selection Method to find PID coefficients for given first-order inertia plant.

* GUI is still in development.
* CLI is ready to use.

```matlab
% Create the Plant Model
plant_model = CreatePlantObject(1.5, 12, 4);
% Find the PID Tuning
[PIDCoeff, FitnessVal, GenCount] = ClassicGeneticAlgorithm(plant_model, 20, 5, -0.4, 25, 0.5, 0.5, 0.5);
 % Simulate the Controller
[y_data, x_data] = SimulatePIDSystem(5, PIDCoeff, plant_model);
plot(x_data, y_data);
```
