![](https://raw.githubusercontent.com/electricalgorithm/ClassicGeneticAlgorithm4PID/main/doc/cga4pid Banner x2.png)

---

It is a library and a GUI application for Classic Genetic Algorithm on the topic PID tuning. It includes functions for analyzing systems, simulating a system consists of PID controller and a plant, and of course Classic Genetic Algorithm.

Classic Genetic Algorithm is a evolutionary programming methodology that does cross-over, mutation and inversion into the each generation to find a solution for particular problems. In our program, problem is finding a controller for given plant. The program firstly creates a initial population from random PID triples. After selection of parents, it applies crossover, mutation and inversion.

Selection of parent chromosomes is done by Elite Selection Method. It is a way to select using a wheel of luck, and put the best ones into the next generation. This method is useful for not losing the best solution that we have.

## Graphical User Interface

It is using the functionality of application interface, with embedded functions inside the GUI code. It is possible to watch how GUI looks in [YouTube](https://youtu.be/aWHEwL5mfgc). Here you can see the features:

* First-order inertia plant creation,
* Inputs for population size, reference value, ending conditions, and chances of biological functions,
* An option to stop simulation,
* Create a text file as a report,
* Save the simulation, and import it to continue later.

<img src="https://raw.githubusercontent.com/electricalgorithm/ClassicGeneticAlgorithm4PID/main/doc/ScreenShot.png" style="zoom:80%;" />

Please note, "Max Same Generation" text-box is an input for the maximum generation of being same "best" individual. So, if your simulation won't stop after `X` amount of generation, don't surprise.  It is going to be finished after the same individual became best for `X` generations.

## Application Interface

You can find the API in the `src/CLI` directory. To start using it, check out a fast-start example in below:

```matlab
% Create the Plant Model
plant_model = CreatePlantObject(1.5, 12, 4);
% Find the PID Tuning
[PIDCoeff, FitnessVal, GenCount] = ClassicGeneticAlgorithm(plant_model, 20, 5, -0.4, 25, 0.5, 0.5, 0.5);
% Simulate the Controller
[y_data, x_data] = SimulatePIDSystem(5, PIDCoeff, plant_model);
plot(x_data, y_data);
```

## License

Application and its functions have licensed with **MIT**. You can directly change the code, update it, or fix something, and publish it. Please contribute to this project to make it better, and free for everyone.
