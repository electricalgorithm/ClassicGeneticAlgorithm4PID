function [y_data, x_data] = SimulatePIDSystem(ref_val, PID_gains, plant_model)
    % SimulatePIDSystem A function to create a first order inertia system
    % with PID controller, and simulate it using a step function.
    % Input arguments: 
    %  - Reference Value for System, 
    %  - Gains=[P, 1/I, 1/D],
    %  - Plant Coefficients=[k, T, tau].

    % Create the PID controller with given gains.
    pid_controller = pid(PID_gains(1), 1/PID_gains(2), 1/PID_gains(3));
    % Create a system to simulate.
    feedback_system = feedback(pid_controller * plant_model, 1);
    % Simulate the system with given reference step function.
    step_options = stepDataOptions('StepAmplitude', ref_val);
    [y_data, x_data] = step(feedback_system, step_options);
end