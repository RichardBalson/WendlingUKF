function Wendling_Simulation(SimulationSettings)

% script created by Richard Balson 21/02/2013

% description
% ~~~~~~~~~~~
% Simulate states and output of the extended neural mass model for
% estimation purposes. Input mean assumed to be 90.

% last edit
% ~~~~~~~~~

%

% Variables required
% ~~~~~~~~~~~~~~~

% stochastic -  Determines level of stochasticity in input
% number_of_sigma_input - Determines the standard deviation of stochastic input
% Parameter_index - Used to decide what parameter should be simulated
% Input_mean_variation - uSed to specifiy whether the mean of the input
% should be different from its described nominal value

% next edit
% ~~~~~~~~~

% Beginning of script
% ~~~~~~~~~~~~~~~~~~~~~

%%

% Model parameters(Variable parameters are to alter parameter values during
% simulation) Matrix specifies parameter values for a specified time period
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Units mV
% ~~~~~~~~~~~~

Max_parameters = SimulationSettings.simulation_time/(SimulationSettings.slope_time+0.5); % Specifiy maximum number of parameters for random simulation
[AV BV GV] = Parameter_choice(SimulationSettings,Max_parameters);

% Solver parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

simulation_changes = length(AV); % Specify the number of times the parameters change in simulation

% Units Hz
% ~~~~~~~~

sampling_frequency = SimulationSettings.fs; % Sampling frequency for solutions and measurements

%%
% Stationery Parameter defintions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~
WendlingModelStationaryParameters


% e = xxxx*randn(NSamples,1);
% Determine stochastic input to the model
% ~~~~~~~~~~~~~~~~~~~~~~~
gauss = randn(1,length(normalise)+1)*std_deviation*SimulationSettings.stochastic;

%         gauss = randn(1)*std_deviation*stochastic; % Determine random fluctuations in input signal with the specified standard deviation

normalised_gaussian_input(1,:) = meanf + gauss(1,:); % Determine value of input for all time steps

%         while (any(normalised_gaussian_input>max_frequency) || any(normalised_gaussian_input <min_frequency))
%
%            gauss(1,k) = randn(1)*std_deviation*stochastic;
%
%             normalised_gaussian_input(k) = meanf + gauss(1,k);
%
%         end

normalised_gaussian_inputSDE(1,:) = meanf + gauss(1,:)*sqrt(dt)*sampling_frequency; % Determine value of input used for calculation purposes for all time steps

%%

% Equations
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% Solver
% ~~~~~~~~~~~~~~~~~~~~~~~
% Iteration to change parameters in function
for k = normalise
    %
    % 8th Order
    % Derivation of equation shown in pdf attached
    
    % Normalised gaussian input
    
    zs1 =  z(k,2)-C(4)*z(k,3)-z(k,4);
    zs2 = C(1)*z(k,1);
    zs3 = C(3)*z(k,1);
    zs4 = C(5)*z(k,1)-C(6)*z(k,3);
    [z(k+1,1) z(k+1,5)] = PSPkernel([z(k,1) z(k,5)],dt,MVI(k,1),a,sigmoid(zs1));
    [z(k+1,2) z(k+1,6)] = PSPkernel([z(k,2) z(k,6)],dt,MVI(k,1),a,normalised_gaussian_inputSDE(k) +C(2)*sigmoid(zs2));
    [z(k+1,3) z(k+1,7)] = PSPkernel([z(k,3) z(k,7)],dt,MVI(k,2),b,sigmoid(zs3));
    [z(k+1,4) z(k+1,8)] = PSPkernel([z(k,4) z(k,8)],dt,MVI(k,3),g,C(7)*sigmoid(zs4));
    
end% End of for loop

% Determine output of model
% ~~~~~~~~~~~~~~~~~~~~~~~~
output8 = z(:,2)-C(4)*z(:,3) -z(:,4);      % 8th order output, Normalised gaussian input

save([SimulationSettings.name,'.mat'],'output8', 'z', 'normalised_gaussian_input', 'sampling_frequency', 'SimulationSettings', 'MVI', 'tcon', 'C', 'frequency_limits', 'meanf');% Save relevant parameters for future use

