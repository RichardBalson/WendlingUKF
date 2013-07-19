% script created by Richard Balson 21/02/2013

% description
% ~~~~~~~~~~~
% this script assignes all static variables, all staes and parameters are
% initialised by realising a gaussin distribution with an assigned number
% of standard deviations.

% last edit
% ~~~~~~~~~


% next edit
% ~~~~~~~~~

% Beginning of script

% ~~~~~~~~~~~~~~~~~~~~~

% Image names and handling
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~`
 sampling_frequency =fs;
 frequency_limits=[30,150];
 a = 100;
 b= 30;
 g=350;
 
 Con = 135;
 
 C= [Con,0.8*Con,0.25*Con,0.25*Con,0.3*Con,0.1*Con,0.8*Con];
 
    simulation_initial_name = [Estimation_Type,'_UKF_WM8_f=',int2str(sampling_frequency),...
        'Data']; % Initaite name for simulation, used for saving purposes

% Assign standard deviation for all variables
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

load VarianceSim meanX stdX

mean_state = (mean(meanX,2))'; % Deterine mean based on numerous previous simulation means

State_std_deviation = (sum(stdX,2)/(size(stdX,2)-1))';

Base_state_uncertainty = (State_std_deviation.^2)*sqrt(dt); % Determine base uncertainty for states


% Determine standard deviation for all slow states
% ~~~~~~~~~~~~~~~~~~~~~~~~~~

State_sigma = State_std_deviation/number_of_sigma;

% Model Input
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Input_mean = (frequency_limits(2)+frequency_limits(1))/2; % Mean input frequency from pyramidal neruons external to model

% Simulated Signal Data mV
% ~~~~~~~~~~~~~~~~~~
EstStart_Sample = round(EstStart*sampling_frequency)+1;
Y = Y(EstStart_Sample:end); % Assign a new variable as the simulated output, data point form time EstStart are used as the beginning of the observations
Number_of_observations = length(Y); % Parameter that defines he number of observations from the simulation

% Number of sigma points
% ~~~~~~~~~~~~~~~~~~~

if kappa > 0
    Sigma_points = 2*Dx+1;  % Specify number of sigma points to generate, if kappa is greater than zero than the mean is propagated
    % through the system and there is therefore one more sigma point generated
else
    Sigma_points = 2*Dx; % When kappa is set to zero sigma points are generated, note there are twice the number of sigma points as there are states
    % the reason for this is that for each mean state
    % value, a sigma point one standard deviation from
    % its mean are propagated through the system. This
    % inclue the sigma points: mean - standard deviation
    % and mean + standard deviation
end

% Estimation states uncertainty
% ~~~~~~~~~~~~~~~~~~~~~~~~~

% Uncertianty for each fast state
% mV

State_uncertainty = ones(1,8).*Base_state_uncertainty; % Specify base
% state uncertainty for all states
% State_uncertainty(1,[2 6]) = (ones(1,2)*stochastic*Variable_state_uncertainty+State_uncertainty(1,[2 6])); % Alter uncertainty of parameter affected directly by stochastic input

% Uncertainty for slow states
% ~~~~~~~~~~~~~~~~~~~~~~~

if Dp >0
    Parameter_uncertainty(1,1) = Base_parameter_uncertainty + Variable_parameter_uncertainty; % Variance to allow for model error and stochastic input effect
    if Dp >1
        Parameter_uncertainty(1,2) = (Base_parameter_uncertainty + Variable_parameter_uncertainty);
        if Dp >2
            Parameter_uncertainty(1,3) = (Base_parameter_uncertainty + Variable_parameter_uncertainty); %[m*1e-1+1e-3 m*2.5e-1+1e-3 m*2.5e-1+1e-3]
        end
    end
else
    Parameter_uncertainty = [];
end

% Uncertainty for model input mean
% ~~~~~~~~~~~~~~~~~

if Dk==1
    Input_uncertainty = Variable_input_uncertainty+Base_input_uncertainty; % Define the uncertainty of the input to the model, 
    Input_uncertainty = Input_uncertainty.^2.*uncertainty_adjustment;
else
    Input_uncertainty =[];
end

% Matrix manipulation
% ~~~~~~~~~~~~~~~~~~~~~~~~

State_sigma = repmat(State_sigma,Ds,1);

State_uncertaintyF = repmat(State_uncertainty,Ds,1);

Parameter_uncertaintyF = repmat(Parameter_uncertainty,Dp,1);% Define a parameter uncertainty matrix that can easily be used for future purposes.


%%%%%%%%%%%%%%%%%%%%%%%%%
% Intialise all parameters
%%%%%%%%%%%%%%%%%%%%%%%%%

tcon =[a,b,g];

X = zeros(Dx,Number_of_observations); % Intialise state estimate matrix
Pxx = zeros(Dx,Dx,Number_of_observations);% Intialise state covariance estimate matrix
Xout = zeros(Dx,Sigma_points,Number_of_observations); % Intialise output states from the point sin the unscented transform matrix
Yout = zeros(Dy,Sigma_points,Number_of_observations); % Intialsie the model output for each set of sigm points matrix
ExpX = zeros(Dx,Number_of_observations);% Intialise athe expected states from the unscented transform matrix
ExpY = zeros(Dy,Number_of_observations);% Intialise the expected states from the unscented transform matrix
Sigma = zeros(Dx,Sigma_points,Number_of_observations);% Intialise the sigma points matrix

State_covariance_matrix = eye(Ds).*State_sigma.^2; % State covariance matrix

State_uncertainty_matrix = eye(Ds).*State_uncertaintyF.^2.*uncertainty_adjustment;

Parameter_uncertainty_matrix = eye(Dp).*Parameter_uncertaintyF.^2.*uncertainty_adjustment;

R = Observation_uncertainty^2;

Q = blkdiag(State_uncertainty_matrix,Input_uncertainty,Parameter_uncertainty_matrix);
