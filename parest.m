% Data
V1=Vdata1;
V2=Vdata2;

% Resampling the data to produce data uniformly sampled in time
% fs = 1/mean(diff(p(:,1)));
% [y, Ty] = resample(p(:,2),p(:,1),fs);

% Creating matlab data object
data = iddata([V1;V2],[],1/0.1721)

% Setting up grey-box model
% free parameters
rc=1;
Kc=1;
lc=1;
rr=1;
Kr=1;
lr=1;
order         = [2 0 2];
parameters    = [rc,Kc,lc,rr,Kr,lr];
initial_states = [1;1];
Ts            = 0;
nonlinear_model = idnlgrey('ode',order,parameters,initial_states,Ts);

% Setting up parameter bounds
nonlinear_model.Parameters(1).Minimum = 0;
nonlinear_model.Parameters(1).Maximum = 1;
nonlinear_model.Parameters(2).Minimum = 0;
nonlinear_model.Parameters(2).Maximum = 1;
nonlinear_model.Parameters(3).Minimum = 0;
nonlinear_model.Parameters(3).Maximum = 1

% Settings for fitting method
opt = nlgreyestOptions('SearchMethod','grad');
opt.SearchOption.Tolerance = 1e-3;
opt.Display='Full';
nonlinear_model = setinit(nonlinear_model, 'Fixed', {false false}); % Estimate the initial states.
% Fit The model
nonlinear_model = nlgreyest(data,nonlinear_model,opt);

% Show fitted parameters
rc=nonlinear_model.Parameters(1).Value
Kc=nonlinear_model.Parameters(2).Value
lc=nonlinear_model.Parameters(3).Value
rr=nonlinear_model.Parameters(4).Value
Kr=nonlinear_model.Parameters(4).Value
lr=nonlinear_model.Parameters(4).Value

% Plot data and ODE fit
compare(data,nonlinear_model);