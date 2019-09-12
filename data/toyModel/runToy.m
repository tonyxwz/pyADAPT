%% initialize

import AMF.*

useMEX = 1;

plotTraj = 1;
plotSim = 1;
plotProfL = 1;
plotPSA = 1;

group_num = 1;
par_fit_num = 1;

model = Model('toyModel',useMEX);
data = DataSet('toyData');

fprintf(['There are ' num2str(numel(data.groups)) ' dataset groups available. Currently showing datagroup: ' num2str(group_num)])

loadGroup(data, data.groups{group_num}); % changed from 'toy'
initiateExperiment(model, data);

model.result

%% config
% THINGS HERE WILL OVERWRITE WHAT IS ASSIGNED IN THE MODEL CLASS

model.options.useMex       = useMEX;
model.options.savePrefix   = '';
model.options.odeTol       = [1e-12 1e-12 100];
model.options.numIter      = 20;
model.options.numTimeSteps = 100;
model.options.parScale     = [2 -2];
model.options.seed         = 1;
model.options.SSTime       = 1000;
model.options.lab1         = 0.1; % .1
model.options.optimset.Display = 'off';

parseAll(model);
compileAll(model);

model.result

%% run

model.result.traj = logical([1 1 1 1 1]);

result = runADAPT(model);
fprintf('after ADAPT')
model.result

% result2 = simulate(model); % wat doe ik hiermee?
fprintf('after simulate')
model.result

%% plot

if plotTraj == 1
    figure;
    plot(result, {'s3','ds3dt','k1'});
end

if plotSim == 1
    plotAll(result, 'states', 'traj');
    plotAll(result, 'parameters', 'traj');

    v1 = getValue(result, 'v1');
    v2 = getValue(result, 'v2');
    s3 = getValue(result, 's3');
    t = result.time;

    figure;
    plot(diff(s3)./repmat(diff(t(:)), 1, size(s3,2)), 'r');
    hold on;
    plot(v1-v2, 'b');
    title('Simulation of states and reactions')
    xlabel('time')
    ylabel('[variable]')
end

if plotPSA == 1
    run_PSA(model);
end

if plotProfL == 1
    run_PLA(model, par_fit_num);
end

model.result



model.result