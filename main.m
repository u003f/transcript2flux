%% MAIN FILE:
%
% Use this to run all the tests in the paper.
% Different tests can be run selectively.
% For convenience all results are stored as .mat files.
%
% Author: Daniel Machado, 2013
% u003f 28-apr-15: added DNApoly option

% no result DNApoly:

% MADE / ishii / sim_intra -> error
% Lee-12 / ishii / sim_intra -> DNF

options.DNApoly = true;

addpath('tests')
addpath('utils')
addpath('wrappers')
addpath('plotting')

%initialize cobra
COBRA_SOLVER = 'gurobi5';
changeCobraSolver(COBRA_SOLVER, 'all');

%initialize cmpi (requires MADE package; see README)
cmpi.init()
cmpi.set_solver('gurobi')

%% PROBLEM SETUP

% select dataset and organism (see load_dataset.m for details)
if (options.DNApoly == false)
    problem_list = {{'ishii', 'ecoli'}, {'holm', 'ecoli'}, {'rintala', 'yeast'}};
else
    problem_list = {{'ishii', 'ecoli_DNApoly'}, {'holm', 'ecoli_DNApoly'}, {'rintala', 'yeast_DNApoly'}};
end

% select experiment type
experiment_type_list = {'sim_all', 'sim_intra', 'sim_secr'};

method_list = {'pFBA', 'GIMME', 'iMAT', 'E-Flux', 'RELATCH', 'GX-FBA', 'MADE', 'Lee-12'};

method_list = {'MADE'};
problem_list = {{'ishii', 'ecoli_DNApoly'}};
experiment_type_list = {'sim_all'};

for method_iter = method_list
    for problem_iter = problem_list
        for experiment_type_iter = experiment_type_list

            method = method_iter{1};
            problem = problem_iter{1};
            experiment_type = experiment_type_iter{1};

            DATASET = problem{1};
            ORGANISM = problem{2};

            %load model
            model = load_model(ORGANISM);

            %load data set
            dataset = load_dataset(DATASET);

            options.experiment_type = experiment_type;
            options.reestimate_data = true; % fit experimental data to model

            %% BENCHMARKING
			benchmark_method(method, model, dataset, options);

		end
	end
end

%% SENSITIVITY ANALYSIS

if (options.DNApoly == false)

    conf = cell(1,6);
    conf{1} = {'GIMME', 'OBJ_FRAC', 0, 1, 'lin'};
    conf{2} = {'GIMME', 'GIMME_LOWER_QUANTILE', 0, 1, 'lin'};
    conf{3} = {'iMAT', 'IMAT_LOWER_QUANTILE', 0, 0.75, 'lin'};
    conf{4} = {'iMAT', 'IMAT_UPPER_QUANTILE', 0.25, 1, 'lin'};
    conf{5} = {'iMAT', 'IMAT_EPS', -2, 2, 'log'};
    conf{6} = {'MADE', 'OBJ_FRAC', 0, 1, 'lin'};

    points = 100;

    for i = 1:6
        [method, parameter, min_val, max_val, scale] = deal(conf{i}{:});
        sensitivity_analysis(model, dataset, method, parameter, min_val, max_val, points, scale, options);
    end
end

%% ROBUSTNESS ANALYSIS

if (options.DNApoly == false)

    methods = {'GIMME', 'iMAT', 'E-Flux', 'Lee-12', 'RELATCH', 'MADE', 'GX-FBA'};
    steps = 2; points = 2;

    for i = 1:length(methods)
         robustness_analysis(model, dataset, methods{i}, dataset.conditions{2}, dataset.conditions{1}, steps, points, options);
    end
end

%% PLOTTING

build_figures(options)
