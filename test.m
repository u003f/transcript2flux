
options.DNApoly = false;
options.SuperDaaaaave = true;
options.SuperDaaaaaveMax = true;
options.display_results = true;

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
% if (options.DNApoly == false)
% %     problem_list = {{'ishii', 'ecoli'}, {'holm', 'ecoli'}, {'rintala', 'yeast'}};
%     problem_list = {{'ishii', 'ecoli'}, {'ishii-protein', 'ecoli'}, {'holm', 'ecoli'}, {'rintala', 'yeast'}, {'rintala-red', 'yeast'}, {'rintala-protein', 'yeast'}};
% else
%     problem_list = {{'ishii', 'ecoli_DNApoly'}, {'holm', 'ecoli_DNApoly'}, {'rintala', 'yeast_DNApoly'}};
% end

problem_list = {{'ishii', 'ecoli'}, {'ishii-protein', 'ecoli'}, {'holm', 'ecoli'}, {'rintala', 'yeast'}, {'rintala-red', 'yeast'}, {'rintala-protein', 'yeast'}};

% select experiment type
experiment_type_list = {'sim_all', 'sim_intra', 'sim_secr'};

% method_list = {'pFBA', 'GIMME', 'iMAT', 'MADE', 'E-Flux', 'Lee-12', 'RELATCH', 'GX-FBA'};
% if (options.SuperDaaaaave == true)
%     method_list{end+1} = 'SuperDaaaaave';
% end
% if (options.SuperDaaaaaveMax == true)
%     method_list{end+1} = 'SuperDaaaaaveMax';
% end

method_list = {'SuperDaaaaaveMax', 'SuperDaaaaave'};

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

%% PLOTTING

build_figures(options)
