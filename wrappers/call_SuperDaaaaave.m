function fluxes = call_SuperDaaaaave(model, gene_names, gene_exp, gene_exp_sd, MaxGrowth)

% Improved version of the Daaaaave (or "Lee et al, 2012")
%
% INPUTS
% model - cobra model
% gene_names - genes ids
% gene_exp - gene expression
% gene_exp_sd - gene expression std
%
% OUTPUTS
% fluxes - flux distribution
%
% u003f, 2015

if nargin < 5
    MaxGrowth = false;
end

% fill in zero entries
unmeasured = setdiff(model.genes, gene_names);
gene_names = [gene_names; unmeasured];
gene_exp = [gene_exp; zeros(length(unmeasured), 1)];
if isempty(gene_exp_sd)
    gene_exp_sd = zeros(size(gene_exp));
    disp('warning: no gene expression std given')
else
    gene_exp_sd = [gene_exp_sd; zeros(length(unmeasured), 1)];
end
gene_min = min(gene_exp_sd(gene_exp_sd>0));
if isempty(gene_min)
    gene_min = min(gene_exp(gene_exp>0));
end
gene_exp_sd(gene_exp_sd == 0) = gene_min/2;

% replace characters in gene identifiers
gene_names = strrep(gene_names, '-', '_');
model.grRules = strrep(model.grRules, '-', '_');

S = model.S;
[nS, nR] = size(S);
L = model.lb;
U = model.ub;
fMaxGrowth = model.c;
f = zeros(size(fMaxGrowth));
b = model.b;
csense = repmat('E', 1, nS);
vartype = repmat('C', 1, nR);

% remove any pseudo-infinites
L(L<-500) = -inf;
U(U>500) = inf;

LPproblem = {};
LPproblem.('A') = S;
LPproblem.('b') = b;
LPproblem.('lb') = L;
LPproblem.('ub') = U;
LPproblem.('osense') = -1;
LPproblem.('csense') = csense;
LPproblem.('c') = fMaxGrowth;
solution = solveCobraLP(LPproblem);
statMaxGrowth = solution.stat;
solnMaxGrowth = zeros(nR,1);
objMaxGrowth = 0;
if statMaxGrowth == 1
    % use max growth if SuperDaaaaave does not converge
    solnMaxGrowth = solution.full;
    objMaxGrowth = floor(solution.obj/eps)*eps; % round down a touch
end

MILPproblem = {};
MILPproblem.('osense') = -1;
MILPproblem.('x0') = [];
params = {}; params.printLevel = 0;

changeCobraSolverParams('MILP','timeLimit', 5*60); % max 5 mins / solve
changeCobraSolverParams('MILP','printLevel', 3); % display all

% create positive FBA problem
S = [S, sparse(nS, 2*nR); speye(nR), -speye(nR), speye(nR)];
L = [L; zeros(2*nR, 1)];
U = [U; inf(2*nR, 1)];
b = [b; zeros(nR, 1)];
f = [f; zeros(2*nR, 1)];
csense = [csense, repmat('E', 1, nR)];
vartype = [vartype, repmat('C', 1, 2*nR)];

% only allow positive or negative flux
M = 1e3*max(abs([U(isfinite(U)); L(isfinite(L)); gene_exp]));
S = [S, sparse(nS+nR, 2*nR)];
L = [L; -inf(2*nR, 1)];
U = [U; inf(2*nR, 1)];
f = [f; zeros(2*nR, 1)];
vartype = [vartype, repmat('B', 1, 2*nR)];

% p <= M * kP -> -p + M*kP >= 0
S = [S; sparse(nR, nR), -speye(nR), sparse(nR, nR), M*speye(nR), sparse(nR, nR)];
b = [b; zeros(nR, 1)];
csense = [csense, repmat('G', 1, nR)];

% n <= M * kN -> -n + M*kN >= 0
S = [S; sparse(nR, nR), sparse(nR, nR), -speye(nR), sparse(nR, nR), M*speye(nR)];
b = [b; zeros(nR, 1)];
csense = [csense, repmat('G', 1, nR)];

% kP + kN = 1
S = [S; sparse(nR, nR), sparse(nR, nR), sparse(nR, nR), speye(nR), speye(nR)];
b = [b; ones(nR, 1)];
csense = [csense, repmat('E', 1, nR)];

% abs(v) variables
S = [S, sparse(nS+4*nR, nR)];
L = [L; -inf(nR, 1)];
U = [U; inf(nR, 1)];
f = [f; zeros(nR, 1)];
vartype = [vartype, repmat('C', 1, nR)];
S = [S; sparse(nR, nR), -speye(nR), -speye(nR), sparse(nR, nR), sparse(nR, nR), speye(nR)];
b = [b; zeros(nR, 1)];
csense = [csense, repmat('E', 1, nR)];
abs_flux_index = containers.Map(model.rxns, 5*nR+1:6*nR);

% add scaling parameter a
scale_index = 6*nR + 1;
L = [L; 0];
U = [U; inf];
f = [f; 0];
vartype = [vartype, 'C'];

% v >= a L -> v - a L >= 0
index = find(not(ismember(L(1:nR), [-inf, 0, inf])));
for i=index(:)';
    NS = size(S, 1) + 1;
    S(NS, i) = 1; %#ok<SPRIX>
    S(NS, scale_index) = - L(i); %#ok<SPRIX>
end
L(index) = -inf;
b = [b; zeros(length(index), 1)];
csense = [csense, repmat('G', 1, length(index))];

% v <= a L -> - v + a L >= 0
index = find(not(ismember(U(1:nR), [-inf, 0, inf])));
for i=index(:)';
    NS = size(S, 1) + 1;
    S(NS, i) = - 1; %#ok<SPRIX>
    S(NS, scale_index) = U(i); %#ok<SPRIX>
end
U(index) = inf;
b = [b; zeros(length(index), 1)];
csense = [csense, repmat('G', 1, length(index))];

% add slack variables for genes
[nS_all, nR_all] = size(S);
nG = length(gene_names);

S = [S, sparse(nS_all, nG)];
L = [L; zeros(nG, 1)];
U = [U; inf(nG, 1)];
f = [f; zeros(nG, 1)];
vartype = [vartype, repmat('C', 1, nG)];

% add genes
S = [S; sparse(nG, nR_all), speye(nG)];
b = [b; gene_exp(:)];
csense = [csense, repmat('E', 1, nG)];
gene_index = containers.Map(gene_names, nS_all+(1:nG));

% add gene associations
gene_exp_sd = containers.Map(gene_names, gene_exp_sd);

rxn_index = containers.Map();
rxn_std = containers.Map();

for ind_rxn = 1:length(model.rxns)
    association = model.grRules{ind_rxn};
    
    if association
        
        rxn_id = model.rxns{ind_rxn};
        [nS_all, nR_all] = size(S);
        rxn_row = nS_all + 1;
        rxn_col = nR_all + 1;
        
        % add reaction
        S(rxn_row, rxn_col) = -1; %#ok<SPRIX>
        L = [L; -inf]; %#ok<AGROW>
        U = [U; inf]; %#ok<AGROW>
        f = [f; 1]; %#ok<AGROW>
        vartype = [vartype, 'C']; %#ok<AGROW>
        b = [b; 0]; %#ok<AGROW>
        csense = [csense, 'E']; %#ok<AGROW>
        
        rxn_index(rxn_id) = rxn_col;
        
        association = char(expand(sym(association))); % to DNF form
        list_of_ors = regexp(association, ' or ', 'split');
        
        std_out = 0;
        
        for ind_or = 1:length(list_of_ors)
            
            % add complex
            [~, nR_all] = size(S);
            complex_col = nR_all + 1;
            S(rxn_row, complex_col) = 1; %#ok<SPRIX>
            L = [L; 0]; %#ok<AGROW>
            U = [U; inf]; %#ok<AGROW>
            f = [f; 0]; %#ok<AGROW>
            vartype = [vartype, 'C']; %#ok<AGROW>
            
            association_or = list_of_ors{ind_or};
            list_of_ands = regexp(association_or, ' and ', 'split');
            
            std_in = inf;
            
            for ind_and = 1:length(list_of_ands)
                
                [nS_all, nR_all] = size(S);
                gene_col = nR_all + 1;
                ineq_row = nS_all + 1;
                % complex > gene -> gene - complex > 0
                S(ineq_row, gene_col) =1; %#ok<SPRIX>
                S(ineq_row, complex_col) = -1; %#ok<SPRIX>
                L = [L; 0]; %#ok<AGROW>
                U = [U; inf]; %#ok<AGROW>
                f = [f; 0]; %#ok<AGROW>
                vartype = [vartype, 'C']; %#ok<AGROW>
                b = [b; 0]; %#ok<AGROW>
                csense = [csense, 'G']; %#ok<AGROW>
                
                gene = list_of_ands{ind_and};
                index = gene_index(gene);
                S(index, gene_col) = 1; %#ok<SPRIX>
                
                std = gene_exp_sd(gene);
                std_in = min([std_in, std]);
            end
            
            std_out = sqrt(std_out^2 + std_in^2);
        end
        
        rxn_std(rxn_id) = std_out;
    end
end

if (statMaxGrowth == 1) && (MaxGrowth)
    % set max growth as constraint
    fMaxGrowthBig = zeros(size(f));
    fMaxGrowthBig(1:length(fMaxGrowth)) = fMaxGrowth;
    fMaxGrowthBig(scale_index) = -objMaxGrowth;
    S = [S; sparse(fMaxGrowthBig')];
    b = [b; 0];
    csense = [csense, 'E'];
end

MILPproblem.('A') = S;
MILPproblem.('b') = b;
MILPproblem.('lb') = L;
MILPproblem.('ub') = U;
MILPproblem.('csense') = csense;
MILPproblem.('vartype') = vartype;
MILPproblem.('c') = f;
solution = solveCobraMILP(MILPproblem, params);

% set objective as constraint
obj = solution.obj;
obj = floor(obj/eps)*eps; % round down a touch
S = [S; sparse(f')];
b = [b; obj];
csense = [csense, 'E'];

% minimise distance from data to flux
f = zeros(size(f));
for ind = 1:nR
    rxn_id = model.rxns{ind};
    if ismember(rxn_id, rxn_index.keys())
        [nS_all, nR_all] = size(S);
        % R - D = P - N -> R - D - P + N = 0
        rxn_ind = nS_all + 1;
        S(rxn_ind, abs_flux_index(rxn_id)) = 1; %#ok<SPRIX>
        S(rxn_ind, rxn_index(rxn_id)) = -1; %#ok<SPRIX>
        S(rxn_ind, nR_all+1) = -1; %#ok<SPRIX>
        S(rxn_ind, nR_all+2) = 1; %#ok<SPRIX>
        L = [L; 0; 0]; %#ok<AGROW>
        U = [U; inf; inf]; %#ok<AGROW>
%         std = rxn_std(rxn_id); f = [f; -1/std; -1/std]; %#ok<AGROW>
        f = [f; -1; -1]; %#ok<AGROW>
        vartype = [vartype, 'C', 'C']; %#ok<AGROW>
        b = [b; 0]; %#ok<AGROW>
        csense = [csense, 'E']; %#ok<AGROW>
    end
end

MILPproblem.('A') = S;
MILPproblem.('b') = b;
MILPproblem.('lb') = L;
MILPproblem.('ub') = U;
MILPproblem.('csense') = csense;
MILPproblem.('vartype') = vartype;
MILPproblem.('c') = f;
solution = solveCobraMILP(MILPproblem, params);
% fprintf('MILP solution:\t%g\n', solution.obj);
soln = solution.full;

% rescale
if not(isempty(soln))
    fluxes = soln(1:nR)/soln(scale_index);
else
    fluxes = solnMaxGrowth;
end
