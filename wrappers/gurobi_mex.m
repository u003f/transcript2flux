function [x, val, flag, output] = gurobi_mex(c, sense, A, b, ctypes, lb, ub, vartypes, options)

% a small function to allow Made (which requires gurobi_mex) to work with 
% Gurobi 6.0 (which has no such function).
% u003f: 27-apr-15

model = struct;
model.A = A;
model.obj = c;
model.sense = ctypes;
model.rhs = b;
model.lb = lb;
model.ub = ub;
model.vtype = vartypes;
switch sense
    case -1
        model.modelsense = 'max';
    otherwise
        model.modelsense = 'min';
end

params = struct;
params.outputflag = 0; % shhhhh

result = gurobi(model, params);

if isfield(result,'x')
	x = result.x;
else
    x = [];
end

if isfield(result,'objval')
	val = result.objval;
else
    val = [];
end

switch result.status
    % from http://www.gurobi.com/documentation/6.0/refman/optimization_status_codes.html
	case 'LOADED'
		flag = 1;
	case 'OPTIMAL'
		flag = 2;
	case 'INFEASIBLE'
		flag = 3;
	case 'INF_OR_UNBD'
		flag = 4;
	case 'UNBOUNDED'
		flag = 5;
	case 'CUTOFF'
		flag = 6;
	case 'ITERATION_LIMIT'
		flag = 7;
	case 'NODE_LIMIT'
		flag = 8;
	case 'TIME_LIMIT'
		flag = 9;
	case 'SOLUTION_LIMIT'
		flag = 10;
	case 'INTERRUPTED'
		flag = 11;
	case 'NUMERIC'
		flag = 12;
	case 'SUBOPTIMAL'
		flag = 13;
	case 'INPROGRESS'
		flag = 14;
end
    
output = '';
