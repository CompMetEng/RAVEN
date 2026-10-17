function [isFeasible, isOptimal]=checkSolution(res)
% checkSolution  Check if a solver result is feasible and optimal.
%
% Parameters
% ----------
% res : struct
%     the output structure from optimizeProb.
%
% Returns
% -------
% isFeasible : logical
%     true if the solution is feasible.
% isOptimal : logical
%     true if the solution is optimal.
%
% Examples
% --------
%     [isFeasible, isOptimal] = checkSolution(res);


isFeasible=false;
isOptimal=false;

switch res.stat
    case 1
        isOptimal=true;
        isFeasible=true;
    case 2
        isOptimal=false;
        isFeasible=true;
    otherwise
        isOptimal=false;
        isFeasible=false;
end
end
