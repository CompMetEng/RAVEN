function nW = parallelWorkersRAVEN(runParallel)
% parallelWorkersRAVEN  Return the parfor worker count for a given setting.
%
% Returns the value for the maximum-workers argument of parfor:
%
%     nW = parallelWorkersRAVEN(runParallel);
%     parfor (i = 1:N, nW)
%         ...
%     end
%
% The function leaves the global parallel pool unchanged. With nW == 0, parfor
% runs as a plain for loop. With nW == Inf, the Parallel Computing Toolbox
% manages the pool size and creates a pool when none exists. When
% runParallel is true but the Parallel Computing Toolbox is not installed, the
% function prints a message and returns 0.
%
% Parameters
% ----------
% runParallel : logical
%     true to use parallel workers, false for sequential execution.
%
% Returns
% -------
% nW : double
%     0   sequential execution (parfor behaves as a for loop, no pool is
%         required or created).
%     Inf available pool workers (the Parallel Computing Toolbox manages
%         the pool size and creates a pool when none exists).

if ~runParallel
    nW = 0;
    return
end
if isempty(ver('parallel'))
    disp('Cannot find MATLAB Parallel Computing Toolbox, process is not parallelized.')
    nW = 0;
else
    nW = Inf;
end
end
