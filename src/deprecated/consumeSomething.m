function [solution, metabolite]=consumeSomething(model,varargin)
% consumeSomething  DEPRECATED. Use findLeakMetabolite(model,'consume',...).
%
% Deprecated wrapper kept so existing scripts keep working; it will be
% removed in the next major release.
%
% The positional argument order is:
%
%     consumeSomething(model, ignoreMets, isNames, minNrFluxes, params, ignoreIntBounds)
%
% findLeakMetabolite takes allowExcretion (a produce-only option) in
% position 4, so the arguments are translated before they are passed on.
%
% See also
% --------
% findLeakMetabolite

deprecationWarning('consumeSomething','findLeakMetabolite(model,''consume'',...)');

p=parseRAVENargs(varargin, {'ignoreMets',[]; ...
    'isNames',false; ...
    'minNrFluxes',false; ...
    'params',[]; ...
    'ignoreIntBounds',false});

[solution, metabolite]=findLeakMetabolite(model,'consume', ...
    'ignoreMets',p.ignoreMets, ...
    'isNames',p.isNames, ...
    'minNrFluxes',p.minNrFluxes, ...
    'params',p.params, ...
    'ignoreIntBounds',p.ignoreIntBounds);
end
