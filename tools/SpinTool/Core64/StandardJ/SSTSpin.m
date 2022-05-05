% Spectroscopic Simulation Tool: Spin object
% Encapsulates a single chemical site in a J-coupled Hamiltonian.
classdef SSTSpin
    properties
        pos
        rho
        cs
        eqMag
    end
    methods
        function obj = SSTSpin(varargin)
            p = inputParser;
            p.addParameter('pos', [0 0 0], @(x) validateattributes(x, {'double'}, {'ncols', 3}));
            p.addParameter('rho', [0 0 1], @(x) validateattributes(x, {'double'}, {'ncols', 3}));
            p.addParameter('cs', 0, @(x) validateattributes(x, {'double'}, {'size', [1 1]}));
            p.addParameter('eqMag', 1, @(x) validateattributes(x, {'double'}, {'size', [1 1]}));
            p.parse(varargin{:});
            
            obj.r = p.Results.r;
            obj.M = p.Results.M;
            obj.cs = p.Results.cs;
            obj.eqMag = p.Results.eqMag;
        end
    end
end