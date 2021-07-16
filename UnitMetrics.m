classdef UnitMetrics < handle
    properties
        missingRate         (1,1)   double = NaN
        rpvRate             (1,1)   double = NaN
        FWHM                (1,1)   double = NaN
        troughToPeak        (1,1)   double = NaN
        meanAC              (1,1)   double = NaN
        ACarea              (1,1)   double = NaN
        repolarizationSlope (1,1)   double = NaN
        recoverySlope       (1,1)   double = NaN
        gmFalsePos          (1,:)   double
        gmFalseNeg          (1,:)   double
        gmUIDs              (1,:)   double
        matchConfidence     (1,:)   double
    end
    
    methods
        function obj = UnitMetrics(varargin)
            allowable = fieldnames(obj);
            if mod(length(varargin),2) ~= 0
                error('Inputs must be in name, value pairs');
            end
            for v = 1:2:length(varargin)
                if find(ismember(allowable,varargin{v}))
                    obj.(varargin{v}) = varargin{v+1};
                else
                    disp([9 'Not assigning ''' varargin{v} ''': not a field in UnitMetrics']);
                end
            end
        end
        
        function fp = falsePositiveRate(obj)
            fp = max(obj.rpvRate, sum(obj.gmFalsePos));
        end
        
        function fn = falseNegativeRate(obj)
            fn = (1 - (1-obj.missingRate)) + sum(obj.gmFalseNeg);
        end
        
    end
end