classdef SOT < handle
   %SOT Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = private)
        expectedTypes = ["Type1","Type2"]
    end

    properties
        Type
        Property1
    end

    methods
       function obj = SOT(varargin)
           %SOT Construct an instance of this class
            DEFAULT.Type = "Type1";
            DEFAULT.Property1 = "Property1";

            p = inputParser;

            cellfun( @(fn) p.addOptional( fn, DEFAULT.(fn) ), ...
                fieldnames( DEFAULT ) );

            p.parse(varargin{:});

            param_names = fieldnames( p.Results );
            for k = 1:numel(param_names)
                obj.(param_names{k}) = p.Results.(param_names{k});
            end

            % Make sure the type of motion model is avaialble
            verify_type(obj)
        end
    end

    methods (Access = private)
        % Make sure the type of motion model is available
        function verify_type(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if ~contains(obj.expectedTypes, obj.Type)
                strJ = strjoin(obj.expectedTypes,", ");
                str = sprintf('%s', strJ);
                error('Wrong type. Input Type \"%s\" must be one of \"%s\".', ...
                    obj.Type, str)
            end
        end
    end
end
% Developed in Matlab 9.12.0.2529717 (R2022a) Update 8 on PCWIN64.
% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2024-05-19 17:36
