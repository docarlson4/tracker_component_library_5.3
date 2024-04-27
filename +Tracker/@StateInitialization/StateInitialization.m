classdef StateInitialization
    %STATEINITIALIZATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Type
    end
    
    methods
        function obj = StateInitialization(inputArg1,inputArg2)
            %STATEINITIALIZATION Construct an instance of this class
            %   Detailed explanation goes here
            obj.Type = inputArg1 + inputArg2;

            obj = init(obj);
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Type + inputArg;
        end
    end
end

function obj = init(obj)
end