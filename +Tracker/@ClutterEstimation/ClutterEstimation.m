classdef ClutterEstimation < handle
    %CLUTTERESTIMATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        expectedTypes = ["Classical","Spatial","Temporal"]
        xc, yc, num_x, num_y, Xc, Yc, numCells
    end
    properties
        Type            % input: ["Classical","Spatial","Temporal"]
        AveragingLength % input: Number of scans to estimate map
        MmtRegion       % input: [[xMin;xMax], [yMin;yMax], ...] 2xNumDim (m)
        CellSize        % input: [delX, delY, ...] 1xNumDim (m)

        Map     % output: 
        Grid    % output: 
        MmtIdx  % output: 
    end
    
    methods
        function obj = ClutterEstimation(varargin)
            %CLUTTERESTIMATION Construct an instance of this class
            %   Detailed explanation goes here
            DEFAULT.Type = "Classical";
            DEFAULT.AveragingLength = 20;
            DEFAULT.MmtRegion = [[-40;40],[-40;40]]*1e3;
            DEFAULT.CellSize = [2,2]*1e3;
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

            % Set the state vector, process noise, and propagator according
            % to type
            obj = update_clutter_model(obj);
        end
        
        function obj = ClutterMap(obj, zClut, curScan)
            %ClutterMap Update clutter map with chosen Type
            switch obj.Type
                case "Classical"
                    obj.Map = classical(obj, zClut, curScan);
                case "Spatial"
                    error("Spatial clutter map not implemented yet")
                case "Temporal"
                    error("Temporal clutter map not implemented yet")
                otherwise
                    warning("Defaulting to Classical clutter map method")
            end
        end

        function PlotClutterMap(obj)
            figure
            km = 1e3;
            imagesc(obj.Grid.x/km, obj.Grid.y/km, obj.Map)
            colorbar
            xlabel("\bfEast (km)")
            ylabel("\bfNorth (km)")
            title("\bf\fontsize{14}"+obj.Type+" Clutter Map")
            axis square xy
        end
    end

    methods (Access = private)
        function map = classical(obj, zClut, curScan)
            map = obj.Map;
            idx = min(curScan, obj.AveragingLength);
            alp = (idx-1)/idx; % Compute alpha outside the loop
            lx = obj.CellSize(1);
            ly = obj.CellSize(2);
            if isempty(zClut)
                lgcl_x = false(obj.numCells,1);
                lgcl_y = false(obj.numCells,1);
            else
                lgcl_x = abs(bsxfun(@minus, zClut(1,:), obj.Xc(:))) < lx/2;
                lgcl_y = abs(bsxfun(@minus, zClut(2,:), obj.Yc(:))) < ly/2;
            end
            lgcl_xy = lgcl_y & lgcl_x; 
            mu = sum(lgcl_xy, 2); % Sum over the second dimension
            [obj.MmtIdx,~] = find(lgcl_xy);
            mu = reshape(mu, [obj.num_y, obj.num_x]);
            map = alp .* map + (1-alp) .* mu/(lx*ly);%sqrt(lx*lx + ly*ly);
        end
        % Make sure the type of motion model is available
        function verify_type(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if ~contains(obj.expectedTypes, obj.Type)
                strJ = strjoin(obj.expectedTypes,", ");
                str = sprintf('%s', strJ);
                error('Wrong type. Input "Type" \"%s\" must be one of \"%s\".', ...
                    obj.Type, str)
            end
        end
        function obj = update_clutter_model(obj)
            obj.xc = obj.MmtRegion(1,1):obj.CellSize(1): ...
                obj.MmtRegion(2,1) - obj.CellSize(1);
            obj.yc = obj.MmtRegion(1,2):obj.CellSize(2): ...
                obj.MmtRegion(2,2) - obj.CellSize(2);
            obj.xc = obj.xc - mean(obj.xc);
            obj.yc = obj.yc - mean(obj.yc);
            obj.num_x = length(obj.xc);
            obj.num_y = length(obj.yc);
            [obj.Xc,obj.Yc] = meshgrid(obj.xc,obj.yc);
            obj.numCells = numel(obj.Xc);
            obj.Grid.x = obj.xc;
            obj.Grid.y = obj.yc;
            switch obj.Type
                case "Classical"
                    obj.Map = zeros(obj.num_y, obj.num_x) ...
                        + 1e-19;%*prod(obj.CellSize);
                case "Spatial"
                    error("Spatial clutter map not implemented yet")
                case "Temporal"
                    error("Temporal clutter map not implemented yet")
                otherwise
                    warning("Defaulting to Classical clutter map method")
            end
        end
    end    
end

