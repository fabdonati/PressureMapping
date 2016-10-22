classdef CenterlineSpline
    % Class implemented to handle the centerline of a vessel with a spline
    %
    % Revision history:
    % 20/12/12: first version
    % TODO: introduce this class every time the spline is used through the
    % code of PC-MRI
    
    properties
        F
        nPointsSpline
        LengthSpline
        PhysLenPts
        tparam
        Directory
    end
    
    methods
        % Initialization:
        function obj = DefaultInitialization(obj)
            obj.F = [];
            obj.nPointsSpline = NaN;
            obj.LengthSpline  = NaN;
            obj.PhysLenPts    = NaN;
            obj.Directory     = '';
        end
        function obj = SetSpline(obj,F,nPointsSpline,LengthSpline)
            obj.F = F;
            obj.nPointsSpline = nPointsSpline;
            obj.LengthSpline  = LengthSpline;
            spacing = (LengthSpline)/(nPointsSpline-1);
            obj.PhysLenPts = 0:spacing:LengthSpline;
            m0 = 0;  m = LengthSpline;
            step=(m-m0)/(nPointsSpline-1);
            obj.tparam = m0:step:m;
        end
        function obj = CenterlineSpline(CaseNumberORDirectory)
            if nargin<1
                obj = obj.DefaultInitialization();
            else
                [F,nPointsSpline,LengthSpline,PhysLenPts,RootDirectory] = LoadSplineFile(CaseNumberORDirectory);
                obj = obj.SetSpline(F,nPointsSpline,LengthSpline);
                obj.Directory     = RootDirectory;
                % Parametric description of the spline: it uses the length of the skeleton
                % as the parameter to "travel":
            end
        end
        % Evaluation:
        function Points = GetAllPoints(obj)
            %Points = ppval(       obj.F , obj.tparam );
            Points = fnval(       obj.F , obj.tparam );
        end
        function Slopes = GetAllSlopes(obj)
            %Slopes = ppval( fnder(obj.F), obj.tparam );
            Slopes = fnval( fnder(obj.F), obj.tparam );
        end
        function Point = GetPoint(obj,iSplinePoint)
            % Function to get the 3D coordinate of a point corresponding to
            % the index iSplinePoint
            if iSplinePoint<0
                fprintf('ERROR (CenterlineSpline.m)! iPoint must be >=0 (current value = %1.2f)\n',iSplinePoint);
            end
            if iSplinePoint>obj.nPointsSpline
                fprintf('ERROR (CenterlineSpline.m)! iPoint must be <=nPointsSpline=%1.2f (current value = %1.2f)\n',obj.nPointsSpline,iSplinePoint);
            end
            t = obj.LinearInterp(iSplinePoint,obj.tparam);
            %Point = ppval(obj.F,t)';
            Point = fnval(obj.F,t)';
        end
        
        function Slope = GetSlope(obj,iSplinePoint)
            % Function to get the 3D coordinate of a slope corresponding to
            % the index iSplinePoint
            if iSplinePoint<0
                fprintf('ERROR (CenterlineSpline.m)! iPoint must be >=0 (current value = %1.2f)\n',iSplinePoint);
            end
            if iSplinePoint>obj.nPointsSpline
                fprintf('ERROR (CenterlineSpline.m)! iPoint must be <=nPointsSpline=%1.2f (current value = %1.2f)\n',obj.nPointsSpline,iSplinePoint);
            end
            t = obj.LinearInterp(iSplinePoint,obj.tparam);
            %Slope = ppval(fnder(obj.F),t)';
            Slope = fnval(fnder(obj.F),t)';
        end        
        
        function Dist  = GetDistance(obj,iSplinePoint)
            % The spline is parametrised by a parameter t, which is
            % the distance in mm from the first point. On the other hand,
            % the spline is subdivided in a set of control points (dividing
            % the cubic interpolation).
            Dist = obj.LinearInterp(iSplinePoint,obj.tparam);
        end
        % Inverse evaluation
        function iPoint = GetiPointAtDist(obj,Dist)
            % Total distance:
            D = obj.LengthSpline;
            % Total number of points:
            N = obj.nPointsSpline;
            iPoint = ((N-1) * Dist/D) + 1;
        end
        %% Search:
        function iPoint = GetValidSplinePoint(obj,ImageName,options)
            % Function to find a first valid point, from both ends of the
            % vessel, where the centerline point is a valid point.
            bDebug = 0;
            [im,hd] = io_ReadMedicalImage(ImageName);
            % Start from the first point of the spline, and look for the last valid
            % point:
            FirstORLast = 1;
            if isfield(options,'FirstORLast')
                FirstORLast = options.FirstORLast;
            end
            
            switch FirstORLast
                case 1
                    bFirst = 1;
                    iP = obj.nPointsSpline + 1;
                case 2
                    bFirst = 0;
                    iP = 0;
                otherwise
                    fprintf('ERROR in GetValidSplinePoint\n');
            end
            
            bNaN = 0;
            iComparison = 1;
            while (~bNaN)&&(iComparison<=obj.nPointsSpline)
                if bFirst, iP= iP - 1;
                else       iP = iP + 1;   end
                point = obj.GetPoint(iP);
                % Get the closest value in the mask of the image loaded:
                VoxCoords = GetVoxelCoordinates(point,hd);
                if (im(VoxCoords(1),VoxCoords(2),VoxCoords(3)))==0
                    bNaN = 1;
                end
                iComparison = iComparison + 1;
            end
            iPoint  = iP;
            if(bDebug)
                figure('color',[1 1 1]);
                show_segment_surface(im,hd);
                hold on;
                obj.PlotSpline();
                obj.PlotSplinePoint(iPoint);
                title(sprintf('iPoint = %i',iPoint));
                pause;
            end
        end
        %% Draw and rendering
        function [] = PlotSpline(obj)
            hold on;
            %Ft=ppval(obj.F,obj.tparam);
            Ft=fnval(obj.F,obj.tparam);
            [a b] = size(Ft); if a<b, Ft = Ft'; a=b; end;
            plot3(Ft(:,1),Ft(:,2),Ft(:,3));
        end
        function [] = PlotSplinePoint(obj,iPoint,options)
            hold on;
            MarkerSize = 10;
            bColor = 0;
            ColorStyle = 'r*';
            bText = 0;
            text2write = sprintf('%1.1f',obj.GetDistance(iPoint));
            if nargin==3
                if isfield(options,'MarkerSize'), MarkerSize = options.MarkerSize; end
                if isfield(options,'ColorStyle'), ColorStyle = options.ColorStyle; end
                if isfield(options,'bText'), bText = options.bText; end
                if isfield(options,'text2write'), text2write = options.text2write; end
                if isfield(options,'color'), color = options.color; bColor = 1; end
            end
            Point = obj.GetPoint(iPoint);
            if(bColor)
                plot3(Point(1),Point(2),Point(3),'.','MarkerSize',30,'color',color);
            else
                plot3(Point(1),Point(2),Point(3),ColorStyle,'MarkerSize',MarkerSize);
            end
            if bText
                dt = obj.LengthSpline/50;
                
                text(Point(1),Point(2)+dt,Point(3)+dt,text2write);
            end
        end
        function [] = PlotDistances(obj,nPoints)
            % Draw a set of points alongside the length of the aorta
            nP = 10;
            if nargin >=2
                nP = nDist;
            end
            for iPoint = 1:nP
                iP = round(iPoint * obj.nPointsSpline / nP);
                optPlot.bText = 1;
                obj.PlotSplinePoint(iP,optPlot);
            end
        end
        function [] = PlotEvery2cm(obj,StartPoint,Values,bColourByValue)
            if nargin<4
                bColourByValue = 0;
            end
            % Distance between plotting points:
            D2P = 20;
            % Distance remaining:
            dist0 = obj.GetDistance(StartPoint);
            % Number of points:
            nPoints = floor((obj.LengthSpline - dist0)/D2P);
            iPoint= zeros(1,nPoints);
            value = zeros(1,nPoints);
            for iP = 1:nPoints
                Dist = dist0 + D2P*iP;
                iPoint(iP) = obj.GetiPointAtDist(Dist);
                value(iP) = obj.LinearInterp(iPoint(iP),Values);
            end
            if(bColourByValue)
                nColorPoints = 256;
                map = colormap(jet(nColorPoints));
                MaxColor = max(abs(value));
            end
            for iP = 1:nPoints
                optPlot.bText = 1;
                optPlot.text2write = sprintf('%1.2f',value(iP));
                if(bColourByValue)
                    iColor = 1 + round((nColorPoints-1) * (MaxColor + value(iP)) / (2*MaxColor));
                    optPlot.color = map(iColor,:);
                    obj.PlotSplinePoint(iPoint(iP),optPlot);
                else
                    obj.PlotSplinePoint(iPoint(iP),optPlot);
                end
            end
        end
    end
    
    methods (Static = true)
        function value = LinearInterp(iPoint,Values)
            % Linear interpolation of the values:
            iP0 = floor(iPoint);
            iP1 = ceil(iPoint);
            R = iPoint - iP0;
            value = Values(iP0)*R + Values(iP1)*(1-R);
        end
    end
end