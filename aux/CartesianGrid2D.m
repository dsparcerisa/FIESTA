% CartesianGrid2D class.
% Auxiliary class to store cartesian objects with spacial coordinates
%
% Constructors:
%  CartesianGrid2D(CartesianGrid2D CGobject)
% --> Copies structure from another CGrid object with empty data
%
%  CartesianGrid2D(minX,maxX,minY,maxY,dx,dy,NX,NY,data*)
% --> data matrix is optional 
% 
% Properties:
%
%  (Bin sizes in cm)
%   dx
%   dy
%
%  (Number of bins) 
%   NX
%   NY
%
%  (Extreme positions in cm at voxel center)
%   minX        
%   maxX
%   minY
%   maxY
%   data
%
% Methods:
%  rebin(newdx,newdy)--> has same size but different dx dy.
%  resize(newSize) --> has size newSize ([NX NY]) and original dx dy.
%  crop(newSize) --> Alias for resize
%  plotSlice(axis, value) --> as in plotSlice('Y',0)
%  passMinMaxTests() --> Verifies dimensions
%  displayCoords() --> Prints coordinates
%  res = posInBounds(x, y) --> Returns TRUE if position is in bounds
%  res = indexInBounds(i, j) --> Returns TRUE if indexes IJK are in bounds
%  subCG = getSubset(cropping) --> Returns a cropped CG object
%  [I,J] = obtainIndexes(inX,inY) --> maps XYZ to indexes
%  data = getData(x, y) --> Indexing by position
%  axisValues = getAxisValues(axis) --> Returns vector of axis positions
%  I = getXindex(X)
%  J = getYindex(Y)

%  X = getXvalue(I)
%  Y = getYvalue(J)
%  plotLinear(freeAxis, value, format) --> Create linear plot
%  [linData, positions] = getLinear(freeAxis, value) --> Extract linear curve
%  increaseResolution(axis, factor) --> Increase resolution in given axis
%
% Overloaded operators:
%  plus/minus --> A + b has max(A,b) size and min(A,b) in dx, dy
%                (individually for each dimension)
%                 A - b = A + (-B);
%  unary minus --> -CG
%  indexing --> CG(i,j) => CG.data(i,j);
%
% FoCa: Forward Dose and LET Calculation tool for proton radiotherapy
% D Sanchez-Parcerisa et al 2014 Phys. Med. Biol. 59 7341 doi:10.1088/0031-9155/59/23/7341
%
% Copyright (C) 2012-2015
% Author: Daniel Sanchez-Parcerisa
% Principal Investigator: Alejandro Carabe-Fernandez
% Contributors: Mark Kondrla, Chris Wang, Alex Shaindlin.
% University of Pennsylvania. All Rights Reserved.
%
% Contact details: daniel.sanchez.parcerisa@gmail.com, a.carabe@uphs.upenn.edu

classdef CartesianGrid2D < matlab.mixin.Copyable
    %CARTESIANGRID Auxiliary class to store cartesian objects with spacial coordinates
    
    properties
        % All distances are in cm
        % min/max values represent the CENTER of the voxels

        % Bin sizes (cm)
        dx
        dy
        
        % Number of bins
        NX
        NY
        
        % Extremal positions (voxel centers)
        minX        
        maxX
        minY
        maxY
        
        % Data matrix
        data
        
    end % properties
    
    methods        
        
        % Constructor
        function this = CartesianGrid2D(varargin)
            if nargin == 1
            % CartesianGrid2D(CGobject)
                    validateattributes(varargin{1}, {'CartesianGrid2D'}, {'nonempty'}, 'CartesianGrid2D Constructor', 'CG2DObject', 1);
                    CGobject = varargin{1};
                    this.dx = CGobject.dx;
                    this.dy = CGobject.dy;
                    this.NX = CGobject.NX;
                    this.NY = CGobject.NY;
                    this.minX = CGobject.minX;
                    this.maxX = CGobject.maxX;
                    this.minY = CGobject.minY;
                    this.maxY = CGobject.maxY;
                    this.passMinMaxTests;                    
                    this.data = zeros(this.NX, this.NY);  
                    
            elseif nargin == 8 || nargin == 9
            % CartesianGrid2D(minX,maxX,minY,maxY,dx,dy,NX,NY,data*)
                
                    validateattributes(varargin{1}, {'numeric'}, {'nonnan', 'real'}, 'CartesianGrid2D Constructor', 'minX', 1);
                    this.minX = varargin{1};
                    validateattributes(varargin{2}, {'numeric'}, {'nonnan', 'real'}, 'CartesianGrid2D Constructor', 'maxX', 2);
                    this.maxX = varargin{2};
                    validateattributes(varargin{3}, {'numeric'}, {'nonnan', 'real'}, 'CartesianGrid2D Constructor', 'minY', 3);
                    this.minY = varargin{3};
                    validateattributes(varargin{4}, {'numeric'}, {'nonnan', 'real'}, 'CartesianGrid2D Constructor', 'maxY', 4);
                    this.maxY = varargin{4};
                    validateattributes(varargin{5}, {'numeric'}, {'nonnan', 'nonnegative', 'real'}, 'CartesianGrid2D Constructor', 'dx', 5);
                    this.dx = varargin{5};
                    validateattributes(varargin{6}, {'numeric'}, {'nonnan', 'nonnegative','real'}, 'CartesianGrid2D Constructor', 'dy', 6);
                    this.dy = varargin{6};
                    validateattributes(varargin{7}, {'numeric'}, {'nonnan', 'positive', 'integer'}, 'CartesianGrid2D Constructor', 'NX', 7);
                    this.NX = varargin{7};
                    validateattributes(varargin{8}, {'numeric'}, {'nonnan', 'positive', 'integer'}, 'CartesianGrid2D Constructor', 'NY', 8);
                    this.NY = varargin{8};
                                        
                    this.passMinMaxTests;
                    
                    if nargin==9
                        validateattributes(varargin{9}, {'numeric'}, {'size',[this.NX this.NY]}, 'CartesianGrid2D Constructor', 'data', 9);
                        this.data = varargin{9};
                    else
                        this.data = zeros(this.NX, this.NY);
                    end
                
            elseif nargin>0
                error('Incorrect number of arguments');
            end
            
        end % constructor
                        
        % rebin(newdx, newdy) has same size but different dx dy
        
        function rebin(this,newdx,newdy)
            validateattributes(newdx, {'numeric'}, {'positive','real'}, 'CartesianGrid2D.rebin', 'newdx', 1);             
            validateattributes(newdy, {'numeric'}, {'positive','real'}, 'CartesianGrid2D.rebin', 'newdy', 2);             
            
            minXB = this.minX - this.dx/2;
            maxXB = this.maxX + this.dx/2;
            minYB = this.minY - this.dy/2;
            maxYB = this.maxY + this.dy/2;
            XposB = minXB:newdx:maxXB;
            YposB = minYB:newdy:maxYB;
            
            Xpositions = 0.5 * (XposB(1:(end-1)) + XposB(2:end));
            Ypositions = 0.5 * (YposB(1:(end-1)) + YposB(2:end));
                        
            %Xpositions = this.minX:newdx:this.maxX;
            %Ypositions = this.minY:newdy:this.maxY;
            
            dataType = class(this.data);
            GI = griddedInterpolant({this.minX:this.dx:this.maxX, ...
                this.minY:this.dy:this.maxY},double(this.data),'linear','none');
            this.data = cast(this.data,dataType);
            this.data = GI({Xpositions, Ypositions});       
            this.NX = numel(Xpositions);
            this.NY = numel(Ypositions);
            this.dx = newdx;
            this.dy = newdy;
            this.minX = Xpositions(1);
            this.minY = Ypositions(1);
            this.maxX = Xpositions(end);
            this.maxY = Ypositions(end);
            
            this.passMinMaxTests;                        
        end % rebin
                
        function shift(this, shiftVector)
            this.minX = this.minX + shiftVector(1);
            this.maxX = this.maxX + shiftVector(1);
            this.minY = this.minY + shiftVector(2);
            this.maxY = this.maxY + shiftVector(2);           
        end
       
        % resize(newSize) has size newSize ([cropping]) and original dx dy.
        
        function resize(this,newSize)
            validateattributes(newSize, {'numeric'}, {'nonnan','real','numel',4}, 'CartesianGrid2D.crop', 'newSize', 1);             
            Xpositions = newSize(1):this.dx:newSize(2);
            Ypositions = newSize(3):this.dy:newSize(4);
            dataType = class(this.data);
            GI = griddedInterpolant({this.minX:this.dx:this.maxX, ...
                this.minY:this.dy:this.maxY},double(this.data),'linear','none');
            this.data = cast(this.data,dataType);
            this.data = GI({Xpositions, Ypositions});       
            this.NX = numel(Xpositions);
            this.NY = numel(Ypositions);
            this.minX = Xpositions(1);
            this.minY = Ypositions(1);
            this.maxX = Xpositions(end);
            this.maxY = Ypositions(end);
            
            this.passMinMaxTests;
        end % resize
        
        % CG: A + b has max(A,b) size and min(A,b) in dx, dy, dz
        % (individually for each dimension)
        % A - b = A + (-B);
        
        function Result = plus(this, CGRight)
        % CG: A + b has max(A,b) size and min(A,b) in dx, dy, dz
        % (individually for each dimension)
        % A - b = A + (-B);            
            validateattributes(this, {'CartesianGrid2D'}, {'nonempty'}, 'CartesianGrid2D.plus', 'CGLeft', 1);                         
            validateattributes(CGRight, {'CartesianGrid2D'}, {'nonempty'}, 'CartesianGrid2D.plus', 'CGRight', 2);     

            minX = min(this.minX, CGRight.minX); 
            minY = min(this.minX, CGRight.minX); 
            maxX = max(this.maxX, CGRight.maxX); 
            maxY = max(this.maxY, CGRight.maxY); 
            dx = min(this.dx, CGRight.dx);
            dy = min(this.dy, CGRight.dy);
            Xpositions = minX:dx:maxX;
            Ypositions = minY:dy:maxY;
            maxX = Xpositions(end);
            maxY = Ypositions(end);
            
            GI_Left = griddedInterpolant({this.minX:this.dx:this.maxX, ...
                this.minY:this.dy:this.maxY},this.data,'linear','none');
            GI_Right = griddedInterpolant({CGRight.minX:CGRight.dx:CGRight.maxX, ...
                CGRight.minY:CGRight.dy:CGRight.maxY},CGRight.data,'linear','none');
            
            NX = numel(Xpositions);
            NY = numel(Ypositions);
            data = GI_Left({Xpositions,Ypositions}) + GI_Right({Xpositions,Ypositions});
            
            Result = CartesianGrid2D(minX,maxX,minY,maxY,dx,dy,NX,NY,data);                 
            
        end % plus        
        
        function Result = uminus(this)
            Result = this.copy;
            Result.data = -this.data;
        end % unary minus
        
        function Result = minus(this, CGRight)
            Result = this + (-CGRight);
        end % minus
        
        function Result = subsindex(this, ind)
            Result = this.data(ind);
        end % subsindex
                
        % Interface
        function crop(this, cropping)
            this.resize(cropping);
        end % crop
                                 
        % Plot a single slice
        function plotSlice(this)
            slice = this.data;
            % TODO: Test
            imagesc(this.getAxisValues('X'),fliplr(this.getAxisValues('Y')),rot90(slice));
            
            set(gca,'xlim',[this.minX-this.dx/2 this.maxX+this.dx/2]);
            set(gca,'ylim',[this.minY-this.dy/2 this.maxY+this.dy/2]);
            
            xlabel('X');
            ylabel('Y');
            
            set(gca,'Ydir','normal');
            set(gca,'Xdir','normal');
            
        end % plotSlice        
        
        % Test minima and maxima        
        function res = passMinMaxTests(this)
            res = false;
            if ( abs(this.maxX - (this.minX + this.dx*double(this.NX-1))) > 1e-6 ) || ...
                    ( abs(this.maxY - (this.minY + this.dy*double(this.NY-1))) > 1e-6 ) 
                error('Minima/maxima tests not passed.\n\tDEBUG INFO:\n\tTest X: %f\n\tTest Y: %f\n', ...
                    ( abs(this.maxX - (this.minX + this.dx*double(this.NX-1))) ), ...
                    ( abs(this.maxY - (this.minY + this.dy*double(this.NY-1))) ) )                            
            else
                res = true;
            end
        end % passMinMaxTests
        
        % Display coordinates
        function displayCoords(this)
            fprintf('## X from %g to %g \n', this.minX, this.maxX);
            fprintf('## Y from %g to %g \n', this.minY, this.maxY);
        end % displayCoords
        
        % Position in bounds
        function res = posInBounds(this, x, y)
            res = x>(this.minX - this.dx/2) & ...
                  x<(this.maxX + this.dx/2) & ...
                  y>(this.minY - this.dy/2) & ...
                  y<(this.maxY + this.dy/2);              
        end % posInBounds
        
        % Index in bounds
        function res = indexInBounds(this, i, j)
            res = i>=1 & j>=1 & i<=this.NX & j<=this.NY;
        end % indexInBounds
                               
        % Get subset
        function subCG = getSubset(this, cropping)
            
            minX = double(max(this.minX,double(cropping(1))));
            minY = double(max(this.minY,double(cropping(3))));
            maxX = double(min(this.maxX,double(cropping(2))));
            maxY = double(min(this.maxY,double(cropping(4))));
            
            if(minX>=maxX || minY>=maxY)
                error('cropping {%4.2f, %4.2f, %4.2f, %4.2f} is incorrect or outside CT dimensions.\nCheck values of cropping.\n', ...
                    minX, maxX, minY, maxY);
            end
            
            % ROUND
            startI = 1 + round((minX-this.minX)/this.dx);
            stopI = 1 + round((maxX-this.minX)/this.dx);
            startJ = 1 + round((minY-this.minY)/this.dy);
            stopJ = 1 + round((maxY-this.minY)/this.dy);
            
            % Finalize
            NX = 1 + stopI - startI;
            NY = 1 + stopJ - startJ;
            maxX = this.minX + double(stopI - 1)*this.dx;
            minX = this.minX + double(startI - 1)*this.dx;
            maxY = this.minY + double(stopJ - 1)*this.dy;
            minY = this.minY + double(startJ - 1)*this.dy;
            
            % Crop matrix
            subCG = CartesianGrid2D(minX, maxX, minY, maxY, ...
                this.dx, this.dy, NX, NY, ...
                this.data(startI:stopI,startJ:stopJ));            
        end % getSubset
        
        % Obtain indexes
        function [I,J] = obtainIndexes(this,inX,inY)
            validateattributes(inX, {'numeric'}, {'nonempty'}, 'CartesianGrid2D.obtainIndexes', 'inX', 1);
            validateattributes(inY, {'numeric'}, {'nonempty'}, 'CartesianGrid2D.obtainIndexes', 'inY', 2);
            
            % To use with only one vectorial argument
            if nargin==2
                inY = inX(:,2);
                inX = inX(:,1);
            end
            
            I = 1 + (inX - this.minX)./this.dx;
            J = 1 + (inY - this.minY)./this.dy;
            
            I(round(I)<1 | round(I)>this.NX) = NaN;
            J(round(J)<1 | round(J)>this.NY) = NaN;
            
        end % obtainIndexes
        
        % Get Data at certain point (not efficient, do not use inside calculation loops)        
        function data = getData(this, x, y)
            [I,J] = this.obtainIndexes(x,y);
            data = this.data(sub2ind(size(this.data),round(I),round(J)));  
        end % getData
        
        % Get Axis values for plotting or other situations
        function axisValues = getAxisValues(this, axis)
            if strcmpi(axis, 'X')
                %axisValues = (this.minX-this.dx/2):this.dx:(this.maxX+this.dx/2);
                axisValues = (this.minX):this.dx:(this.maxX);                
            elseif strcmpi(axis, 'Y')
                %axisValues = (this.minY-this.dy/2):this.dy:(this.maxY+this.dy/2);                        
                axisValues = (this.minY):this.dy:(this.maxY);                                              
            else
                error('Incorrect axis. Use only uppercase X or Y');
            end                
        end % getAxisValues
        
        % IJK / XYZ conversion
        function I = getXindex(this, X)
            I = 1 + round((X - this.minX)./this.dx);
        end
        function J = getYindex(this, Y)
            J = 1 + round((Y - this.minY)./this.dy);
        end

        function X = getXvalue(this, I)
            X = this.minX + (I-1)*this.dx;
        end
        function Y = getYvalue(this, J)
            Y = this.minY + (J-1)*this.dy;
        end
           
        % Plot a 1D linear plot
        function plotLinear(this, freeAxis, value, format)
%            if nargin>5
%                validateattributes(format, {''}, {}, 'ResultGrids.plotLinear', 'format', 5);              
%            else
%                format='b';
%            end
            % Obtain indexes
            [linData, positions] = this.getLinear(freeAxis, value);
            if nargin==3 || isempty(format) 
                plot(positions,linData);
            else
                plot(positions,linData, format);
            end
                        
        end % plotLinear
        
        % Extract a 1D linear plot                
        function [linData, positions] = getLinear(this, freeAxis, value)
            validateattributes(freeAxis, {'char'}, {'nonempty'}, 'CartesianDose.getLinear', 'freeAxis', 1);
            validateattributes(value, {'numeric'}, {'nonempty'}, 'CartesianDose.getLinear', 'value', 2);                       
                                
            % Select data
            if strcmp(freeAxis,'X')
                index_X = this.getXindex(value);
                if index_X <1 || index_X > this.NX
                    error('The specified position is out of bounds. Cannot produce the plot');
                end
                positions = this.minY:this.dy:this.maxY;
                linData = squeeze(this.data(index_X,:));
            elseif strcmp(freeAxis,'Y')          
                index_Y = this.getYindex(value);                
                if index_Y <1 || index_Y > this.NY
                    error('The specified position is out of bounds. Cannot produce the plot');
                end                
                positions = this.minX:this.dx:this.maxX;
                linData = squeeze(this.data(:,index_Y));
            else
                error('Incorrect axis. Use uppercase X or Y.');
            end
            
        end % getLinear        
        
        % Increase resolution of calculated grid (for example to allow exporting to DICOM)
        function increaseResolution(this,axis,factor)
            validateattributes(axis, {'char'}, {'nonempty'}, 'CartesianGrid2D.increaseResolution', 'axis', 1);
            validateattributes(factor, {'numeric'}, {'nonnan', 'real', 'positive'}, ...
                'CartesianGrid2D.increaseResolution', 'factor', 2);     
            
            factor = double(factor);
            
            if strcmp(axis,'X')

                this.dx = this.dx./factor;
                this.NX = this.NX.*factor;
                this.data = this.data(ceil((1:(factor*size(this.data,2)))/factor), :);
                fprintf('Result matrix resized in X direction. New dx is %3.3f\n',this.dx);
                
            elseif strcmp(axis,'Y')
                
                this.dy = this.dy./factor;
                this.NY = this.NY.*factor;
                this.data = this.data(:, ceil((1:(factor*size(this.data,2)))/factor));
                fprintf('Result matrix resized in Y direction. New dy is %3.3f\n',this.dy);               
                                
            else
                error('Incorrect axis. Use uppercase X, Y or Z.');
            end            
            
        end % increaseResolution        
                
    end % methods
    
end % classdef

