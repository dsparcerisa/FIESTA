% CartesianGrid class.
% Auxiliary class to store cartesian objects with spacial coordinates
%
% Constructors:
%  CartesianGrid(CartesianGrid CGobject)
% --> Copies structure from another CGrid object with empty data
%
%  CartesianGrid(minX,maxX,minY,maxY,minZ,maxZ,dx,dy,dz,NX,NY,NZ,data*)
% --> data matrix is optional 
% 
% Properties:
%
%  (Bin sizes in cm)
%   dx
%   dy
%   dz
%
%  (Number of bins) 
%   NX
%   NY
%   NZ
%
%  (Extreme positions in cm at voxel center)
%   minX        
%   maxX
%   minY
%   maxY
%   minZ
%   maxZ
%   data
%
% Methods:
%  rebin(newdx,newdy,newdz)--> has same size but different dx dy dz.
%  resize(newSize) --> has size newSize ([NX NY NZ]) and original dx dy dz.
%  crop(newSize) --> Alias for resize
%  plotSlice(axis, value) --> as in plotSlice('Y',0)
%  passMinMaxTests() --> Verifies dimensions
%  displayCoords() --> Prints coordinates
%  res = posInBounds(x, y, z) --> Returns TRUE if position is in bounds
%  res = indexInBounds(i, j, k) --> Returns TRUE if indexes IJK are in bounds
%  slice = getSlice(axis, value) --> as in getSlice('Y',0)
%  subCG = getSubset(cropping) --> Returns a cropped CG object
%  [I,J,K] = obtainIndexes(inX,inY,inZ) --> maps XYZ to indexes
%  data = getData(x, y, z) --> Indexing by position
%  axisValues = getAxisValues(axis) --> Returns vector of axis positions
%  I = getXindex(X)
%  J = getYindex(Y)
%  K = getZindex(Z)
%  X = getXvalue(I)
%  Y = getYvalue(J)
%  Z = getZvalue(K)
%  plotLinear(freeAxis, X, Y, Z, format) --> Create linear plot
%  [linData, positions] = getLinear(freeAxis, X, Y, Z) --> Extract linear curve
%  increaseResolution(axis, factor) --> Increase resolution in given axis
%
% Overloaded operators:
%  plus/minus --> A + b has max(A,b) size and min(A,b) in dx, dy, dz
%                (individually for each dimension)
%                 A - b = A + (-B);
%  unary minus --> -CG
%  indexing --> CG(i,j,k) => CG.data(i,j,k);
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

classdef CartesianGrid < matlab.mixin.Copyable
    %CARTESIANGRID Auxiliary class to store cartesian objects with spacial coordinates
    
    properties
        % All distances are in cm
        % min/max values represent the CENTER of the voxels

        % Bin sizes (cm)
        dx
        dy
        dz
        
        % Number of bins
        NX
        NY
        NZ
        
        % Extremal positions (voxel centers)
        minX        
        maxX
        minY
        maxY
        minZ
        maxZ
        
        % Data matrix
        data
        
    end % properties
    
    methods        
        
        % Constructor
        function this = CartesianGrid(varargin)
            if nargin == 1
            % CartesianGrid(CGobject)
                    validateattributes(varargin{1}, {'CartesianGrid'}, {'nonempty'}, 'CartesianGrid Constructor', 'CGObject', 1);
                    CGobject = varargin{1};
                    this.dx = CGobject.dx;
                    this.dy = CGobject.dy;
                    this.dz = CGobject.dz;
                    this.NX = CGobject.NX;
                    this.NY = CGobject.NY;
                    this.NZ = CGobject.NZ;
                    this.minX = CGobject.minX;
                    this.maxX = CGobject.maxX;
                    this.minY = CGobject.minY;
                    this.maxY = CGobject.maxY;
                    this.minZ = CGobject.minZ;
                    this.maxZ = CGobject.maxZ;
                    this.passMinMaxTests;                    
                    this.data = zeros(this.NX, this.NY, this.NZ);  
                    
            elseif nargin == 12 || nargin == 13
            % CartesianGrid(minX,maxX,minY,maxY,minZ,maxZ,dx,dy,dz,NX,NY,NZ,data*)
                
                    validateattributes(varargin{1}, {'numeric'}, {'nonnan', 'real'}, 'CartesianGrid Constructor', 'minX', 1);
                    this.minX = varargin{1};
                    validateattributes(varargin{2}, {'numeric'}, {'nonnan', 'real'}, 'CartesianGrid Constructor', 'maxX', 2);
                    this.maxX = varargin{2};
                    validateattributes(varargin{3}, {'numeric'}, {'nonnan', 'real'}, 'CartesianGrid Constructor', 'minY', 3);
                    this.minY = varargin{3};
                    validateattributes(varargin{4}, {'numeric'}, {'nonnan', 'real'}, 'CartesianGrid Constructor', 'maxY', 4);
                    this.maxY = varargin{4};
                    validateattributes(varargin{5}, {'numeric'}, {'nonnan', 'real'}, 'CartesianGrid Constructor', 'minZ', 5);
                    this.minZ = varargin{5};
                    validateattributes(varargin{6}, {'numeric'}, {'nonnan', 'real'}, 'CartesianGrid Constructor', 'maxZ', 6);
                    this.maxZ = varargin{6};
                    validateattributes(varargin{7}, {'numeric'}, {'nonnan', 'nonnegative', 'real'}, 'CartesianGrid Constructor', 'dx', 7);
                    this.dx = varargin{7};
                    validateattributes(varargin{8}, {'numeric'}, {'nonnan', 'nonnegative','real'}, 'CartesianGrid Constructor', 'dy', 8);
                    this.dy = varargin{8};
                    validateattributes(varargin{9}, {'numeric'}, {'nonnan', 'nonnegative','real'}, 'CartesianGrid Constructor', 'dz', 9);
                    this.dz = varargin{9};
                    validateattributes(varargin{10}, {'numeric'}, {'nonnan', 'positive', 'integer'}, 'CartesianGrid Constructor', 'NX', 10);
                    this.NX = varargin{10};
                    validateattributes(varargin{11}, {'numeric'}, {'nonnan', 'positive', 'integer'}, 'CartesianGrid Constructor', 'NY', 11);
                    this.NY = varargin{11};
                    validateattributes(varargin{12}, {'numeric'}, {'nonnan', 'positive', 'integer'}, 'CartesianGrid Constructor', 'NZ', 12);
                    this.NZ = varargin{12};
                                        
                    this.passMinMaxTests;
                    
                    if nargin==13
                        validateattributes(varargin{13}, {'numeric'}, {'size',[this.NX this.NY this.NZ]}, 'CartesianGrid Constructor', 'data', 13);
                        this.data = varargin{13};
                    else
                        this.data = zeros(this.NX, this.NY, this.NZ);
                    end
                
            elseif nargin>0
                error('Incorrect number of arguments');
            end
            
        end % constructor
                        
        % rebin(newdx, newdy, newdz) has same size but different dx dy dz
        
        function rebin(this,newdx,newdy,newdz)
            validateattributes(newdx, {'numeric'}, {'positive','real'}, 'CartesianGrid.rebin', 'newdx', 1);             
            validateattributes(newdy, {'numeric'}, {'positive','real'}, 'CartesianGrid.rebin', 'newdy', 2);             
            validateattributes(newdz, {'numeric'}, {'positive','real'}, 'CartesianGrid.rebin', 'newdz', 3);             
            
            Xpositions = this.minX:newdx:this.maxX;
            Ypositions = this.minY:newdy:this.maxY;
            Zpositions = this.minZ:newdz:this.maxZ;
            dataType = class(this.data);
            GI = griddedInterpolant({this.minX:this.dx:this.maxX, ...
                this.minY:this.dy:this.maxY, ...
                this.minZ:this.dz:this.maxZ},double(this.data),'linear','none');
            this.data = cast(this.data,dataType);
            this.data = GI({Xpositions, Ypositions, Zpositions});       
            this.NX = numel(Xpositions);
            this.NY = numel(Ypositions);
            this.NZ = numel(Zpositions);
            this.dx = newdx;
            this.dy = newdy;
            this.dz = newdz;
            this.minX = Xpositions(1);
            this.minY = Ypositions(1);
            this.minZ = Zpositions(1);
            this.maxX = Xpositions(end);
            this.maxY = Ypositions(end);
            this.maxZ = Zpositions(end);    
            
            this.passMinMaxTests;                        
        end % rebin
                
        % resize(newSize) has size newSize ([NX NY NZ]) and original dx dy dz.
        
        function resize(this,newSize)
            validateattributes(newSize, {'numeric'}, {'nonnan','real','numel',6}, 'CartesianGrid.crop', 'newSize', 1);             
            Xpositions = newSize(1):this.dx:newSize(2);
            Ypositions = newSize(3):this.dy:newSize(4);
            Zpositions = newSize(5):this.dz:newSize(6);
            dataType = class(this.data);
            GI = griddedInterpolant({this.minX:this.dx:this.maxX, ...
                this.minY:this.dy:this.maxY, ...
                this.minZ:this.dz:this.maxZ},double(this.data),'linear','none');
            this.data = cast(this.data,dataType);
            this.data = GI({Xpositions, Ypositions, Zpositions});       
            this.NX = numel(Xpositions);
            this.NY = numel(Ypositions);
            this.NZ = numel(Zpositions);
            this.minX = Xpositions(1);
            this.minY = Ypositions(1);
            this.minZ = Zpositions(1);
            this.maxX = Xpositions(end);
            this.maxY = Ypositions(end);
            this.maxZ = Zpositions(end);    
            
            this.passMinMaxTests;
        end % resize
        
        % CG: A + b has max(A,b) size and min(A,b) in dx, dy, dz
        % (individually for each dimension)
        % A - b = A + (-B);
        
        function Result = plus(this, CGRight)
        % CG: A + b has max(A,b) size and min(A,b) in dx, dy, dz
        % (individually for each dimension)
        % A - b = A + (-B);            
            validateattributes(this, {'CartesianGrid'}, {'nonempty'}, 'CartesianGrid.plus', 'CGLeft', 1);                         
            validateattributes(CGRight, {'CartesianGrid'}, {'nonempty'}, 'CartesianGrid.plus', 'CGRight', 2);     

            minX = min(this.minX, CGRight.minX); 
            minY = min(this.minX, CGRight.minX); 
            minZ = min(this.minX, CGRight.minX); 
            maxX = max(this.maxX, CGRight.maxX); 
            maxY = max(this.maxY, CGRight.maxY); 
            maxZ = max(this.maxZ, CGRight.maxZ); 
            dx = min(this.dx, CGRight.dx);
            dy = min(this.dy, CGRight.dy);
            dz = min(this.dz, CGRight.dz);
            Xpositions = minX:dx:maxX;
            Ypositions = minY:dy:maxY;
            Zpositions = minZ:dz:maxZ;
            maxX = Xpositions(end);
            maxY = Ypositions(end);
            maxZ = Zpositions(end);
            
            GI_Left = griddedInterpolant({this.minX:this.dx:this.maxX, ...
                this.minY:this.dy:this.maxY, ...
                this.minZ:this.dz:this.maxZ},this.data,'linear','none');
            GI_Right = griddedInterpolant({CGRight.minX:CGRight.dx:CGRight.maxX, ...
                CGRight.minY:CGRight.dy:CGRight.maxY, ...
                CGRight.minZ:CGRight.dz:CGRight.maxZ},CGRight.data,'linear','none');
            
            NX = numel(Xpositions);
            NY = numel(Ypositions);
            NZ = numel(Zpositions);            
            data = GI_Left({Xpositions,Ypositions,Zpositions}) + GI_Right({Xpositions,Ypositions,Zpositions});
            
            Result = CartesianGrid(minX,maxX,minY,maxY,minZ,maxZ,dx,dy,dz,NX,NY,NZ,data);                 
            
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
        function plotSlice(this, axis, value)
            validateattributes(axis, {'char'}, {'nonempty'}, 'CartesianGrid.plotSlice', 'axis', 1);
            validateattributes(value, {'numeric'}, {'nonnan', 'real'}, 'CartesianGrid.plotSlice', 'value', 2);                
            
            % Get slice
            slice = this.getSlice(axis,value);
                   
            % Plot
            if strcmp(axis,'X')
                imagesc(this.getAxisValues('Z'),fliplr(this.getAxisValues('Y')),flipud(slice));
                
                set(gca,'xlim',[this.minZ this.maxZ]);
                set(gca,'ylim',[this.minY this.maxY]);
                
                xlabel('Z');
                ylabel('Y');
                
            elseif strcmp(axis,'Y')
                imagesc(this.getAxisValues('X'),fliplr(this.getAxisValues('Z')),rot90(slice));  
              
                set(gca,'xlim',[this.minX this.maxX]);
                set(gca,'ylim',[this.minZ this.maxZ]);                
                
                xlabel('X');
                ylabel('Z');
                
            elseif strcmp(axis,'Z')
                imagesc(this.getAxisValues('X'),fliplr(this.getAxisValues('Y')),rot90(slice));                
             
                set(gca,'xlim',[this.minX this.maxX]);
                set(gca,'ylim',[this.minY this.maxY]);                    
                
                xlabel('X');
                ylabel('Y');                
            else
                error('Incorrect axis. Use uppercase X, Y or Z.');
            end     
            
            set(gca,'Ydir','normal');
            set(gca,'Xdir','normal');
            
        end % plotSlice        
        
        % Test minima and maxima        
        function passMinMaxTests(this)
            if ( abs(this.maxX - (this.minX + this.dx*double(this.NX-1))) > 1e-6 ) || ...
                    ( abs(this.maxY - (this.minY + this.dy*double(this.NY-1))) > 1e-6 ) || ...
                    ( abs(this.maxZ - (this.minZ + this.dz*double(this.NZ-1))) > 1e-6 )
                error('Minima/maxima tests not passed.\n\tDEBUG INFO:\n\tTest X: %f\n\tTest Y: %f\n\tTest Z: %f', ...
                    ( abs(this.maxX - (this.minX + this.dx*double(this.NX-1))) ), ...
                    ( abs(this.maxY - (this.minY + this.dy*double(this.NY-1))) ), ...
                    ( abs(this.maxZ - (this.minZ + this.dz*double(this.NZ-1))) ) )                            
            end
        end % passMinMaxTests
        
        % Display coordinates
        function displayCoords(this)
            fprintf('## X from %g to %g \n', this.minX, this.maxX);
            fprintf('## Y from %g to %g \n', this.minY, this.maxY);
            fprintf('## Z from %g to %g \n\n', this.minZ, this.maxZ);
        end % displayCoords
        
        % Position in bounds
        function res = posInBounds(this, x, y, z)
            res = x>(this.minX - this.dx/2) & ...
                  x<(this.maxX + this.dx/2) & ...
                  y>(this.minY - this.dy/2) & ...
                  y<(this.maxY + this.dy/2) & ...
                  z>(this.minZ - this.dz/2) & ...
                  z<(this.maxZ + this.dz/2);              
        end % posInBounds
        
        % Index in bounds
        function res = indexInBounds(this, i, j, k)
            res = i>=1 & j>=1 & k>=1 & i<=this.NX & j<=this.NY & k<=this.NZ;
        end % indexInBounds
        
        % Get dose slice
        function slice = getSlice(this, axis, value)
            validateattributes(axis, {'char'}, {'nonempty'}, 'CartesianGrid.getSlice', 'axis', 1);
            validateattributes(value, {'numeric'}, {'nonnan', 'real'}, 'CartesianGrid.getSlice', 'value', 2);
            
            if strcmp(axis,'X')
                validateattributes(value, {'numeric'}, {'>=', this.minX-this.dx/2, '<=', this.maxX+this.dx/2}, 'CartesianGrid.getSlice', 'X value', 2);                
                X = value;                
                index_X=1+round((X-this.minX)/(this.dx));
                fprintf('Selecting dose slice at x=%f with index=%i\n',this.minX+(index_X-1)*this.dx, index_X);                
                slice=squeeze(this.data(index_X,1:this.NY,1:this.NZ));
            elseif strcmp(axis,'Y')
                validateattributes(value, {'numeric'}, {'>=', this.minY-this.dy/2, '<=', this.maxY+this.dy/2}, 'CartesianGrid.getSlice', 'Y value', 2);                                
                Y = value;
                index_Y=1+round((Y-this.minY)/(this.dy));
                fprintf('Selecting dose slice at y=%f with index=%i\n',this.minY+(index_Y-1)*this.dy, index_Y);
                slice=squeeze(this.data(1:this.NX,index_Y,1:this.NZ));
            elseif strcmp(axis,'Z')
                validateattributes(value, {'numeric'}, {'>=', this.minZ-this.dz/2, '<=', this.maxZ+this.dz/2}, 'CartesianGrid.getSlice', 'Z value', 2);                                
                Z = value;
                index_Z=1+round((Z-this.minZ)/(this.dz));
                fprintf('Selecting dose slice at z=%f with index=%i\n',this.minZ+(index_Z-1)*this.dz, index_Z);
                slice=squeeze(this.data(1:this.NX,1:this.NY,index_Z));
            else
                error('Incorrect axis. Use uppercase X, Y or Z.');
            end
        end % getSlice
                       
        % Get subset
        function subCG = getSubset(this, cropping)
            
            minX = double(max(this.minX,double(cropping(1))));
            minY = double(max(this.minY,double(cropping(3))));
            minZ = double(max(this.minZ,double(cropping(5))));
            maxX = double(min(this.maxX,double(cropping(2))));
            maxY = double(min(this.maxY,double(cropping(4))));
            maxZ = double(min(this.maxZ,double(cropping(6))));
            
            if(minX>=maxX || minY>=maxY || minZ>=maxZ)
                error('cropping {%4.2f, %4.2f, %4.2f, %4.2f, %4.2f, %4.2f} is incorrect or outside CT dimensions.\nCheck values of cropping and userOrigin ([%4.2f, %4.2f, %4.2f])', ...
                    minX, maxX, minY, maxY, minZ, maxZ, this.userOrigin(1), this.userOrigin(2), this.userOrigin(3));
            end
            
            % ROUND
            startI = 1 + round((minX-this.minX)/this.dx);
            stopI = 1 + round((maxX-this.minX)/this.dx);
            startJ = 1 + round((minY-this.minY)/this.dy);
            stopJ = 1 + round((maxY-this.minY)/this.dy);
            startK = 1 + round((minZ-this.minZ)/this.dz);
            stopK = 1 + round((maxZ-this.minZ)/this.dz);
            
            % Finalize
            NX = 1 + stopI - startI;
            NY = 1 + stopJ - startJ;
            NZ = 1 + stopK - startK;
            maxX = this.minX + double(stopI - 1)*this.dx;
            minX = this.minX + double(startI - 1)*this.dx;
            maxY = this.minY + double(stopJ - 1)*this.dy;
            minY = this.minY + double(startJ - 1)*this.dy;
            maxZ = this.minZ + double(stopK - 1)*this.dz;
            minZ = this.minZ + double(startK - 1)*this.dz;
            
            % Crop matrix
            subCG = CartesianGrid(minX, maxX, minY, maxY, minZ, maxZ, ...
                this.dx, this.dy, this.dz, NX, NY, NZ, ...
                this.data(startI:stopI,startJ:stopJ,startK:stopK));            
        end % getSubset
        
        % Obtain indexes
        function [I,J,K] = obtainIndexes(this,inX,inY,inZ)
            validateattributes(inX, {'numeric'}, {'nonempty'}, 'CartesianGrid.obtainIndexes', 'inX', 1);
            validateattributes(inY, {'numeric'}, {'nonempty'}, 'CartesianGrid.obtainIndexes', 'inY', 2);
            validateattributes(inZ, {'numeric'}, {'nonempty'}, 'CartesianGrid.obtainIndexes', 'inZ', 3);
            
            % To use with only one vectorial argument
            if nargin==2
                inZ = inX(:,3);
                inY = inX(:,2);
                inX = inX(:,1);
            end
            
            I = 1 + (inX - this.minX)./this.dx;
            J = 1 + (inY - this.minY)./this.dy;
            K = 1 + (inZ - this.minZ)./this.dz;
            
            I(round(I)<1 | round(I)>this.NX) = NaN;
            J(round(J)<1 | round(J)>this.NY) = NaN;
            K(round(K)<1 | round(K)>this.NZ) = NaN;                                  
            
        end % obtainIndexes
        
        % Get Data at certain point (not efficient, do not use inside calculation loops)        
        function data = getData(this, x, y, z)
            [I,J,K] = this.obtainIndexes(x,y,z);
            data = this.data(sub2ind(size(this.data),round(I),round(J),round(K)));  
        end % getData
        
        % Get Axis values for plotting or other situations
        function axisValues = getAxisValues(this, axis)
            if strcmp(axis, 'X')
                axisValues = this.minX:this.dx:this.maxX;
            elseif strcmp(axis, 'Y')
                axisValues = this.minY:this.dy:this.maxY;                
            elseif strcmp(axis, 'Z')
                axisValues = this.minZ:this.dz:this.maxZ;                
            else
                error('Incorrect axis. Use only uppercase X, Y or Z');
            end                
        end % getAxisValues
        
        % IJK / XYZ conversion
        function I = getXindex(this, X)
            I = 1 + round((X - this.minX)./this.dx);
        end
        function J = getYindex(this, Y)
            J = 1 + round((Y - this.minY)./this.dy);
        end
        function K = getZindex(this, Z)
            K = 1 + round((Z - this.minZ)./this.dz);
        end
        function X = getXvalue(this, I)
            X = this.minX + (I-1)*this.dx;
        end
        function Y = getYvalue(this, J)
            Y = this.minY + (J-1)*this.dy;
        end
        function Z = getZvalue(this, K)
            Z = this.minZ + (K-1)*this.dz;
        end
           
        % Plot a 1D linear plot
        function plotLinear(this, freeAxis, X, Y, Z, format)
%            if nargin>5
%                validateattributes(format, {''}, {}, 'ResultGrids.plotLinear', 'format', 5);              
%            else
%                format='b';
%            end
            % Obtain indexes
            [linData, positions] = this.getLinear(freeAxis, X, Y, Z);
            if isempty(format) || nargin==5
                plot(positions,linData)
            else
                plot(positions,linData, format);
            end
                        
        end % plotLinear
        
        % Extract a 1D linear plot                
        function [linData, positions] = getLinear(this, freeAxis, X, Y, Z)
            validateattributes(freeAxis, {'char'}, {'nonempty'}, 'CartesianDose.getLinear', 'freeAxis', 1);
            validateattributes(X, {'numeric'}, {'nonempty'}, 'CartesianDose.getLinear', 'X', 2);            
            validateattributes(Y, {'numeric'}, {'nonempty'}, 'CartesianDose.getLinear', 'Y', 3);            
            validateattributes(Z, {'numeric'}, {'nonempty'}, 'CartesianDose.getLinear', 'Z', 4);               
            
            if ~this.posInBounds(X,Y,Z)
                error('The specified position is out of bounds. Cannot produce the plot');
            end            
            
            % Obtain indexes
            index_X = this.getXindex(X);
            index_Y = this.getYindex(Y);
            index_Z = this.getZindex(Z);
            
            % Plot
            if strcmp(freeAxis,'X')
                positions = this.minX:this.dx:this.maxX;
                linData = squeeze(this.data(:,index_Y,index_Z));
            elseif strcmp(freeAxis,'Y')
                positions = this.minY:this.dy:this.maxY;
                linData = squeeze(this.data(index_X,:,index_Z));
            elseif strcmp(freeAxis,'Z')
                positions = this.minZ:this.dz:this.maxZ;
                linData = squeeze(this.data(index_X,index_Y,:));
            else
                error('Incorrect axis. Use uppercase X, Y or Z.');
            end
            
        end % getLinear        
        
        % Increase resolution of calculated grid (for example to allow exporting to DICOM)
        function increaseResolution(this,axis,factor)
            validateattributes(axis, {'char'}, {'nonempty'}, 'CartesianGrid.increaseResolution', 'axis', 1);
            validateattributes(factor, {'numeric'}, {'nonnan', 'real', 'positive'}, ...
                'CartesianGrid.increaseResolution', 'factor', 2);     
            
            factor = double(factor);
            
            if strcmp(axis,'X')

                this.dx = this.dx./factor;
                this.NX = this.NX.*factor;
                this.data = this.data(ceil((1:(factor*size(this.data,2)))/factor), :, :);
                fprintf('Result matrix resized in X direction. New dx is %3.3f\n',this.dx);
                
            elseif strcmp(axis,'Y')
                
                this.dy = this.dy./factor;
                this.NY = this.NY.*factor;
                this.data = this.data(:, ceil((1:(factor*size(this.data,2)))/factor), :);
                fprintf('Result matrix resized in Y direction. New dy is %3.3f\n',this.dy);

                
            elseif strcmp(axis,'Z')
                
                this.dz = this.dz./factor;
                this.NZ = this.NZ.*factor;
                this.data = this.data(:, :, ceil((1:(factor*size(this.data,2)))/factor));
                fprintf('Result matrix resized in Z direction. New dz is %3.3f\n',this.dz);
                                
            else
                error('Incorrect axis. Use uppercase X, Y or Z.');
            end            
            
        end % increaseResolution        
                
    end % methods
    
end % classdef

