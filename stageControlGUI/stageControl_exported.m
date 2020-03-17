classdef stageControl_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                  matlab.ui.Figure
        LimitsPanel               matlab.ui.container.Panel
        XMinEditFieldLabel        matlab.ui.control.Label
        XMinEditField             matlab.ui.control.NumericEditField
        XMaxEditFieldLabel        matlab.ui.control.Label
        XMaxEditField             matlab.ui.control.NumericEditField
        YMinEditFieldLabel        matlab.ui.control.Label
        YMinEditField             matlab.ui.control.NumericEditField
        YMaxEditFieldLabel        matlab.ui.control.Label
        YMaxEditField             matlab.ui.control.NumericEditField
        ZMinEditFieldLabel        matlab.ui.control.Label
        ZMinEditField             matlab.ui.control.NumericEditField
        ZMaxEditFieldLabel        matlab.ui.control.Label
        ZMaxEditField             matlab.ui.control.NumericEditField
        ClearButton               matlab.ui.control.Button
        SetLimitsButton           matlab.ui.control.Button
        AutoCenterButton          matlab.ui.control.Button
        GotoPanel                 matlab.ui.container.Panel
        XEditFieldLabel           matlab.ui.control.Label
        XEditField                matlab.ui.control.NumericEditField
        YEditFieldLabel           matlab.ui.control.Label
        YEditField                matlab.ui.control.NumericEditField
        ZEditFieldLabel           matlab.ui.control.Label
        ZEditField                matlab.ui.control.NumericEditField
        GoButton                  matlab.ui.control.Button
        StepmovementPanel         matlab.ui.container.Panel
        Image                     matlab.ui.control.Image
        XButton_minus             matlab.ui.control.Button
        ZButton_plus              matlab.ui.control.Button
        YButton_plus              matlab.ui.control.Button
        YButton_minus             matlab.ui.control.Button
        ZButton_minus             matlab.ui.control.Button
        XButton_plus              matlab.ui.control.Button
        StepsizecmEditFieldLabel  matlab.ui.control.Label
        StepsizecmEditField       matlab.ui.control.NumericEditField
        RefreshButton             matlab.ui.control.Button
        StopButton                matlab.ui.control.Button
        CurrentpositionPanel      matlab.ui.container.Panel
        XLabel                    matlab.ui.control.Label
        currentXvalue             matlab.ui.control.NumericEditField
        minEditField_2Label       matlab.ui.control.Label
        minXvalue                 matlab.ui.control.NumericEditField
        maxEditField_4Label       matlab.ui.control.Label
        maxXvalue                 matlab.ui.control.NumericEditField
        YLabel                    matlab.ui.control.Label
        currentYvalue             matlab.ui.control.NumericEditField
        EditFieldLabel_2          matlab.ui.control.Label
        maxYvalue                 matlab.ui.control.NumericEditField
        minEditField_3Label       matlab.ui.control.Label
        minYvalue                 matlab.ui.control.NumericEditField
        ZLabel                    matlab.ui.control.Label
        currentZvalue             matlab.ui.control.NumericEditField
        EditFieldLabel_3          matlab.ui.control.Label
        maxZvalue                 matlab.ui.control.NumericEditField
        minEditField_4Label       matlab.ui.control.Label
        minZvalue                 matlab.ui.control.NumericEditField
        ZerohereButton            matlab.ui.control.Button
        Image2                    matlab.ui.control.Image
        Image3                    matlab.ui.control.Image
        DSanchezSEspaaLabel       matlab.ui.control.Label
        LogTextAreaLabel          matlab.ui.control.Label
        LogTextArea               matlab.ui.control.TextArea
        StagestatusPanel          matlab.ui.container.Panel
        ConnectedLampLabel        matlab.ui.control.Label
        ConnectedLamp             matlab.ui.control.Lamp
        SaveButton                matlab.ui.control.Button
    end

    
    properties (Access = private)
        currentX % current X position in cm
        currentY % current Y position in cm
        currentZ % current Z position in cm
        Xmin
        Xmax
        Ymin
        Ymax
        Zmin
        Zmax
        notInitialized  % Description
        status      % 0: Off, 1: Connected, 2: Moving
        velmexPort  % Serial object connecting to stage controller
    end
    
    methods (Static)
        
        function singleObj = getInstance(COM)
            persistent localObj
            if isempty(localObj) || ~isvalid(localObj)
                localObj = stageControl_exported(COM);
            end
            singleObj = localObj;
        end
        
    end
    
    methods (Access = private)
        
        function updateValues(app)
            
            app.XMaxEditField.Value = app.Xmax;
            app.XMinEditField.Value = app.Xmin;
            app.YMaxEditField.Value = app.Ymax;
            app.YMinEditField.Value = app.Ymin;
            app.ZMaxEditField.Value = app.Zmax;
            app.ZMinEditField.Value = app.Zmin;
            
            app.maxXvalue.Value = app.Xmax;
            app.minXvalue.Value = app.Xmin;
            app.maxYvalue.Value = app.Ymax;
            app.minYvalue.Value = app.Ymin;
            app.maxZvalue.Value = app.Zmax;
            app.minZvalue.Value = app.Zmin;
            
            app.currentXvalue.Value = app.currentX;
            app.currentYvalue.Value = app.currentY;
            app.currentZvalue.Value = app.currentZ;
            
            
        end        
        
       
        function moveTo(app, posVector)  
            try
            app.autoMonitorStatus
            if app.status==1
                app.setStatus(2);
                stageControl_setPosition(app, posVector);
                while (true)
                    app.autoMonitorStatus;
                    if (app.status==1)
                        break                      
                    end
                end
                app.setStatus(1);
            elseif app.status==2
                app.logLine('WARNING: Wait until movement finishes...');
            elseif app.status==0
                app.logLine('ERROR: Stage is offline');
            end
            app.autoMonitorStatus
            message = sprintf('Stage is now at (X: %3.2f  Y: %3.2f  Z: %3.2f)\n', app.currentX, app.currentY, app.currentZ);
            app.logLine(message);    
            catch
                app.logLine('WARNING: Movement stopped unexpectedly')
                app.autoMonitorStatus;
                message = sprintf('Stage is now at (X: %3.2f  Y: %3.2f  Z: %3.2f)\n', app.currentX, app.currentY, app.currentZ);
                app.logLine(message);                
            end
        end
        
        function autoMonitorStatus(app)
            [readStatus, posX, posY, posZ] = monitorStatus(app.velmexPort);
            app.setStatus(readStatus);
            switch readStatus
                case 1
                    app.updatePosition([posX posY posZ]);
                case 2
                    app.updatePosition([posX posY app.currentZ]);
            end
        end
    end
    
    methods (Access = public)
        
        function result = isWithinLimits(app, position)
            result = true;
            if position(1) < app.Xmin
                result = false;
                return
            end
            if position(1) > app.Xmax
                result = false;
                return
            end
            if position(2) < app.Ymin
                result = false;
                return
            end
            if position(2) > app.Ymax
                result = false;
                return
            end
            if position(3) < app.Zmin
                result = false;
                return
            end
            if position(3) > app.Zmax
                result = false;
                return
            end
        end
        
        function updatePosition(app, position)
            app.currentX = position(1);
            app.currentY = position(2);
            app.currentZ = position(3);
            app.updateValues;            
        end
        
        function updateLimits(app, limitVector)
            app.Xmin = limitVector(1);
            app.Xmax = limitVector(2);
            app.Ymin = limitVector(3);
            app.Ymax = limitVector(4);
            app.Zmin = limitVector(5);
            app.Zmax = limitVector(6);
            app.updateValues
        end
        
        function logLine(app, newLine)
            newLine = [datestr(now,'[HH:MM:SS]') ' ' newLine];
            app.LogTextArea.Value = { newLine, app.LogTextArea.Value{:} };
        end
        
        function setStatus(app, status)
            
            switch (status)
                case 0
                    app.status = 0;
                    app.ConnectedLampLabel = 'Offline';
                    app.ConnectedLamp.Color = [1 0 0];
                case 1
                    app.status = 1;
                    app.ConnectedLampLabel.Text = 'Ready';
                    app.ConnectedLamp.Color = [0 1 0];
                case 2
                    app.status = 2;
                    app.ConnectedLampLabel.Text = 'Moving';
                    app.ConnectedLamp.Color = [1 1 0];
                    
                otherwise
                    error('ERROR. Status %i not defined', status);
            end
            app.status = status;
        end
        
        function currentPos = getCurrentPos(app)
            currentPos = [app.currentX app.currentY app.currentZ];
        end
        
        function port = getPort(app)
            port = app.velmexPort;
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, port)
            if isempty(app.notInitialized)
                app.notInitialized = false;
                app.velmexPort = port;
                app.logLine('Initializing...');
                app.currentX = 0;
                app.currentY = 0;
                app.currentZ = 0;
                app.updateLimits([-5 5 -5 5 -5 5]);
                app.setStatus(1); % 0: Offline; 1: Connected; 2: Moving
                app.updateValues;
                app.StepsizecmEditField.Value = 0.2; % Default
                app.autoMonitorStatus
            end
        end

        % Button pushed function: ClearButton
        function ClearButtonPushed(app, event)
            app.logLine('Clear limits')
            app.Xmax = 0;
            app.Xmin = 0;
            app.Ymax = 0;
            app.Ymin = 0;
            app.Zmax = 0;
            app.Zmin = 0;
            
            app.updateValues;
        end

        % Button pushed function: GoButton
        function GoButtonPushed(app, event)
            posVector = [app.XEditField.Value app.YEditField.Value app.ZEditField.Value];
            if app.isWithinLimits(posVector)
                if posVector(1) == app.currentX && ...
                        posVector(2) == app.currentY && ...
                        posVector(3) == app.currentZ
                    app.logLine('Already in target position. Not moved.');
                else
                    app.moveTo([posVector(1) app.currentY app.currentZ]);
                    while(true)
                        readSt = monitorStatus(app.velmexPort);
                        if readSt==1
                            break
                        end
                    end
                    app.moveTo([posVector(1) posVector(2) app.currentZ]);
                    while(true)
                        readSt = monitorStatus(app.velmexPort);
                        if readSt==1
                            break
                        end
                    end
                    app.moveTo([posVector]);
                    
                end
            else
                app.logLine('ERROR: Cannot move outside of defined limits');
            end
        end

        % Button pushed function: SetLimitsButton
        function SetLimitsButtonPushed(app, event)
            limitVector = [app.XMinEditField.Value, ...
                app.XMaxEditField.Value, ...
                app.YMinEditField.Value, ...
                app.YMaxEditField.Value, ...
                app.ZMinEditField.Value, ...
                app.ZMaxEditField.Value];
            stageControl_setLimits(app, limitVector);
        end

        % Button pushed function: ZerohereButton
        function ZerohereButtonPushed(app, event)
            readStatus(app.velmexPort, 'N');
            app.logLine('Position zeroed without moving stage.')
            app.updatePosition([0 0 0]);
            app.autoMonitorStatus
        end

        % Button pushed function: ZButton_plus
        function ZButton_plusPushed(app, event)
            step = app.StepsizecmEditField.Value;
            newPosition = [app.currentX app.currentY (app.currentZ + step)];
            if app.isWithinLimits(newPosition)
                app.moveTo(newPosition);
            else
                app.logLine('ERROR: Stage already at its limit')
            end
        end

        % Button pushed function: ZButton_minus
        function ZButton_minusPushed(app, event)
            step = app.StepsizecmEditField.Value;
            newPosition = [app.currentX app.currentY (app.currentZ - step)];
            if app.isWithinLimits(newPosition)
                app.moveTo(newPosition);
            else
                app.logLine('ERROR: Stage already at its limit')
            end
        end

        % Button pushed function: YButton_plus
        function YButton_plusPushed(app, event)
            step = app.StepsizecmEditField.Value;
            newPosition = [app.currentX (app.currentY+step) app.currentZ];
            if app.isWithinLimits(newPosition)
                app.moveTo(newPosition);
            else
                app.logLine('ERROR: Stage already at its limit')
            end
        end

        % Button pushed function: YButton_minus
        function YButton_minusPushed(app, event)
            step = app.StepsizecmEditField.Value;
            newPosition = [app.currentX (app.currentY-step) app.currentZ];
            if app.isWithinLimits(newPosition)
                app.moveTo(newPosition);
            else
                app.logLine('ERROR: Stage already at its limit')
            end
        end

        % Button pushed function: XButton_minus
        function XButton_minusPushed(app, event)
            step = app.StepsizecmEditField.Value;
            newPosition = [(app.currentX-step) app.currentY app.currentZ];
            if app.isWithinLimits(newPosition)
                app.moveTo(newPosition);
            else
                app.logLine('ERROR: Stage already at its limit')
            end
        end

        % Button pushed function: XButton_plus
        function XButton_plusPushed(app, event)
            step = app.StepsizecmEditField.Value;
            newPosition = [(app.currentX+step) app.currentY app.currentZ];
            if app.isWithinLimits(newPosition)
                app.moveTo(newPosition);
            else
                app.logLine('ERROR: Stage already at its limit')
            end
        end

        % Button pushed function: SaveButton
        function SaveButtonPushed(app, event)
            logFileName = datestr(now,'dd_mmmm_yyyy_HH_MM_SS.log');
            
            filePh = fopen(logFileName,'w');
            logText = flip(app.LogTextArea.Value);
            fprintf(filePh,'%s\n',logText{:});
            fclose(filePh);
            
            app.logLine(['Log saved at ' logFileName]);
        end

        % Button pushed function: RefreshButton
        function RefreshButtonPushed(app, event)
            app.autoMonitorStatus
            message = sprintf('At -> (X: %3.2f  Y: %3.2f  Z: %3.2f)\n', app.currentX, app.currentY, app.currentZ);
            app.logLine(message)            
        end

        % Button pushed function: StopButton
        function StopButtonPushed(app, event)
            readStatus(app.velmexPort, 'D');
            readStatus(app.velmexPort, 'K');
            readStatus(app.velmexPort, 'C');            
            app.autoMonitorStatus           
            message = sprintf('Manual STOP at -> (X: %3.2f  Y: %3.2f  Z: %3.2f)\n', app.currentX, app.currentY, app.currentZ);
            app.logLine(message)             
        end

        % Button pushed function: AutoCenterButton
        function AutoCenterButtonPushed(app, event)
            app.logLine('Starting auto center procedure...')
            app.setStatus(2);
            limits = autoCenter(app.velmexPort);
            while(true)
                app.autoMonitorStatus;
                if(app.status == 1)
                    break;
                end
                pause(0.5);
            end
            app.updateLimits(limits);            
            app.autoMonitorStatus;
            app.updateValues;
            app.logLine('Finished auto center procedure.')
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [0.651 0.651 0.651];
            app.UIFigure.Position = [100 100 648 559];
            app.UIFigure.Name = 'UI Figure';

            % Create LimitsPanel
            app.LimitsPanel = uipanel(app.UIFigure);
            app.LimitsPanel.TitlePosition = 'centertop';
            app.LimitsPanel.Title = 'Limits';
            app.LimitsPanel.FontWeight = 'bold';
            app.LimitsPanel.Position = [137 394 364 152];

            % Create XMinEditFieldLabel
            app.XMinEditFieldLabel = uilabel(app.LimitsPanel);
            app.XMinEditFieldLabel.HorizontalAlignment = 'right';
            app.XMinEditFieldLabel.Position = [11 95 36 22];
            app.XMinEditFieldLabel.Text = 'X Min';

            % Create XMinEditField
            app.XMinEditField = uieditfield(app.LimitsPanel, 'numeric');
            app.XMinEditField.Position = [65 95 53 22];

            % Create XMaxEditFieldLabel
            app.XMaxEditFieldLabel = uilabel(app.LimitsPanel);
            app.XMaxEditFieldLabel.HorizontalAlignment = 'right';
            app.XMaxEditFieldLabel.Position = [11 68 39 22];
            app.XMaxEditFieldLabel.Text = 'X Max';

            % Create XMaxEditField
            app.XMaxEditField = uieditfield(app.LimitsPanel, 'numeric');
            app.XMaxEditField.Position = [65 68 53 22];

            % Create YMinEditFieldLabel
            app.YMinEditFieldLabel = uilabel(app.LimitsPanel);
            app.YMinEditFieldLabel.HorizontalAlignment = 'right';
            app.YMinEditFieldLabel.Position = [128 95 36 22];
            app.YMinEditFieldLabel.Text = 'Y Min';

            % Create YMinEditField
            app.YMinEditField = uieditfield(app.LimitsPanel, 'numeric');
            app.YMinEditField.Position = [182 95 53 22];

            % Create YMaxEditFieldLabel
            app.YMaxEditFieldLabel = uilabel(app.LimitsPanel);
            app.YMaxEditFieldLabel.HorizontalAlignment = 'right';
            app.YMaxEditFieldLabel.Position = [127 68 40 22];
            app.YMaxEditFieldLabel.Text = 'Y Max';

            % Create YMaxEditField
            app.YMaxEditField = uieditfield(app.LimitsPanel, 'numeric');
            app.YMaxEditField.Position = [182 68 53 22];

            % Create ZMinEditFieldLabel
            app.ZMinEditFieldLabel = uilabel(app.LimitsPanel);
            app.ZMinEditFieldLabel.HorizontalAlignment = 'right';
            app.ZMinEditFieldLabel.Position = [244 95 36 22];
            app.ZMinEditFieldLabel.Text = 'Z Min';

            % Create ZMinEditField
            app.ZMinEditField = uieditfield(app.LimitsPanel, 'numeric');
            app.ZMinEditField.Position = [298 95 53 22];

            % Create ZMaxEditFieldLabel
            app.ZMaxEditFieldLabel = uilabel(app.LimitsPanel);
            app.ZMaxEditFieldLabel.HorizontalAlignment = 'right';
            app.ZMaxEditFieldLabel.Position = [244 68 39 22];
            app.ZMaxEditFieldLabel.Text = 'Z Max';

            % Create ZMaxEditField
            app.ZMaxEditField = uieditfield(app.LimitsPanel, 'numeric');
            app.ZMaxEditField.Position = [298 68 53 22];

            % Create ClearButton
            app.ClearButton = uibutton(app.LimitsPanel, 'push');
            app.ClearButton.ButtonPushedFcn = createCallbackFcn(app, @ClearButtonPushed, true);
            app.ClearButton.Position = [34 13 62 22];
            app.ClearButton.Text = 'Clear';

            % Create SetLimitsButton
            app.SetLimitsButton = uibutton(app.LimitsPanel, 'push');
            app.SetLimitsButton.ButtonPushedFcn = createCallbackFcn(app, @SetLimitsButtonPushed, true);
            app.SetLimitsButton.Position = [263 13 68 22];
            app.SetLimitsButton.Text = 'Set Limits';

            % Create AutoCenterButton
            app.AutoCenterButton = uibutton(app.LimitsPanel, 'push');
            app.AutoCenterButton.ButtonPushedFcn = createCallbackFcn(app, @AutoCenterButtonPushed, true);
            app.AutoCenterButton.Position = [128 13 98 22];
            app.AutoCenterButton.Text = 'Auto Center';

            % Create GotoPanel
            app.GotoPanel = uipanel(app.UIFigure);
            app.GotoPanel.TitlePosition = 'centertop';
            app.GotoPanel.Title = 'Go to...';
            app.GotoPanel.FontWeight = 'bold';
            app.GotoPanel.Position = [514 394 124 152];

            % Create XEditFieldLabel
            app.XEditFieldLabel = uilabel(app.GotoPanel);
            app.XEditFieldLabel.HorizontalAlignment = 'right';
            app.XEditFieldLabel.Position = [12 104 25 22];
            app.XEditFieldLabel.Text = 'X';

            % Create XEditField
            app.XEditField = uieditfield(app.GotoPanel, 'numeric');
            app.XEditField.Position = [55 104 53 22];

            % Create YEditFieldLabel
            app.YEditFieldLabel = uilabel(app.GotoPanel);
            app.YEditFieldLabel.HorizontalAlignment = 'right';
            app.YEditFieldLabel.Position = [12 75 25 22];
            app.YEditFieldLabel.Text = 'Y';

            % Create YEditField
            app.YEditField = uieditfield(app.GotoPanel, 'numeric');
            app.YEditField.Position = [55 75 53 22];

            % Create ZEditFieldLabel
            app.ZEditFieldLabel = uilabel(app.GotoPanel);
            app.ZEditFieldLabel.HorizontalAlignment = 'right';
            app.ZEditFieldLabel.Position = [12 46 25 22];
            app.ZEditFieldLabel.Text = 'Z';

            % Create ZEditField
            app.ZEditField = uieditfield(app.GotoPanel, 'numeric');
            app.ZEditField.Position = [55 46 53 22];

            % Create GoButton
            app.GoButton = uibutton(app.GotoPanel, 'push');
            app.GoButton.ButtonPushedFcn = createCallbackFcn(app, @GoButtonPushed, true);
            app.GoButton.FontWeight = 'bold';
            app.GoButton.Position = [12 13 100 23];
            app.GoButton.Text = 'Go';

            % Create StepmovementPanel
            app.StepmovementPanel = uipanel(app.UIFigure);
            app.StepmovementPanel.TitlePosition = 'centertop';
            app.StepmovementPanel.Title = 'Step movement';
            app.StepmovementPanel.BackgroundColor = [1 1 1];
            app.StepmovementPanel.FontWeight = 'bold';
            app.StepmovementPanel.Position = [13 80 438 301];

            % Create Image
            app.Image = uiimage(app.StepmovementPanel);
            app.Image.Position = [53 5 316 280];
            app.Image.ImageSource = 'axis.png';

            % Create XButton_minus
            app.XButton_minus = uibutton(app.StepmovementPanel, 'push');
            app.XButton_minus.ButtonPushedFcn = createCallbackFcn(app, @XButton_minusPushed, true);
            app.XButton_minus.FontWeight = 'bold';
            app.XButton_minus.Position = [53 92 40 22];
            app.XButton_minus.Text = 'X-';

            % Create ZButton_plus
            app.ZButton_plus = uibutton(app.StepmovementPanel, 'push');
            app.ZButton_plus.ButtonPushedFcn = createCallbackFcn(app, @ZButton_plusPushed, true);
            app.ZButton_plus.FontWeight = 'bold';
            app.ZButton_plus.Position = [329 93 40 22];
            app.ZButton_plus.Text = 'Z+';

            % Create YButton_plus
            app.YButton_plus = uibutton(app.StepmovementPanel, 'push');
            app.YButton_plus.ButtonPushedFcn = createCallbackFcn(app, @YButton_plusPushed, true);
            app.YButton_plus.FontWeight = 'bold';
            app.YButton_plus.Position = [192 10 40 22];
            app.YButton_plus.Text = 'Y+';

            % Create YButton_minus
            app.YButton_minus = uibutton(app.StepmovementPanel, 'push');
            app.YButton_minus.ButtonPushedFcn = createCallbackFcn(app, @YButton_minusPushed, true);
            app.YButton_minus.FontWeight = 'bold';
            app.YButton_minus.Position = [192 257 40 22];
            app.YButton_minus.Text = 'Y-';

            % Create ZButton_minus
            app.ZButton_minus = uibutton(app.StepmovementPanel, 'push');
            app.ZButton_minus.ButtonPushedFcn = createCallbackFcn(app, @ZButton_minusPushed, true);
            app.ZButton_minus.FontWeight = 'bold';
            app.ZButton_minus.Position = [82 185 40 22];
            app.ZButton_minus.Text = 'Z-';

            % Create XButton_plus
            app.XButton_plus = uibutton(app.StepmovementPanel, 'push');
            app.XButton_plus.ButtonPushedFcn = createCallbackFcn(app, @XButton_plusPushed, true);
            app.XButton_plus.FontWeight = 'bold';
            app.XButton_plus.Position = [304 182 40 22];
            app.XButton_plus.Text = 'X+';

            % Create StepsizecmEditFieldLabel
            app.StepsizecmEditFieldLabel = uilabel(app.StepmovementPanel);
            app.StepsizecmEditFieldLabel.HorizontalAlignment = 'right';
            app.StepsizecmEditFieldLabel.FontWeight = 'bold';
            app.StepsizecmEditFieldLabel.Position = [245 43 85 22];
            app.StepsizecmEditFieldLabel.Text = 'Step size (cm)';

            % Create StepsizecmEditField
            app.StepsizecmEditField = uieditfield(app.StepmovementPanel, 'numeric');
            app.StepsizecmEditField.FontWeight = 'bold';
            app.StepsizecmEditField.Position = [340 43 48 22];

            % Create RefreshButton
            app.RefreshButton = uibutton(app.StepmovementPanel, 'push');
            app.RefreshButton.ButtonPushedFcn = createCallbackFcn(app, @RefreshButtonPushed, true);
            app.RefreshButton.Position = [12 10 100 22];
            app.RefreshButton.Text = 'Refresh';

            % Create StopButton
            app.StopButton = uibutton(app.StepmovementPanel, 'push');
            app.StopButton.ButtonPushedFcn = createCallbackFcn(app, @StopButtonPushed, true);
            app.StopButton.FontWeight = 'bold';
            app.StopButton.Position = [10 39 100 23];
            app.StopButton.Text = 'Stop';

            % Create CurrentpositionPanel
            app.CurrentpositionPanel = uipanel(app.UIFigure);
            app.CurrentpositionPanel.TitlePosition = 'centertop';
            app.CurrentpositionPanel.Title = 'Current position';
            app.CurrentpositionPanel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.CurrentpositionPanel.FontWeight = 'bold';
            app.CurrentpositionPanel.Position = [467 80 172 301];

            % Create XLabel
            app.XLabel = uilabel(app.CurrentpositionPanel);
            app.XLabel.HorizontalAlignment = 'right';
            app.XLabel.FontSize = 18;
            app.XLabel.Position = [36 248 25 22];
            app.XLabel.Text = 'X:';

            % Create currentXvalue
            app.currentXvalue = uieditfield(app.CurrentpositionPanel, 'numeric');
            app.currentXvalue.Editable = 'off';
            app.currentXvalue.FontSize = 18;
            app.currentXvalue.Position = [76 247 63 23];

            % Create minEditField_2Label
            app.minEditField_2Label = uilabel(app.CurrentpositionPanel);
            app.minEditField_2Label.HorizontalAlignment = 'right';
            app.minEditField_2Label.FontSize = 15;
            app.minEditField_2Label.FontColor = [0.6353 0.0784 0.1843];
            app.minEditField_2Label.Position = [-15 216 51 22];
            app.minEditField_2Label.Text = 'min';

            % Create minXvalue
            app.minXvalue = uieditfield(app.CurrentpositionPanel, 'numeric');
            app.minXvalue.Editable = 'off';
            app.minXvalue.FontSize = 15;
            app.minXvalue.FontColor = [0.6353 0.0784 0.1843];
            app.minXvalue.Position = [43 215 38 23];

            % Create maxEditField_4Label
            app.maxEditField_4Label = uilabel(app.CurrentpositionPanel);
            app.maxEditField_4Label.HorizontalAlignment = 'right';
            app.maxEditField_4Label.FontSize = 15;
            app.maxEditField_4Label.FontColor = [0.6353 0.0784 0.1843];
            app.maxEditField_4Label.Position = [64 215 51 22];
            app.maxEditField_4Label.Text = 'max';

            % Create maxXvalue
            app.maxXvalue = uieditfield(app.CurrentpositionPanel, 'numeric');
            app.maxXvalue.Editable = 'off';
            app.maxXvalue.FontSize = 15;
            app.maxXvalue.FontColor = [0.6353 0.0784 0.1843];
            app.maxXvalue.Position = [125 215 32 22];

            % Create YLabel
            app.YLabel = uilabel(app.CurrentpositionPanel);
            app.YLabel.HorizontalAlignment = 'right';
            app.YLabel.FontSize = 18;
            app.YLabel.Position = [37 171 25 22];
            app.YLabel.Text = 'Y:';

            % Create currentYvalue
            app.currentYvalue = uieditfield(app.CurrentpositionPanel, 'numeric');
            app.currentYvalue.Editable = 'off';
            app.currentYvalue.FontSize = 18;
            app.currentYvalue.Position = [77 170 62 23];

            % Create EditFieldLabel_2
            app.EditFieldLabel_2 = uilabel(app.CurrentpositionPanel);
            app.EditFieldLabel_2.HorizontalAlignment = 'right';
            app.EditFieldLabel_2.FontSize = 15;
            app.EditFieldLabel_2.FontColor = [0.6353 0.0784 0.1843];
            app.EditFieldLabel_2.Position = [65 139 51 22];
            app.EditFieldLabel_2.Text = 'max';

            % Create maxYvalue
            app.maxYvalue = uieditfield(app.CurrentpositionPanel, 'numeric');
            app.maxYvalue.Editable = 'off';
            app.maxYvalue.FontSize = 15;
            app.maxYvalue.FontColor = [0.6353 0.0784 0.1843];
            app.maxYvalue.Position = [126 139 32 22];

            % Create minEditField_3Label
            app.minEditField_3Label = uilabel(app.CurrentpositionPanel);
            app.minEditField_3Label.HorizontalAlignment = 'right';
            app.minEditField_3Label.FontSize = 15;
            app.minEditField_3Label.FontColor = [0.6353 0.0784 0.1843];
            app.minEditField_3Label.Position = [-14 139 51 22];
            app.minEditField_3Label.Text = 'min';

            % Create minYvalue
            app.minYvalue = uieditfield(app.CurrentpositionPanel, 'numeric');
            app.minYvalue.Editable = 'off';
            app.minYvalue.FontSize = 15;
            app.minYvalue.FontColor = [0.6353 0.0784 0.1843];
            app.minYvalue.Position = [44 138 37 23];

            % Create ZLabel
            app.ZLabel = uilabel(app.CurrentpositionPanel);
            app.ZLabel.HorizontalAlignment = 'right';
            app.ZLabel.FontSize = 18;
            app.ZLabel.Position = [37 93 25 22];
            app.ZLabel.Text = 'Z:';

            % Create currentZvalue
            app.currentZvalue = uieditfield(app.CurrentpositionPanel, 'numeric');
            app.currentZvalue.Editable = 'off';
            app.currentZvalue.FontSize = 18;
            app.currentZvalue.Position = [77 92 62 23];

            % Create EditFieldLabel_3
            app.EditFieldLabel_3 = uilabel(app.CurrentpositionPanel);
            app.EditFieldLabel_3.HorizontalAlignment = 'right';
            app.EditFieldLabel_3.FontSize = 15;
            app.EditFieldLabel_3.FontColor = [0.6353 0.0784 0.1843];
            app.EditFieldLabel_3.Position = [65 61 51 22];
            app.EditFieldLabel_3.Text = 'max';

            % Create maxZvalue
            app.maxZvalue = uieditfield(app.CurrentpositionPanel, 'numeric');
            app.maxZvalue.Editable = 'off';
            app.maxZvalue.FontSize = 15;
            app.maxZvalue.FontColor = [0.6353 0.0784 0.1843];
            app.maxZvalue.Position = [126 61 32 22];

            % Create minEditField_4Label
            app.minEditField_4Label = uilabel(app.CurrentpositionPanel);
            app.minEditField_4Label.HorizontalAlignment = 'right';
            app.minEditField_4Label.FontSize = 15;
            app.minEditField_4Label.FontColor = [0.6353 0.0784 0.1843];
            app.minEditField_4Label.Position = [-14 61 51 22];
            app.minEditField_4Label.Text = 'min';

            % Create minZvalue
            app.minZvalue = uieditfield(app.CurrentpositionPanel, 'numeric');
            app.minZvalue.Editable = 'off';
            app.minZvalue.FontSize = 15;
            app.minZvalue.FontColor = [0.6353 0.0784 0.1843];
            app.minZvalue.Position = [44 60 37 23];

            % Create ZerohereButton
            app.ZerohereButton = uibutton(app.CurrentpositionPanel, 'push');
            app.ZerohereButton.ButtonPushedFcn = createCallbackFcn(app, @ZerohereButtonPushed, true);
            app.ZerohereButton.Position = [35 17 100 22];
            app.ZerohereButton.Text = 'Zero here';

            % Create Image2
            app.Image2 = uiimage(app.UIFigure);
            app.Image2.Position = [19 470 100 100];
            app.Image2.ImageSource = 'LogoGFN.png';

            % Create Image3
            app.Image3 = uiimage(app.UIFigure);
            app.Image3.Position = [23 385 91 88];
            app.Image3.ImageSource = 'LogoUCM.jpg';

            % Create DSanchezSEspaaLabel
            app.DSanchezSEspaaLabel = uilabel(app.UIFigure);
            app.DSanchezSEspaaLabel.HorizontalAlignment = 'center';
            app.DSanchezSEspaaLabel.FontSize = 10;
            app.DSanchezSEspaaLabel.Position = [6 470 125 22];
            app.DSanchezSEspaaLabel.Text = 'D Sanchez / S España';

            % Create LogTextAreaLabel
            app.LogTextAreaLabel = uilabel(app.UIFigure);
            app.LogTextAreaLabel.HorizontalAlignment = 'right';
            app.LogTextAreaLabel.FontWeight = 'bold';
            app.LogTextAreaLabel.Position = [148 43 27 22];
            app.LogTextAreaLabel.Text = 'Log';

            % Create LogTextArea
            app.LogTextArea = uitextarea(app.UIFigure);
            app.LogTextArea.Editable = 'off';
            app.LogTextArea.Position = [190 7 449 60];

            % Create StagestatusPanel
            app.StagestatusPanel = uipanel(app.UIFigure);
            app.StagestatusPanel.TitlePosition = 'centertop';
            app.StagestatusPanel.Title = 'Stage status ';
            app.StagestatusPanel.FontWeight = 'bold';
            app.StagestatusPanel.Position = [13 7 123 60];

            % Create ConnectedLampLabel
            app.ConnectedLampLabel = uilabel(app.StagestatusPanel);
            app.ConnectedLampLabel.HorizontalAlignment = 'right';
            app.ConnectedLampLabel.Position = [12 12 65 22];
            app.ConnectedLampLabel.Text = 'Connected';

            % Create ConnectedLamp
            app.ConnectedLamp = uilamp(app.StagestatusPanel);
            app.ConnectedLamp.Position = [92 12 20 20];

            % Create SaveButton
            app.SaveButton = uibutton(app.UIFigure, 'push');
            app.SaveButton.ButtonPushedFcn = createCallbackFcn(app, @SaveButtonPushed, true);
            app.SaveButton.Position = [144 18 40 22];
            app.SaveButton.Text = 'Save';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = private)

        % Construct app
        function app = stageControl_exported(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

            if nargout == 0
                clear app
            end
        end
    end
    
    methods (Access = public)

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end