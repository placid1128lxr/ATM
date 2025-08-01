%% Stanford-ESA Integrity Diagram
%  Version : 1.0.0.0
%  Author  : E. Ogier
%  Release : 25th january 2021
classdef StanfordDiagram < matlab.mixin.SetGet
    
    % Private properties
    properties (SetAccess = 'private')
        
        % Histogram data
        TotalEpochs   = 0;                                  % Total epochs, all domains included
        Epochs        = zeros(1,7);                         % Epochs by domain
        Statistics    = zeros(1,7);                         % Statistics by domain (normalized epochs, between 0 and 1)
        MatrixEpochs  = [];                                 % Matrix of epochs (data of histogram)
        
        % Histogram setting
        CategoryLimit = [];                                 % Category limit (if relevant)
        AlertLimit    = [];                                 % Alert limit (if relevant)
        Minimum       = 0;                                  % Minimum value of PE/PL
        Maximum       = NaN;                                % Maximum value of PE/PL
        Number        = NaN;                                % Number of edges of PE/PL 
        Step          = NaN;                                % Step of PE/PL
        Edges         = NaN;                                % Edges of PE/PL
        
        % Graphical elements
        Axes          = [];                                 % Axes of diagram (new axes if not defined)
        Scale         = 'auto';                             % Scale of color bar ('auto', 'log' or 'linear')
        Title         = 'Stanford-ESA Integrity Diagram';   % Title of diagram
        Xlabel        = 'Position error [m]';               % Label of diagram X-axis
        Ylabel        = 'Protection level [m]';             % Label of diagram Y-axis
        Xtick         = [];                                 % Ticks of diagram X-axis
        Xticklabel    = '';                                 % Tick labels of diagram X-axis
        Ytick         = [];                                 % Ticks of diagram Y-axis
        Yticklabel    = '';                                 % Tick labels of diagram Y-axis
        FontSize      = 9;                                  % Fonsize of all labels
        
        % Domains settings
        Strings       = ...                                                                         % Strings of diagram domains:
            {'Nominal Operations #1',...                                                            % - String of domain #1: nominal operations #1
             'Nominal Operations #2',...                                                            % - String of domain #2: nominal operations #2
             'System Unavailable',...                                                               % - String of domain #3: system unavailable
             'Misleading Informations #1',...                                                       % - String of domain #4: misleading informations #1
             'Misleading Informations #2',...                                                       % - String of domain #5: misleading informations #2
             ['System Unavailable and',char(10),'Misleading Informations'],...                      % - String of domain #6: system unavailable and misleading informations
             'Hazardous Operations'};                                                               % - String of domain #7: hazardous operations
        StringFcns    = ...                                                                         % String functions of diagram domains
            {@(e,s)[char(10) sprintf('epochs: %u',e),char(10),sprintf('stat.: %.5f%%',100*s)],...   % - String function of domain #1: nominal operations #1
             @(e,s)[char(10) sprintf('epochs: %u',e),char(10),sprintf('stat.: %.5f%%',100*s)],...   % - String function of domain #2: nominal operations #2
             @(e,s)[char(10) sprintf('epochs: %u',e),char(10),sprintf('stat.: %.5f%%',100*s)],...   % - String function of domain #3: system unavailable
             @(e,s)[char(10) sprintf('epochs: %u',e),char(10),sprintf('stat.: %.5f%%',100*s)],...   % - String function of domain #4: misleading informations #1
             @(e,s)[char(10) sprintf('epochs: %u',e),char(10),sprintf('stat.: %.5f%%',100*s)],...   % - String function of domain #5: misleading informations #2
             @(e,s)[char(10) sprintf('epochs: %u',e),char(10),sprintf('stat.: %.5f%%',100*s)],...   % - String function of domain #6: system unavailable and misleading informations
             @(e,s)[char(10) sprintf('epochs: %u',e),char(10),sprintf('stat.: %.5f%%',100*s)]};     % - String function of domain #7: hazardous operations  
        Colors        = ...                                                                         % Colors of diagram domains:
            {[255 255 255]/255,...                                                                  % - Colors of domain #1: nominal operations #1
             [255 255 255]/255,...                                                                  % - Colors of domain #2: nominal operations #1
             [255 255   0]/255,...                                                                  % - Colors of domain #3: system unavailable
             [255 204 153]/255,...                                                                  % - Colors of domain #4: misleading informations #1
             [255 204 153]/255,...                                                                  % - Colors of domain #5: misleading informations #2
             [254 153   0]/255,...                                                                  % - Colors of domain #6: system unavailable and misleading informations
             [255   0   0]/255};                                                                    % - Colors of domain #7: hazardous operations  
        
    end
    
    % Hidden properties
    properties (Hidden)
        
        % Graphical elements
        Histogram     = [];                             % Graphical histogram
        Domains         = [];                           % Graphical domains
        Labels        = [];                             % Labels of domains
        Colorbar      = [];                             % color bar of histogram
        
    end
    
    % Constant properties
    properties (Constant)
        
        % Constants
        NumberOfDomains  = 7;                           % Number of domains of diagram
        FaceAlpha1     = 0.9;                           % Transparency factor of histogram
        FaceAlpha2     = 0.5;                           % Transparency of domain patchs
        DefaultMaximum = 50;                            % Default value of maximum PE/PL
        DefaultNumber  = 100;                           % Default number of steps of PE/PL
        
    end
    
    % Public methods
    methods (Access = 'public')
        
        % Constructor
        function Object = StanfordDiagram(varargin)
            
            % Default position error and protection level
            PositionError   = [];
            ProtectionLevel = [];
            
            % Update of properties
            for v = 1:2:length(varargin)
                
                switch lower(varargin{v})
                    case {'error','positionerror'}
                        PositionError = varargin{v+1};
                    case {'level','protectionlevel'}
                        ProtectionLevel = varargin{v+1};
                    otherwise
                        set(Object,varargin{v},varargin{v+1});
                end
                
            end
            
            % Application of default values for maximum, number and step
            if isnan(Object.Number)
                Object.Number = Object.DefaultNumber;
            end
            if isnan(Object.Maximum)
                Object.Maximum = max([max(PositionError),max(ProtectionLevel)]);
                if isempty(Object.Maximum)
                    Object.Maximum = Object.DefaultMaximum;
                    warning('Maximum value for both position error and protection level arbitrary set to 50 meters.');
                end
                Object.Step = Object.Maximum/Object.Number;
            end
            
            % Computation of edges
            Object.Edges = Object.Minimum:Object.Step:Object.Maximum;  
            
            % Rectification in case of inconstancy between step and
            % boundary
            if ~isnan(Object.Maximum)
                if gt(Object.Maximum,Object.Edges(end))
                    Object.Edges(end+1) = Object.Edges(end)+Object.Step;
                    Object.Maximum = Object.Edges(end);
                end
            end
            
            % Show if available data
            if ~isempty(PositionError) && ~isempty(ProtectionLevel)
                Object.show(PositionError,ProtectionLevel);
            end
            
        end
        
        % Method 'set'
        function Object = set(Object,varargin)
            
            Properties = varargin(1:2:end);
            Values     = varargin(2:2:end);
            
            PropertiesObject = properties(Object);
            
            for n = 1:length(Properties)
                
                I = strcmpi(Properties{n},PropertiesObject);
                
                if any(I)
                    
                    % Control of values
                    switch lower(Properties{n})
                        
                        case {'totalepochs','epochs','statistics','matrixepochs'}
                            
                            % Warning
                            if ~isempty(Object.MatrixEpochs)
                                warning('Property ''%s'' can be reset through method ''reset''.',PropertiesObject{I});
                            end
                            
                        case {'categorylimit','alertlimit','minimum','maximum','step'}
                            
                            % Empty, numeric, not NaN and not infinite value
                            if ~isempty(Values{n})
                                if ~isnumeric(Values{n}) || isnan(Values{n}) || isinf(Values{n})
                                    error('Property ''%s'' shall be numeric.',PropertiesObject{I});
                                end
                            end
                            
                        case 'edges'
                            
                            % Valid edges i.e. vector of linearly increasing values
                            Steps = unique(diff(Values{n}));
                            if isscalar(Steps)
                                error('Property ''Edges'' shall be a linearly increasing vector.');
                            end
                            Object.Minimum = Values{n}(1);
                            Object.Maximum = Values{n}(2);
                            Object.Step    = Steps;                         
                        
                        case 'axes'
                            
                            % Valid axes
                            if ~isa(Values{n},'matlab.graphics.axis.Axes')
                                error('Property ''Axes'' shall be valid axes.');                            
                            end
                            
                        case 'scale'
                            
                            % Scale 'Log' or 'Linear'
                            if ~ischar(Values{n})
                                error('Property ''Scale'' shall be ''Log'' or ''Linear''.');
                            end
                            if ~ismember(lower(Values{n}),{'auto','log','linear'})
                                error('Property ''Scale'' shall be ''Log'' or ''Linear''.');
                            end
                            
                        case {'title','xlabel','ylabel'}
                            
                            % String
                            if ~ischar(Values{n})
                                error('Property ''%s'' shall be a string.',PropertiesObject{I});
                            end
                            
                        case {'xtick','ytick'}
                            
                            % Monotonic numeric vector
                            if ~isnumeric(Values{n})
                                error('Property ''%s'' shall be a numeric vector.',PropertiesObject{I});                                
                            elseif ~all(gt(Values{n},0))                     
                                error('Property ''%s'' shall be a monotonic vector.',PropertiesObject{I});          
                            end
                            
                        case {'xticklabel','yticklabel'}
                            
                            % Cell
                            if ~iscell(Values{n})
                                error('Property ''%s'' shall be a cell of strings',PropertiesObject{I});
                            end
                            if ~all(strcmp(cellfun(@class,Values{n},'UniformOutput',false),'char'))
                                error('Property ''%s'' shall be a cell of strings',PropertiesObject{I});
                            end
                            
                        case {'strings','stringfcns','colors'}
                            
                            % Cell
                            if ~iscell(Values{n})
                                error('Property ''%s'' shall be a cell.',PropertiesObject{I});
                            end
                            
                            switch lower(Properties{n})
                                
                                case 'strings'
                                    
                                    % Cell of strings
                                    if ~all(strcmp(cellfun(@class,Values{n},'UniformOutput',false),'char'))                                        
                                        error('Property ''Strings'' shall be a cell of strings (use char(10) for line break.');
                                    end
                                    
                                case 'stringfcns'
                                    
                                    % Cell of anonymous functions or empty values
                                    if ~all(ismember(cellfun(@class,Values{n},'UniformOutput',false),{'char','function_handle'}))                                        
                                        error('Property ''StringFcns'' shall be a cell of anonymous functions (inputs: epochs and statistics) or empty values.');
                                    end
                                    
                                case 'colors'  
                                    
                                    % Cell of RGB values or colors names
                                    if ~all(strcmp(cellfun(@class,Values{n},'UniformOutput',false),'char') |...
                                            (cell2mat(cellfun(@isnumeric,Values{n},'UniformOutput',false)) & ...
                                             cell2mat(cellfun(@(c)eq(numel(c),3),Values{n},'UniformOutput',false)))) 
                                         error('Property ''StringFcns'' shall be a cell of vectors (3 RGB values) or strings (color names).');
                                    end
                                    
                            end
                            
                            % Reaffectation of the values relative to domains in a
                            % cell of 7 values, in the case only 2 or 5 domains are
                            % defined
                            switch numel(Values{n})
                                case 2
                                    Values2 = Object.(PropertiesObject{I});
                                    Values2{1} = Values{n}{1};
                                    Values2{4} = Values{n}{2};
                                case 5
                                    Values2 = Object.(PropertiesObject{I});
                                    if ~isempty(Object.AlertLimit)
                                        if isnumeric(Object.AlertLimit) && ~isnan(Object.AlertLimit) && ~isinf(Object.AlertLimit)
                                            % No category
                                            Values2{1} = Values{n}{1};
                                            Values2{3} = Values{n}{2};
                                            Values2{4} = Values{n}{3};
                                            Values2{6} = Values{n}{4};
                                            Values2{7} = Values{n}{5};
                                        else              
                                            % No alert limit
                                            Values2{1} = Values{n}{1};
                                            Values2{2} = Values{n}{2};
                                            Values2{4} = Values{n}{3};
                                            Values2{5} = Values{n}{4};
                                            Values2{7} = Values{n}{5};                                            
                                        end
                                    else                                        
                                        % No alert limit
                                        Values2{1} = Values{n}{1};
                                        Values2{2} = Values{n}{2};
                                        Values2{4} = Values{n}{3};
                                        Values2{5} = Values{n}{4};
                                        Values2{7} = Values{n}{5};
                                    end
                                case 7
                                    Values2 = Values{n};
                                otherwise
                                    error('Property ''%s'' shall be a cell including 2, 5 or 7 values.',PropertiesObject{I});
                            end
                            Values{n} = Values2;                       
                            
                    end
                                        
                    % Update of properties
                    Object.(PropertiesObject{I}) = Values{n};
                    
                    % Update of graphical properties
                    switch lower(PropertiesObject{I})
                        case 'xlabel'                            
                            Object.Axes.XLabel.String = Values{n};
                        case 'ylabel'                            
                            Object.Axes.YLabel.String = Values{n};
                        case 'title'                            
                            Object.Axes.Title.String = Values{n};
                        case {'xtick','xticklabel','ytick','yticklabel'}
                            set(Object.Axes,PropertiesObject{I},Values{n});
                        case 'fontsize'
                            set([Object.Axes;Object.Colorbar;findall(Object.Axes,'type','text')],'fontsize',Object.FontSize);
                        case {'scale','strings','stringfcns','colors'}                                                                               
                            if ~isempty(Object.Domains)
                                Object.updateHistogramAndColorBar();
                            end
                    end
                                            
                end
                
            end
            
            % Computation of edges
            Object.Edges = Object.Minimum:Object.Step:Object.Maximum; 
            
        end
        
        % Method 'show'
        function Object = show(Object,PositionError,ProtectionLevel)
            
            % Return
            if ~isempty(Object.Histogram)
                warning('Stanford diagram already created.');
                return
            end
            
            % Update of default inputs
            if lt(nargin,3)
                PositionError   = [];
                ProtectionLevel = [];
            end
            
            % Update of epochs and statistics
            Object.updateEpochsAndStatistics(PositionError,ProtectionLevel);
            
            % Update of axes
            if isempty(Object.Axes)
                Object.Axes = axes();
            end            
            set(Object.Axes,'NextPlot','Add');
            
            % Creation of domains
            Object = CreationDomains(Object);
            
            % Creation of histogram
            Object.Histogram = ...
                imagesc('Parent',    Object.Axes,...
                        'Xdata',     Object.Edges,...
                        'Ydata',     Object.Edges,...
                        'Cdata',     Object.MatrixEpochs,...
                        'AlphaData', ~~Object.MatrixEpochs*Object.FaceAlpha1);
            
            % Creation of frames
            for n = 1:Object.NumberOfDomains
                Xdata = get(Object.Domains(n),'Xdata');
                Ydata = get(Object.Domains(n),'Ydata');
                plot([Xdata; Xdata(1)],[Ydata; Ydata(1)],'k:');
            end
            
            % Creation of the limit protection level = error
            plot([0 Object.Edges(end)],[0 Object.Edges(end)],'k:');
            
            % Creation of labels       
            Object = CreationLabels(Object);
            
            % Creation of colormap
            Map = fliplr(hsv);
            Map = Map(1:45,:);            
            colormap(Object.Axes,Map);
            
            % Creation of color bar
            Object.Colorbar = colorbar(Object.Axes,'Location','eastoutside');
            
            % Update of histogram and color bar
            Object.updateHistogramAndColorBar();
            
            % Creation of titles
            title(Object.Axes,Object.Title);
            xlabel(Object.Axes,Object.Xlabel);
            ylabel(Object.Axes,Object.Ylabel);
            
            % Update of axes
            set([Object.Axes;Object.Colorbar;findall(Object.Axes,'type','text')],'fontsize',Object.FontSize);
            set(Object.Axes,...
                'Box',  'on',...
                'Xgrid','on',...
                'Ygrid','on',...
                'Xlim', [Object.Minimum Object.Maximum],...
                'Ylim', [Object.Minimum Object.Maximum]);
            
            % Update of X ticks
            if ~isempty(Object.Xtick)
                set(Object.Axes,'Xtick',Object.Xtick)
            end
            if ~isempty(Object.Xticklabel)
                set(Object.Axes,'Xticklabel',Object.Xticklabel);
            end
            
            % Update of Y ticks
            if ~isempty(Object.Ytick)
                set(Object.Axes,'Ytick',Object.Ytick)
            end
            if ~isempty(Object.Yticklabel)
                set(Object.Axes,'Yticklabel',Object.Yticklabel);
            end
            
            % Creation of domains
            function Object = CreationDomains(Object)
               
                Xdata{1} = [Object.Minimum Object.CategoryLimit Object.Minimum ];
                Ydata{1} = [Object.Minimum Object.CategoryLimit Object.CategoryLimit];
                
                Xdata{2} = [Object.Minimum       Object.CategoryLimit Object.AlertLimit Object.Minimum];
                Ydata{2} = [Object.CategoryLimit Object.CategoryLimit Object.AlertLimit Object.AlertLimit];
                
                Xdata{3} = [Object.Minimum    Object.AlertLimit Object.Maximum Object.Minimum];
                Ydata{3} = [Object.AlertLimit Object.AlertLimit Object.Maximum Object.Maximum];
                
                Xdata{4} = [Object.Minimum Object.CategoryLimit Object.CategoryLimit];
                Ydata{4} = [Object.Minimum Object.Minimum       Object.CategoryLimit];
                
                Xdata{5} = [Object.CategoryLimit Object.AlertLimit    Object.AlertLimit];
                Ydata{5} = [Object.CategoryLimit Object.CategoryLimit Object.AlertLimit];
                
                Xdata{6} = [Object.AlertLimit Object.Maximum    Object.Maximum];
                Ydata{6} = [Object.AlertLimit Object.AlertLimit Object.Maximum];
                
                Xdata{7} = [Object.CategoryLimit Object.Maximum   Object.Maximum    Object.AlertLimit Object.AlertLimit    Object.CategoryLimit];
                Ydata{7} = [Object.Minimum       Object.Minimum   Object.AlertLimit Object.AlertLimit Object.CategoryLimit Object.CategoryLimit ];
                
                for z = 1:Object.NumberOfDomains
                    
                    Object.Domains(z) = ...
                        patch('Parent',    Object.Axes,...
                              'Xdata',     Xdata{z},...
                              'Ydata',     Ydata{z},...
                              'FaceColor', Object.Colors{z},...
                              'FaceAlpha', Object.FaceAlpha2);
                    
                end
                
            end
            
            % Creation of labels
            function Object = CreationLabels(Object)
               
                for z = 1:Object.NumberOfDomains
                    
                    Xdata = get(Object.Domains(z),'Xdata');
                    Ydata = get(Object.Domains(z),'Ydata');
                    
                    if lt(z,7)
                        x = mean(Xdata);
                        y = mean(Ydata);
                    else
                        if eq(Object.AlertLimit,Object.Maximum)
                            x = mean([Object.CategoryLimit Object.Maximum]);
                        else
                            x = mean([Object.AlertLimit Object.Maximum]);
                        end
                        y = mean([Object.Minimum Object.CategoryLimit]);
                    end
                    
                    if all(eq(Xdata,x)) || all(eq(Ydata,y))
                        Visibility = 'off';
                    else
                        Visibility = 'on';
                    end
                    
                    String = Object.Strings{z};
                    if ~isa(String,'char')
                        String = '';
                    end
                    if isa(Object.StringFcns{z},'function_handle')
                        String = [String Object.StringFcns{z}(0,0)]; %#ok<AGROW>
                    end
                    
                    Object.Labels(z) = ...
                        text(x,y,String,...
                            'Parent',              Object.Axes,...
                            'FontSize',            Object.FontSize,...
                            'HorizontalAlignment', 'Center',...
                            'VerticalAlignment',   'Middle',...
                            'Visible',             Visibility);
                    
                end
                
            end
            
        end
        
        % Method 'update'
        function Object = update(Object,PositionError,ProtectionLevel)
            
            % Update of epochs and statistics
            Object.updateEpochsAndStatistics(PositionError,ProtectionLevel);
            
            % Update of histogram and color bar
            if ~isempty(Object.Histogram)
                Object.updateHistogramAndColorBar();
            end
            
            drawnow();
            
        end        
        
        % Method 'reset'
        function Object = reset(Object)
            
            % Reinitialization of data
            Object.MatrixEpochs = zeros(numel(Object.Edges),numel(Object.Edges));
            Object.Epochs       = zeros(1,Object.NumberOfDomains);     
            Object.Statistics   = zeros(1,Object.NumberOfDomains);     
            Object.TotalEpochs  = 0;  
            
            % Reset of labels
            for n = 1:Object.NumberOfDomains
                String = Object.Strings{n};
                if ~isa(String,'char')
                    String = '';
                end
                if isa(Object.StringFcns{n},'function_handle')
                    String = [String Object.StringFcns{n}(Object.Epochs(n),Object.Statistics(n))]; %#ok<AGROW>
                end
                set(Object.Labels(n),'String',String);
            end
            
            % Reset of color bar
            Object.Colorbar.Ticks = [];
            Object.Colorbar.Label.String = sprintf('Histogram (total epochs: %u)',Object.TotalEpochs);
            
            % Reset of histogram
            set(Object.Histogram,'AlphaData',~~Object.MatrixEpochs*Object.FaceAlpha1);
            
            % Reinitialization of the matrix of epochs
            Object.MatrixEpochs = [];
            
        end
        
        % Method 'report'
        function Report = report(Object)
                                 
            if eq(Object.CategoryLimit,Object.AlertLimit)
                if eq(Object.AlertLimit,Object.Maximum)
                    I = [1 4];
                else
                    I = [1 3 4 6 7];
                end
            else
                I = 1:Object.NumberOfDomains;
            end
            
            RowNames = strrep(Object.Strings,char(10),' ');      
           
            Table = ...
                table(Object.Epochs(I)',Object.Statistics(I)',...
                      'VariableNames',{'Epochs','Statistics'},...
                      'RowNames',Renaming(RowNames(I)));
            
            if nargout
                Report = Table;
            else                
                disp(Table);
                Report = [];
            end
            
            % Renaming of domains
            function List2 = Renaming(List1)
                
                List2 = {};
                N = 1:numel(List1);
                
                for n = N
                    if ismember(List1{n},List1(setdiff(N,n)))
                        i = 1;
                        while ismember(sprintf('%s #%u',List1{n},i),List1(setdiff(N,n))) || ...
                                ismember(sprintf('%s #%u',List1{n},i),List2)
                            i = i+1;
                        end
                        List2{n} = sprintf('%s #%u',List1{n},i); %#ok<AGROW>
                    else
                        List2{n} = List1{n}; %#ok<AGROW>
                    end
                    
                end
                
            end
            
        end        
        
    end
    % Private methods
    methods (Access = 'private')
      
        % Update of epochs and statistics
        function Object = updateEpochsAndStatistics(Object,PositionError,ProtectionLevel)
            
            % Update of default value of alert limit
            if isempty(Object.AlertLimit) || isnan(Object.AlertLimit)
                Object.AlertLimit = Object.Maximum;
            end
            
            % Update of default value of alert limit
            if isempty(Object.CategoryLimit) || isnan(Object.CategoryLimit)
                Object.CategoryLimit = Object.AlertLimit;
            end
            
            % Filtering
            I = or(isnan(PositionError),isnan(ProtectionLevel));
            PositionError(I)   = [];
            ProtectionLevel(I) = [];
            
            % Saturation
            PositionError   = max(PositionError,Object.Minimum);
            PositionError   = min(PositionError,Object.Maximum);
            ProtectionLevel = max(ProtectionLevel,Object.Minimum);
            ProtectionLevel = min(ProtectionLevel,Object.Maximum);
                        
            % Discretization
            IndicesPositionError   = discretize(PositionError,  Object.Edges,'IncludedEdge','right');
            IndicesProtectionLevel = discretize(ProtectionLevel,Object.Edges,'IncludedEdge','right');
            
            % Update of matrix
            if isempty(Object.MatrixEpochs)
                Object.MatrixEpochs = zeros(numel(Object.Edges),numel(Object.Edges));
            end
            for n = 1:numel(IndicesPositionError)
                Object.MatrixEpochs(IndicesProtectionLevel(n),IndicesPositionError(n)) = ...
                    Object.MatrixEpochs(IndicesProtectionLevel(n),IndicesPositionError(n))+1;   
            end
            
            % Update of epochs
            I = lt(PositionError,ProtectionLevel);            
            Counts(1) = sum(le(ProtectionLevel(I),Object.CategoryLimit));
            Counts(3) = sum(gt(ProtectionLevel(I),Object.AlertLimit));
            Counts(2) = sum(I)-Counts(1)-Counts(3);
            I = ~I;           
            Counts(4) = sum(lt(PositionError(I),Object.CategoryLimit));            
            Counts(5) = sum(lt(PositionError(I),Object.AlertLimit)&...
                            ge(ProtectionLevel(I),Object.CategoryLimit));            
            Counts(6) = sum(ge(ProtectionLevel(I),Object.AlertLimit));            
            Counts(7) = sum(I)-Counts(4)-Counts(5)-Counts(6);            
            for z = 1:Object.NumberOfDomains
                Object.Epochs(z) = Object.Epochs(z)+Counts(z);
            end            
            Object.TotalEpochs = sum(Object.Epochs);
            
            % Statistics
            Object.Statistics = Object.Epochs/Object.TotalEpochs;
            
        end
        
        % Method 'updateHistogramAndColorBar'
        function Object = updateHistogramAndColorBar(Object)
            
            if ~ishandle(Object.Axes)
                return
            end
            
            % Display function
            Max = max(max(Object.MatrixEpochs));
            Scale2 = lower(Object.Scale);
            switch Scale2
                case 'auto'
                    if Max < 100
                        Scale2 = 'linear';
                        fcn = @(x)x;
                    else
                        Scale2 = 'log';
                        fcn = @log10;
                    end
                case 'linear', fcn = @(x)x;
                case 'log',    fcn = @log10;
            end
            
            % Update of histogram
            set(Object.Histogram,...
                'Cdata',     fcn(Object.MatrixEpochs),...
                'AlphaData', ~~Object.MatrixEpochs*Object.FaceAlpha1);  
            
            % Update of counts
            for n = 1:Object.NumberOfDomains
                String = Object.Strings{n};
                if isa(Object.StringFcns{n},'function_handle')
                    String = [String Object.StringFcns{n}(Object.Epochs(n),Object.Statistics(n))]; %#ok<AGROW>
                end
                set(Object.Labels(n),'String',String);
            end
            
            % Update of colormap
            switch Scale2
                case 'linear'
                    if lt(Max,10)
                        I = 0:Max;
                    else
                        I = 0:(Max/10):Max;
                    end
                    Object.Colorbar.Ticks = I;
                    Object.Colorbar.TickLabels = arrayfun(@(e)sprintf('%.0f',e),I,'UniformOutput',false);
                    if gt(Max,0)
                        Object.Colorbar.Limits = [0 Max];
                    end
                case 'log'
                    I = 0:floor(log10(Max));
                    if ~isscalar(I)                        
                        Object.Colorbar.Ticks = I;
                        Object.Colorbar.TickLabels = arrayfun(@(e)sprintf('10^%.0f',e),I,'UniformOutput',false);
                    else                        
                        Object.Colorbar.Ticks = [0 log10(Max)];
                        Object.Colorbar.TickLabels = {'10^0' sprintf('10^{%0.2f}',log10(Max))};
                    end
                    if gt(Max,1)
                        Object.Colorbar.Limits = [0 log10(Max)];
                    end
            end
            Object.Colorbar.Label.String = sprintf('Histogram (total epochs: %u)',Object.TotalEpochs);
            
        end
            
    end
    
end


