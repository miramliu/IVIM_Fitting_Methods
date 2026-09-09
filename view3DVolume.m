function view3DVolume(vol)
    % VIEW3DVOLUME View a 3D image volume with slice and independent min/max contrast sliders and labels.
    % Now supports Left and Right arrow keys for changing slices.
    %
    % Syntax:
    %   view3DVolume(vol)
    %
    % Input:
    %   vol - A 3D matrix of size (nx, ny, slices)
    
    % Check input dimensions
    if ndims(vol) ~= 3
        error('Input volume must be a 3D array with dimensions (nx, ny, slices).');
    end
    
    [nx, ny, num_slices] = size(vol);
    
    % Determine the absolute data range for the sliders
    absolute_min = double(min(vol(:)));
    absolute_max = double(max(vol(:)));
    if absolute_min == absolute_max
        absolute_max = absolute_min + 1; % Prevent flat range errors
    end
    
    % Set initial states
    current_slice = round(num_slices / 2);
    current_min_clim = absolute_min;
    current_max_clim = absolute_max;
    
    % Create the main figure window
    fig = figure('Name', '3D Volume Viewer', ...
                 'NumberTitle', 'off', ...
                 'MenuBar', 'none', ...
                 'ToolBar', 'figure', ...
                 'WindowKeyPressFcn', @handleKeyPress); % <-- Added keypress listener
             
    % Create the image axes (Adjusted margins to fit left and right sliders and texts)
    ax = axes('Parent', fig, ...
              'Units', 'normalized', ...
              'Position', [0.15, 0.18, 0.70, 0.72]);
    
    % Display the initial slice
    img_handle = imagesc(vol(:, :, current_slice), 'Parent', ax);
    colormap(ax, 'jet'); 
    colorbar(ax);
    axis(ax, 'image'); 
    clim(ax, [current_min_clim, current_max_clim]); 
    title(ax, sprintf('Slice: %d / %d', current_slice, num_slices), 'FontSize', 12);
    
    % -------------------------------------------------------------------------
    % --- Bottom Slider: Slice Selection ---
    % -------------------------------------------------------------------------
    slice_slider = uicontrol('Parent', fig, ...
        'Style', 'slider', ...
        'Units', 'normalized', ...
        'Position', [0.15, 0.08, 0.70, 0.04], ...
        'Min', 1, 'Max', num_slices, 'Value', current_slice, ...
        'SliderStep', [1/(num_slices-1), 10/(num_slices-1)], ...
        'Callback', @updateSliceSlider);
    
    % Slice Value Labels (Left, Center/Current, Right)
    uicontrol('Parent', fig, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.12, 0.08, 0.03, 0.03], 'String', '1', 'HorizontalAlignment', 'right');
    
    uicontrol('Parent', fig, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.85, 0.08, 0.04, 0.03], 'String', num2str(num_slices), 'HorizontalAlignment', 'left');
    
    slice_text = uicontrol('Parent', fig, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.15, 0.03, 0.70, 0.04], ...
        'String', sprintf('Current Slice: %d (Use Left/Right Arrows)', current_slice), 'FontWeight', 'bold');

    % -------------------------------------------------------------------------
    % --- Left Slider: Minimum Contrast Adjustment ---
    % -------------------------------------------------------------------------
    min_contrast_slider = uicontrol('Parent', fig, ...
        'Style', 'slider', ...
        'Units', 'normalized', ...
        'Position', [0.05, 0.18, 0.04, 0.72], ...
        'Min', absolute_min, 'Max', absolute_max, 'Value', current_min_clim, ...
        'Callback', @updateMinContrast);
    
    % Left Slider Labels (Title, Top/Max, Bottom/Min, Current)
    uicontrol('Parent', fig, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.02, 0.92, 0.10, 0.03], 'String', 'Min Contrast', 'FontWeight', 'bold');
    
    uicontrol('Parent', fig, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.01, 0.87, 0.04, 0.03], 'String', num2str(absolute_max, '%.2g'), 'HorizontalAlignment', 'right');
    
    uicontrol('Parent', fig, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.01, 0.18, 0.04, 0.03], 'String', num2str(absolute_min, '%.2g'), 'HorizontalAlignment', 'right');
    
    min_text = uicontrol('Parent', fig, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.02, 0.13, 0.10, 0.03], ...
        'String', sprintf('Val: %.2f', current_min_clim), 'ForegroundColor', 'blue');

    % -------------------------------------------------------------------------
    % --- Right Slider: Maximum Contrast Adjustment ---
    % -------------------------------------------------------------------------
    max_contrast_slider = uicontrol('Parent', fig, ...
        'Style', 'slider', ...
        'Units', 'normalized', ...
        'Position', [0.91, 0.18, 0.04, 0.72], ...
        'Min', absolute_min, 'Max', absolute_max, 'Value', current_max_clim, ...
        'Callback', @updateMaxContrast);
    
    % Right Slider Labels (Title, Top/Max, Bottom/Min, Current)
    uicontrol('Parent', fig, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.88, 0.92, 0.10, 0.03], 'String', 'Max Contrast', 'FontWeight', 'bold');
    
    uicontrol('Parent', fig, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.95, 0.87, 0.04, 0.03], 'String', num2str(absolute_max, '%.2g'), 'HorizontalAlignment', 'left');
    
    uicontrol('Parent', fig, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.95, 0.18, 0.04, 0.03], 'String', num2str(absolute_min, '%.2g'), 'HorizontalAlignment', 'left');
    
    max_text = uicontrol('Parent', fig, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.88, 0.13, 0.10, 0.03], ...
        'String', sprintf('Val: %.2f', current_max_clim), 'ForegroundColor', 'red');

    % -------------------------------------------------------------------------
    % --- Callback Functions ---
    % -------------------------------------------------------------------------
    function updateImageDisplay()
        % Centralized function to refresh image display and label
        set(img_handle, 'CData', vol(:, :, current_slice));
        title(ax, sprintf('Slice: %d / %d', current_slice, num_slices));
        set(slice_text, 'String', sprintf('Current Slice: %d (Use Left/Right Arrows)', current_slice));
    end

    function updateSliceSlider(src, ~)
        % Slider adjustment callback
        current_slice = round(get(src, 'Value'));
        updateImageDisplay();
    end

    function handleKeyPress(~, event)
        % Keyboard arrow key callback
        switch event.Key
            case 'leftarrow'
                if current_slice > 1
                    current_slice = current_slice - 1;
                    set(slice_slider, 'Value', current_slice); % Sync slider handle
                    updateImageDisplay();
                end
            case 'rightarrow'
                if current_slice < num_slices
                    current_slice = current_slice + 1;
                    set(slice_slider, 'Value', current_slice); % Sync slider handle
                    updateImageDisplay();
                end
        end
    end

    function updateMinContrast(src, ~)
        val = get(src, 'Value');
        % Safeguard: prevent Min from exceeding or equaling Max
        if val >= current_max_clim
            val = current_max_clim - 0.001 * (absolute_max - absolute_min);
            set(src, 'Value', val); 
        end
        current_min_clim = val;
        clim(ax, [current_min_clim, current_max_clim]);
        set(min_text, 'String', sprintf('Val: %.2f', current_min_clim));
    end

    function updateMaxContrast(src, ~)
        val = get(src, 'Value');
        % Safeguard: prevent Max from falling below or equaling Min
        if val <= current_min_clim
            val = current_min_clim + 0.001 * (absolute_max - absolute_min);
            set(src, 'Value', val); 
        end
        current_max_clim = val;
        clim(ax, [current_min_clim, current_max_clim]);
        set(max_text, 'String', sprintf('Val: %.2f', current_max_clim));
    end
end