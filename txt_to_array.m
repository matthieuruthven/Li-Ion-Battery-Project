function surface_t_array = txt_to_array(dirpath, idx)
    % Function to read TXT files of pouch cell surface temperature values
    % (and corresponding coordinates), arrange these as an array and then
    % save the array as a MAT file
    %
    % Author: Matthieu Ruthven (matthieu.ruthven@uni.lu)
    % Last modified: 9th February 2024
    %
    % Input arguments:
    % 1) dirpath (string): path to folder containing data from battery surface
    % temperature simulation
    % 2) idx (integer): the simulation time point
    %
    % Output:
    % 1) MAT file of surface temperature values

    % Read TXT files of temperature data
    b_body = readmatrix(fullfile(dirpath, 'Raw_Surface_T', sprintf('Surface_T_Battery_%04d.txt', idx)));
    neg_tab = readmatrix(fullfile(dirpath, 'Raw_Surface_T', sprintf('Surface_T_Neg_Tab_%04d.txt', idx)));
    pos_tab = readmatrix(fullfile(dirpath, 'Raw_Surface_T', sprintf('Surface_T_Pos_Tab_%04d.txt', idx)));
    
    % Read TXT file of model parameter values
    S = readlines(fullfile(dirpath, 'Model_Parameter_values.txt'));

    % Extract relevant information from TXT file of model parameter values
    for tmp_string = S'
        if startsWith(tmp_string, 'W_cell')
            b_body_w = split(tmp_string, '"');
            b_body_w = b_body_w{2};
            b_body_w = split(b_body_w, '[');
            b_body_w = str2num(b_body_w{1});
        end
        if startsWith(tmp_string, 'H_cell')
            b_body_l = split(tmp_string, '"');
            b_body_l = b_body_l{2};
            b_body_l = split(b_body_l, '[');
            b_body_l = str2num(b_body_l{1});
        end
        if startsWith(tmp_string, 'W_tab')
            b_tab_w = split(tmp_string, '"');
            b_tab_w = b_tab_w{2};
            b_tab_w = split(b_tab_w, '[');
            b_tab_w = str2num(b_tab_w{1});
        end
        if startsWith(tmp_string, 'H_tab')
            b_tab_l = split(tmp_string, '"');
            b_tab_l = b_tab_l{2};
            b_tab_l = split(b_tab_l, '[');
            b_tab_l = str2num(b_tab_l{1});
        end
        if startsWith(tmp_string, 'Itab')
            b_tab_off = split(tmp_string, '"');
            b_tab_off = b_tab_off{2};
            b_tab_off = split(b_tab_off, '[');
            b_tab_off = str2num(b_tab_off{1});
        end
        if startsWith(tmp_string, 'T0')
            T_amb = split(tmp_string);
            T_amb = split(T_amb{2}, '[');
            T_amb = str2num(T_amb{1});
        end
    end
    
    % For battery body
    
    % Extract x-coordinates
    x_coords = b_body(:, 1);
    
    % Extract y-coordinates
    y_coords = b_body(:, 2);
    
    % Extract temperature values
    T_vals = b_body(:, 3);
    
    % Find minimum and maximum coordinate values
    xmin = min(x_coords);
    xmax = max(x_coords);
    ymin = min(y_coords);
    ymax = max(y_coords);
    
    % Check consistency of coordinates
    assert(round(xmax - xmin) == b_body_w, 'Inconsistency in battery width')
    assert(round(ymax - ymin) == b_body_l, 'Inconsistency in battery length')
    
    % Create mesh grid
    [xq, yq] = meshgrid(round(xmin):round(xmax), round(ymin):round(ymax));
    
    % Interpolate 2D scattered data
    b_body = griddata(x_coords, y_coords, T_vals, xq, yq, "nearest");
    b_body = rot90(b_body, 2);
    
    % Plot image
    imagesc(b_body); axis image off
    
    % For negative tab
    
    % Extract x-coordinates
    x_coords = neg_tab(:, 1);
    
    % Extract y-coordinates
    y_coords = neg_tab(:, 2);
    
    % Extract temperature values
    T_vals = neg_tab(:, 3);
    
    % Find minimum and maximum coordinate values
    xmin = min(x_coords);
    xmax = max(x_coords);
    ymin = min(y_coords);
    ymax = max(y_coords);
    
    % Check consistency of coordinates
    if xmax < 0
        assert(abs(round(xmin - xmax)) == b_tab_w, 'Inconsistency in negative tab width')
    else
        assert(round(xmax - xmin) == b_tab_w, 'Inconsistency in negative tab width')
    end
    assert(round(ymax - ymin) == b_tab_l, 'Inconsistency in negative tab length')
    
    % Create mesh grid
    [xq, yq] = meshgrid(round(xmin):round(xmax), round(ymin):round(ymax));
    
    % Interpolate 2D scattered data
    neg_tab = griddata(x_coords, y_coords, T_vals, xq, yq, "nearest");
    neg_tab = rot90(neg_tab, 2);
    
    % Plot image
    % imagesc(neg_tab); axis image off
    
    % For positive tab
    
    % Extract x-coordinates
    x_coords = pos_tab(:, 1);
    
    % Extract y-coordinates
    y_coords = pos_tab(:, 2);
    
    % Extract temperature values
    T_vals = pos_tab(:, 3);
    
    % Find minimum and maximum coordinate values
    xmin = min(x_coords);
    xmax = max(x_coords);
    ymin = min(y_coords);
    ymax = max(y_coords);
    
    % Check consistency of coordinates
    if xmax < 0
        assert(abs(round(xmin - xmax)) == b_tab_w, 'Inconsistency in positive tab width')
    else
        assert(round(xmax - xmin) == b_tab_w, 'Inconsistency in positive tab width')
    end
    assert(round(ymax - ymin) == b_tab_l, 'Inconsistency in positive tab length')
    
    % Create mesh grid
    [xq, yq] = meshgrid(round(xmin):round(xmax), round(ymin):round(ymax));
    
    % Interpolate 2D scattered data
    pos_tab = griddata(x_coords, y_coords, T_vals, xq, yq, "nearest");
    pos_tab = rot90(pos_tab, 2);
    
    % Plot image
    % imagesc(pos_tab); axis image off
    
    % Dimensions of battery body and tabs
    [b_body_l, b_body_w] = size(b_body);
    [b_tab_l, b_tab_w] = size(pos_tab);
    
    % Create image of battery
    b_full = T_amb * ones(b_body_l + b_tab_l, b_body_w);
    b_full((b_tab_l + 1):end, :) = b_body;
    b_full(1:b_tab_l, (b_tab_off + 1):(b_tab_w + b_tab_off)) = neg_tab;
    b_full(1:b_tab_l, (b_body_w - b_tab_w - b_tab_off + 1):(end - b_tab_off)) = pos_tab;
    
    % Show image
    % imagesc(b_full); axis off image; colobar
    
    % Pad image
    % b_full = padarray(b_full, [1 1], T_amb, 'both');

    % Surface temperature array
    surface_t_array = b_full;