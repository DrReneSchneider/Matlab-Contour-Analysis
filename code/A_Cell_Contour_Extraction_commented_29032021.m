%% Clear workspace
    clear all
    close all

%% Ask for plasma membrane (e.g. LTi6B-GFP marker) and microtubule (mCHERRY-TUA5) image or z-stack
    [FileName1 PathName1] = uigetfile('*.TIF', 'Select TIFF containing cell outlines');
    info = imfinfo([PathName1 '/' FileName1]);
    num_images = numel(info);
    [FileName2 PathName2] = uigetfile('*.TIF', 'Select TIFF containing microtubules');
    info2 = imfinfo([PathName2 '/' FileName2]);
    num_images2 = numel(info2);
    
%% Import image/z-stack
    Image1_raw = zeros(info(1).Height, info(1).Width);
    Image2_raw = zeros(info(1).Height, info(1).Width);
    for k = 1:1:num_images
        I1 = imread([PathName1 FileName1], k);
        I2 = imread([PathName2 FileName2], k);
        Image1_raw(:,:,k)= I1;
        Image2_raw(:,:,k)= I2;
        k
    end
    [rows, columns, nSlices] = size(Image1_raw);

%% Use 16-bit scale
    I1 = round(im2double(I1)*2^16);
    I2 = round(im2double(I2)*2^16);

%% Calculate maximum intensity projection (MIP) of both images/z-stacks
    MIP_Image1 = max(Image1_raw, [],3);
    MIP_Image2 = max(Image2_raw, [],3);

%% Save plasma membrane and microtubule MIPs
    imwrite(uint16(MIP_Image1), [PathName1 '\A_' num2str(FileName1(1:end-4)) '_MIP.TIF']);
    imwrite(uint16(MIP_Image2), [PathName1 '\C_' num2str(FileName2(1:end-4)) '_MIP.TIF']);

%% Create watershed image
    L = watershed(Image1_raw);
    Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');

%% Save watersheded cell outlines (initial detection)
    imwrite(uint16(L), [PathName1 '\E_' num2str(FileName2(1:end-4)) '_WATERSHED_INIT.TIF']);

%% Figure B - Watershed figure with plasma membrane overlaid
    figure(3)
    imshow(MIP_Image1, [1*mean(mean(MIP_Image1)) 3*mean(mean(MIP_Image1))]);
    hold on
    himage = imshow(Lrgb);
    himage.AlphaData = 0.33;
    title('Wathershed superimposed on plasma membrane signal')
    saveas(gcf,[PathName1 '\B_' num2str(FileName1(1:end-4)) '_Watershed'], 'jpeg');
    close figure 3

%% Figure D - Watershed figure with microtubules overlaid
    figure(4)
    imagesc(MIP_Image2);
    hold on
    himage = imshow(Lrgb);
    himage.AlphaData = 0.33;
    title('Wathershed superimposed on microtubule signal')
    saveas(gcf,[PathName1 '\D_' num2str(FileName2(1:end-4)) '_Watershed'], 'jpeg');
    close figure 4

%% Remove falsely detected borders in initial watershedded cell outlines
    imagesc(L)
    msgbox('Use zoom to navigate in- and out of image. Select two cells that should be merged! Do one pair at a time','Advice!');
    L_new_1 = L;
    button = 1;
    while button == 1;
        imagesc(L_new_1)
        zoom on;
        pause() % you can zoom with your mouse and when your image is okay, you press any key
        zoom off; % to escape the zoom mode
        [ex_y ex_x] = ginput(2);
        cell_lyse_1 = L_new_1(round(ex_x(1,1)), round(ex_y(1,1)));
        cell_fix_2 = L_new_1(round(ex_x(2,1)), round(ex_y(2,1)));
        combine_1 = L_new_1 == cell_lyse_1;
        combine_2 = L_new_1 == cell_fix_2;
        combine = logical(combine_1 + combine_2);

        %% Merge two neighboring binary areas (cells)
        se = strel('disk',3);
        combine_new = imclose(combine, se);
        new = [];
        new = double(combine_new) * double(max([cell_lyse_1 cell_fix_2]));
        imagesc(new)
        Merge = zeros(rows, columns, 2);
        Merge(:,:,1) = L_new_1;
        Merge(:,:,2) = new;
        L_new_1 = max(Merge, [], 3);

        % User input: How to go on?
        choice = questdlg('How to proceeed?', 'Merge neighboring areas/cells', 'Fix more outlines', 'End this nightmare, please!', 'Fix more outlines');
        switch choice
        case 'Fix more outlines'
        button = 1;
        case 'End this nightmare, please!'
        button = 0; 
        end
        close all
    end
    close all

%% Reject cell from analysis (stomata, dividing cells, partly captured cells, etc.)
    imagesc(L_new_1)
    msgbox('Delete individual cells if they should not be included in the analysis. Click the cell you want to destroy. Leave stomata in the dataset for now, we will deal with stomata later.','Advice!');
    L_new_2 = L_new_1;
    button = 2;
    while button == 2;
        imagesc(L_new_2)
        % User input: x-y-coordinates of rejected cells
        [ex_y ex_x] = ginput(1);
        cell_lyse_1 = L_new_2(round(ex_x(1,1)), round(ex_y(1,1)));
        lyse_1 = L_new_2 == cell_lyse_1;

        % Remove selected cell (set to zero)
        new = [];
        new = (-1)*double(lyse_1)*cell_lyse_1;
        L_new_2 = L_new_2 + new;

        % User input: How to go on?
        choice = questdlg('How to proceeed?', 'Delete/lyse cells', 'Kill/lyse/destroy more cells!!!', 'Please! I had enough!', 'Fix more outlines');
        switch choice
        case 'Lyse more cells'
            button = 2;
        case 'Please! I had enough!'
            button = 0; 
        end
        close all
    end
    close all

%% Save 'cleaned-up' watersheded cell contours
    imwrite(uint16(L_new_2), [PathName1 '\F_' num2str(FileName2(1:end-4)) '_WATERSHED_CLEAN.TIF']);

%% Display cleaned watershed image
    figure(1)
    imagesc(L_new_2)
    hold on
    title('Wathershed image');

%% Select the cells you want to analyze and export cell outlines into excel file (x-,y- coordinates, cell number, etc.)
    cell_outlines = [];
    data = [];
    m=1;
    L_new_3 = L_new_2;
    imagesc(L_new_3)
    hold on
    msgbox('Select the cells you want to analyze. Now really only take the cells that you can track from 0 to 96h throughout','Advice!')
    button = 3;
    while button == 3
        % User input: x-y-coordinates of selected cells
        [take_y take_x] = ginput(1);
        cell_take_1 = L_new_3(round(take_x(1,1)), round(take_y(1,1)));
        take_1 = L_new_3 == cell_take_1;
        cell_i = cell_take_1;
        cell_number = num2str(m);
        contour = L_new_3 == cell_take_1;
        contour_0 = bwboundaries(bwperim(contour));
        cell_xpos = contour_0{1}(:,2)';
        cell_ypos = contour_0{1}(:,1)';
        cell_outlines{m} = {cell_i cell_number cell_xpos cell_ypos};

        % Plot cell number onto image
        CoM = regionprops(contour, 'Centroid');
        str = num2str(m);
        text(CoM.Centroid(1,1), CoM.Centroid(1,2), str, 'Color', 'red', 'FontSize', 10, 'HorizontalAlignment', 'center')
        hold on

        % User Input: How to continue?
        choice = questdlg('How to proceeed?', 'Select cells', 'Select more cells.', 'I am done', 'Select more cells.');
        switch choice
        case 'Select more cells.'
        button = 3;
        case 'I am done'
        button = 0; 
        end
        m = m+1;

    end

%% Save image of the selected cells
    saveas(gcf,[PathName1 'G_Selected Cells_' num2str(FileName2(1:end-4))], 'tiff');    
    close figure 1

%% Select stomata so that cell borders in contact with them can be rejected 
    L_new_4 = L_new_3;
    imagesc(L_new_4)
    hold on

    % User input: x-y-coordinates of stomata
        msgbox('Please select all stomata. Walls in contact with stomata will be ignored during the analysis. End selection by pressing ENTER.','Advice!')
        [stomata_y stomata_x] = ginput();
        if isempty(stomata_y) == 1;
            STOMATA_ALL = zeros(rows, columns);
        else
            STOMATA_ALL = zeros(rows, columns);
            cell_stomata = diag(L_new_4(round(stomata_x(:,1)), round(stomata_y(:,1))));
            for i = 1:1:length(cell_stomata)
                stomata = L_new_4 == cell_stomata(i,1);
                STOMATA_ALL = STOMATA_ALL + stomata; 
            end
        end
        close figure 1

%% Print selected cell outlines to excel file
    t_str = ['Cell contours for ' FileName1(8:end-8)];
    tab = char(9);
    newline = char(10);
    file_id = fopen([PathName1 'Cell_Contours_' num2str(m-1) '_Cells_' FileName1(8:end-9) '.txt'], 'w');
    strA = ['Cell number' tab 'cell x-position [pixels]' newline];
    strB = ['Watershed index' tab 'cell y-position [pixels]' newline];
    str = [t_str newline strA strB newline];
    fprintf(file_id, str);
    for p = 1:length(cell_outlines)
        foo_xpos = num2str(cell_outlines{1,p}{1,3},4);
        foo2_xpos = strsplit(foo_xpos);
        foo3_xpos = strjoin(foo2_xpos, '\t');
        str1 = [num2str(cell_outlines{1,p}{1,2},0) tab foo3_xpos newline];
        foo_ypos = num2str(cell_outlines{1,p}{1,4},4);
        foo2_ypos = strsplit(foo_ypos);
        foo3_ypos = strjoin(foo2_ypos, '\t');
        str2 = [num2str(cell_outlines{1,p}{1,1}) tab foo3_ypos newline];
        str_loop = [str1 str2 newline];
        fprintf(file_id, str_loop);
    end
    fclose(file_id);

%% Save workspace for further analyses ([...]_Workspace_Contours.mat)
    save([PathName1 'workspace_Contours.mat']);