%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  This script requires the following matlab scripts:                                                          %%
%%  densityplot.m (https://de.mathworks.com/matlabcentral/fileexchange/65166-densityplot-x-y-varargin)          %%
%%  circfit.m (https://de.mathworks.com/matlabcentral/fileexchange/5557-circle-fit)                             %%
%%  LineNormals2D.m (https://de.mathworks.com/matlabcentral/fileexchange/32696-2d-line-curvature-and-normals)   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear workspace
clear all
close all

%% Ask for location of workspace.mat file (generated by the A_Cell_Contour_Extraction.m script)
[FileName_WS PathName_WS] = uigetfile('*.MAT', 'Select workspace_Contours file');
PathName3 = PathName_WS;
mkdir([PathName3 'SingleCell-Data\']);
mkdir([PathName3 'Ensemble-Data\']);
load([PathName_WS FileName_WS])

%% Define the size of the contour margin (in pixels)
MARGIN_START = 1;                           %margin begins 1 pixel from the contour (towards inside)
MARGIN_END = 8;                             %margin ends 8 pixels from the contour (towards inside)
MARGIN = MARGIN_END - MARGIN_START +1;      %total depth of margin

%% Generate a mask to exclude walls around stomata
se = strel('disk', MARGIN_END+1);
STOMATA_ALL_DILATED = imdilate(STOMATA_ALL, se);
MIP_Image_Opt_2 = MIP_Image2 .* imcomplement(STOMATA_ALL_DILATED);

%% Go through all cells (index i) and collect data along all contour pixels (index j)
data_single = {};
data_multi = {};
for i = 1:length(cell_outlines);
    
    %% Define the cell's outline and measure geometric properties
        cell = L_new_2(round(mean(cell_outlines{1,i}{1,4})), round(mean(cell_outlines{1,i}{1,3})));
        contour = L_new_2 == cell;
        contour_dist = bwdist(~contour);
        contour_0 = bwboundaries(bwperim(bwmorph(contour, 'spur', Inf)));
        cell_props = regionprops(contour, 'Area', 'Perimeter', 'ConvexImage');
        cell_Area = cell_props.Area;                                        %Measure Area
        cell_Perimeter = cell_props.Perimeter;                              %Measure Perimeter
        cell_ConvexHull = regionprops(cell_props.ConvexImage, 'Perimeter'); %Measure Convex Hull
        cell_Lobyness = cell_Perimeter/cell_ConvexHull.Perimeter;           %Calculate Lobeyness (Sappalla et al., eLife, 2018)

    % Create an intensity profile along contour (obtain coordinates cx_0, cy_0)
        [cx_0, cy_0, I_signal_0] = improfile(MIP_Image_Opt_2, contour_0{1}(:,2), contour_0{1}(:,1));
    
    % Generate smoothed cell outline using smooth() function; stitch three
    % full circumnavigations of the cell together to avoid end-effects
    % (take the middle of the three)
        cx_sgfilt = [cx_0(1:end-1)' cx_0(1:end-1)' cx_0(1:end-1)']';
        cx_smooth = smooth(cx_sgfilt);
        cy_sgfilt = [cy_0(1:end-1)' cy_0(1:end-1)' cy_0(1:end-1)']';
        cy_smooth = smooth(cy_sgfilt);

    % Measure curvature & normal vectors along smoothed contour
        Vertices = [cx_smooth cy_smooth];
        Normals = LineNormals2D(Vertices);
    
    %% Create operators for orientation and anisotropy measurements (this
    % part is based on the Fiji-script from FibrilTool by Arezki Boudaoud
    % et al., Nature Protocols, 2014)
        imgI = MIP_Image_Opt_2;
        thresh = 2;

    % Compute gradient in x- and y-direction
        x = imtranslate(imgI, [-0.5, 0], 'cubic');
        x1 = imtranslate(x, [1, 0]);
        X = imsubtract(x, x1);
        y = imtranslate(imgI, [0, -0.5], 'cubic');
        y1 = imtranslate(y, [0, 1]);
        Y = imsubtract(y, y1);

    % Compute norm of gradient in g
        g = immultiply(X,X);
        gp = immultiply(Y,Y);
        G = sqrt(imadd(g, gp));
        % set the effect of the gradient to 1/255 when too low ; threshold = thresh
        G(G<=thresh) = 255;

    % Normalize "x" and "y" to components of normal
        x_norm = imdivide(X,G);
        y_norm = imdivide(Y,G);

    % Compute nxx, nxy, nyy
        nxx = immultiply(x_norm,x_norm);
        nxy = immultiply(x_norm,y_norm);
        nyy = immultiply(y_norm,y_norm);
    
    %% Go along contour and fit a circle (to get the local radius/curvature)
    % Use a subset of the contour to define local environment (here: window = MARGIN)
        Curv = [];
        CrossProd = [];
        Window = MARGIN;
        Anisotropy = [];
        Orientation = [];
        localLEC_list = [];
        Aniso_Norm = [];
        th = 0:pi/50:2*pi;
    
        for index = length(cy_0):1:2*length(cy_0)-1
            sample_x = [];
            sample_y = [];
            inner_sample_x = [];
            inner_sample_y = [];
    
            %% Curvature determination
            %get the query positions (x,y) along the contour
            %include pixels before and after the focus position (MARGIN/2)
            %for improved circle fitting
                sample_x = cx_smooth(index-floor(Window/2):index+floor(Window/2));
                sample_y = cy_smooth(index-floor(Window/2):index+floor(Window/2));
                
            % Option to have fitted circle plotted (for checking the functionality of the script)
%                 plot(cx_sgfilt, cy_sgfilt, 'b')
%                 hold on
%                 plot(sample_x, sample_y, 'bo');
%                 hold on
                [xc, yc, R] = circfit(sample_x, sample_y);
%                 plot(xc, yc, 'rx')
%                 hold on
%                 th = 0:pi/50:2*pi;
%                 plot(R*cos(th)+xc, R*sin(th)+yc, 'r-');
%                 (1/R);
                Curv = [Curv' (1/R)']';
                s = [cx_sgfilt(index+1,1)-cx_sgfilt(index,1) cy_sgfilt(index+1,1)-cy_sgfilt(index,1) 0];    %Vector component parallel to path
                k = [xc-cx_sgfilt(index,1) yc-cy_sgfilt(index,1) 0];                                        %Vector component perpendicular to path (in image plane)  
                SP = cross(s,k);                                                                            %Out-of-plane vector component (cross product)
                CrossProd = [CrossProd' SP']';
%                 axis([330 540 70 280])
%                 close

            % Option to use an additional separation between the cell contour and the analysis marging
                AnticlinalWallDistance = 2;
        
            % Define the outline of the analysis box     
                outer_sample_x = cx_smooth(index-floor(2*Window/2):index+floor(2*Window/2)) - (MARGIN_START-1+AnticlinalWallDistance)*Normals(index-floor(2*Window/2):index+floor(2*Window/2),1);
                outer_sample_y = cy_smooth(index-floor(2*Window/2):index+floor(2*Window/2)) - (MARGIN_START-1+AnticlinalWallDistance)*Normals(index-floor(2*Window/2):index+floor(2*Window/2),2);
                inner_sample_x = outer_sample_x - 2*(MARGIN+AnticlinalWallDistance)*Normals(index-floor(2*Window/2):index+floor(2*Window/2),1);
                inner_sample_y = outer_sample_y - 2*(MARGIN+AnticlinalWallDistance)*Normals(index-floor(2*Window/2):index+floor(2*Window/2),2);
                hull_x = [outer_sample_x' inner_sample_x' outer_sample_x(1)];
                hull_y = [outer_sample_y' inner_sample_y' outer_sample_y(1)];
                K = convhull(hull_x, hull_y);
                
            %% Anisotropy/Orientation Measurements (Fibril Tool)
            % Prepare Fibril Tool measurements in analysis box
                segment_x = round(hull_x(K));
                segment_y = round(hull_y(K));
                segment_ROI = roipoly(MIP_Image_Opt_2, segment_x, segment_y).*MIP_Image_Opt_2;
                %segment_ROI(segment_ROI==0) = NaN;
                
            % Option to display the analysis box for each query point       
%                 imagesc(segment_ROI)
%                 hold on
%                 plot(cx_smooth(length(cy_0):1:2*length(cy_0)-1,1), cy_smooth(length(cy_0):1:2*length(cy_0)-1,1));
%                 hold on
%                 plot([cx_smooth(length(cy_0):1:2*length(cy_0)-1,1) cx_smooth(length(cy_0):1:2*length(cy_0)-1,1)-MARGIN*Normals(length(cy_0):1:2*length(cy_0)-1,1)]', [cy_smooth(length(cy_0):1:2*length(cy_0)-1,1) cy_smooth(length(cy_0):1:2*length(cy_0)-1,1)-MARGIN*Normals(length(cy_0):1:2*length(cy_0)-1,2)]');
%                 hold on
%                 plot(sample_x, sample_y, 'ro')
%                 hold on
%                 plot(inner_sample_x, inner_sample_y, 'rx')
%                 axis equal
%                 hold on
%                 plot(hull_x(K), hull_y(K), 'g-')
%                 hold on
%                 axis([340 570 170 370])
%                 close
        
            % Fibril Tool : Crop part of nxx, nxy, and nyy within ROI
                segment_nxx = roipoly(nxx, segment_x, segment_y).*nxx;
                segment_nxx(segment_nxx==0) = NaN;
                segment_nxy = roipoly(nxy, segment_x, segment_y).*nxy;
                segment_nxy(segment_nxy==0) = NaN;
                segment_nyy = roipoly(nyy, segment_x, segment_y).*nyy;
                segment_nyy(segment_nyy==0) = NaN;
        
            % Fibril Tool : Measure the mean value of nxx, nxy, and nyy in ROI
                xx = nanmean(sort(segment_nxx(:)));
                xy = nanmean(sort(segment_nxy(:)));
                yy = nanmean(sort(segment_nyy(:)));

            % Fibril Tool: Eigenvalues and eigenvector of texture tensor
                m = (xx + yy)/2;
                d = (xx - yy)/2;
                v1 = m + sqrt(xy * xy + d * d);
                v2 = m - sqrt(xy * xy + d * d);

            % Fibril Tool : Direction/Orientation rel. to image frame
                tn_abs = - atan((v2 - xx) / xy);
                orient_abs = rad2deg(tn_abs);
        
            % Fibril Tool : Anisotropy
                anisotropy = abs((v1-v2) / 2 / m);
        
            % Fibril Tool : Direction/Orientation rel. to contour normal
                n_vector = [-Normals(index,1) -Normals(index,2) 0];
                orient_abs_vector = [anisotropy*cos(tn_abs) -anisotropy*sin(tn_abs) 0];
                tn_rel = dot(n_vector, orient_abs_vector)/(norm(n_vector)*norm(orient_abs_vector));
                orient_rel = acosd(tn_rel);
                aniso_norm = abs(dot(n_vector/norm(n_vector), orient_abs_vector/norm(orient_abs_vector)))*norm(orient_abs_vector);
        
            % Option to display Fibril Tool measurements in analysis box
            % for each query point
%                 imagesc(segment_ROI)
%                 colormap gray
%                 hold on
%                 plot(cx_smooth(length(cy_0):1:2*length(cy_0)-1,1), cy_smooth(length(cy_0):1:2*length(cy_0)-1,1));
%                 hold on
%                 plot([cx_sgfilt(index) cx_sgfilt(index)-MARGIN*Normals(index,1)], [cy_sgfilt(index) cy_sgfilt(index)-MARGIN*Normals(index,2)], 'r')
%                 hold on
%                 plot([cx_sgfilt(index) cx_sgfilt(index)+100*anisotropy*cos(tn_abs)], [cy_sgfilt(index) cy_sgfilt(index)-100*anisotropy*sin(tn_abs)], 'g') 
                segment_ROI(segment_ROI == 0) = NaN;    
%                 imagesc(segment_ROI)
%                 colormap gray
%                 hold on
%                 plot(cx_smooth(length(cy_0):1:2*length(cy_0)-1,1), cy_smooth(length(cy_0):1:2*length(cy_0)-1,1));
%                 hold on
%                 plot([cx_smooth(length(cy_0):1:2*length(cy_0)-1,1) cx_smooth(length(cy_0):1:2*length(cy_0)-1,1)-MARGIN*Normals(length(cy_0):1:2*length(cy_0)-1,1)]', [cy_smooth(length(cy_0):1:2*length(cy_0)-1,1) cy_smooth(length(cy_0):1:2*length(cy_0)-1,1)-MARGIN*Normals(length(cy_0):1:2*length(cy_0)-1,2)]');
%                 hold on
%                 plot(sample_x, sample_y, 'ro')
%                 hold on
%                 plot(inner_sample_x, inner_sample_y, 'rx')
%                 axis equal
%                 hold on
%                 plot(hull_x(K), hull_y(K), 'g-')
%                 hold on
%                 axis([cx_smooth(index)-30 cx_smooth(index)+30 cy_smooth(index)-30 cy_smooth(index)+30])
%                 axis ij
%                 plot([cx_smooth(index) cx_smooth(index)+100*anisotropy*cos(tn_abs)], [cy_smooth(index) cy_smooth(index)-100*anisotropy * sin(tn_abs)], '-', 'LineWidth', 5)
%                 hold off

            %% local LEC determination
                inside_list = [];
                inside = 0;
                R = 1;
                while inside < 2
                    xc = cx_smooth(index,1)-(R)*Normals(index,1);
                    yc = cy_smooth(index,1)-(R)*Normals(index,2);
                    % Option to plot fitted circle
%                     plot(xc, yc, 'rx')
%                     hold on
%                     plot(R*cos(th)+xc, R*sin(th)+yc, 'r-');
%                     hold on
                    inside = sum(inpolygon(cx_smooth, cy_smooth, R*cos(th)+xc, R*sin(th)+yc));
                    inside_list = [inside_list inside];
                    R = R+1;
                end
                localLEC_list = [localLEC_list length(inside_list)];

            % Collect data into smoothed vectors    
            Anisotropy = [Anisotropy anisotropy];
            Orientation = [Orientation tn_rel]; %tn_rel = 1 --> anisotropy || normal; tn_rel ~ 0 --> anisotropy perpendicular to normal
            Aniso_Norm = [Aniso_Norm aniso_norm]; %normal component of anisotropy (normal to contour)

        end
        
        %% Collect measurements into vectors and smooth data
        Curv2 = (CrossProd./abs(CrossProd)).*Curv;
        Curv2 = smooth(Curv2(:,3));
        Anisotropy2 = smooth(Anisotropy);
        Orientation2 = smooth(abs(Orientation));
        localLEC2 = smooth(localLEC_list);
        Aniso_Norm2 = smooth(Aniso_Norm);
        data_single{i} = [cx_0 cy_0 Curv2 Normals(length(cy_0):1:2*length(cy_0)-1,1) Normals(length(cy_0):1:2*length(cy_0)-1,2) Anisotropy2 Orientation2 localLEC2 Aniso_Norm2];
    
        %% Display data
        % Figure 4 - Print analyzed margin of cell{j} into projection image of MTs
            figure(4)
            MASK = ones(info(1).Height, info(1).Width);
            MASK_OUT = contour_dist > (MARGIN_START - 1);
            MASK_IN = contour_dist <= MARGIN-1;
            MASK = MASK_IN.*MASK_OUT;
            imagesc(MIP_Image_Opt_2.*MASK)
            axis equal
            axis off
            colormap gray
            imwrite(uint16(MIP_Image_Opt_2.*MASK), [PathName3 'SingleCell-Data\' 'H_Cell_' num2str(i) '_' num2str(FileName2(1:end-4)) '_Contours&MTs.TIF']);
            imwrite(uint16(MIP_Image_Opt_2.*MASK), [PathName3 'Ensemble-Data\' 'H_AllCells_' num2str(FileName2(1:end-4)) '_Contours&MTs.TIF'], 'WriteMode','append');
            close figure 4
    
        % Figure 5 - print mean intensity of margin around the contour
            figure(5)
            %Measure the intensities in a box perpendicular to the cell's outline
            %and depth defined by MARGIN_START and MARGIN_END
            mean_I_signal = [];
            for l = 1:1:length(cx_0)
                I_signal_normal = improfile(MIP_Image_Opt_2, [cx_0(l,1)-(MARGIN_START-1)*Normals(l,1) cx_0(l,1)-(MARGIN_END)*Normals(l,1)]', [cy_0(l,1)-(MARGIN_START-1)*Normals(l,2) cy_0(l,1)-(MARGIN_END)*Normals(l,2)]');
                if min(I_signal_normal) == 0
                    I_signal_normal = [];
                    I_signal_normal = 0;
                else
                    I_signal_normal = I_signal_normal(1:end-1);
                end
                mean_I_signal = [mean_I_signal mean(I_signal_normal(I_signal_normal>0))];
            end
    
        % Ignore walls in contact with a stomata
            STOMATA_INDEX = ~isnan(mean_I_signal);
            mean_I_signal = mean_I_signal .* STOMATA_INDEX;
            Curv_signal = Curv2;
            Curv_signal(STOMATA_INDEX == 0) = -1;
            Anisotropy_signal = Anisotropy2;
            Anisotropy_signal(STOMATA_INDEX == 0) = NaN;
            Orientation_signal = Orientation2;
            Orientation_signal(STOMATA_INDEX == 0) = NaN;
            localLEC_signal = localLEC2;
            localLEC_signal(STOMATA_INDEX == 0) = NaN;
            Aniso_Norm_signal = Aniso_Norm2;
            Aniso_Norm_signal(STOMATA_INDEX == 0) = NaN;
    
        % Prepare image output files
            meanMTplot = zeros(info(1).Height, info(1).Width);
            meanCurvplot = (-1)*ones(info(1).Height, info(1).Width);
            meanAnisotropyplot = (-1)*ones(info(1).Height, info(1).Width);
            meanOrientationplot = (-1)*ones(info(1).Height, info(1).Width);
            meanlocalLECplot = zeros(info(1).Height, info(1).Width);
            Aniso_Normplot = (-1)*ones(info(1).Height, info(1).Width);
            for l = 1:1:length(cx_0)
                meanMTplot(round(cy_0(l,1)):round(cy_0(l,1))+1, round(cx_0(l,1)):round(cx_0(l,1))+1) = mean_I_signal(l);
                meanCurvplot(round(cy_0(l,1)):round(cy_0(l,1))+1, round(cx_0(l,1)):round(cx_0(l,1))+1) = Curv_signal(l);
                meanAnisotropyplot(round(cy_0(l,1)):round(cy_0(l,1))+1, round(cx_0(l,1)):round(cx_0(l,1))+1) = Anisotropy_signal(l);
                meanOrientationplot(round(cy_0(l,1)):round(cy_0(l,1))+1, round(cx_0(l,1)):round(cx_0(l,1))+1) = Orientation_signal(l);
                meanlocalLECplot(round(cy_0(l,1)):round(cy_0(l,1))+1, round(cx_0(l,1)):round(cx_0(l,1))+1) = localLEC_signal(l);
                Aniso_Normplot(round(cy_0(l,1)):round(cy_0(l,1))+1, round(cx_0(l,1)):round(cx_0(l,1))+1) = Aniso_Norm_signal(l);
            end
            clims = [min(mean_I_signal) max(mean_I_signal)];
            imagesc(uint16(meanMTplot), clims);
            axis equal
            axis off
            colormap gray
            hold off
            imwrite(uint16(meanMTplot), [PathName3 'SingleCell-Data\I_Cell_' num2str(i) '_' num2str(FileName2(1:end-4)) '_Contours_meanIntensity.TIF']);
            imwrite(uint16(meanMTplot), [PathName3 'Ensemble-Data\I_AllCells_' num2str(FileName2(1:end-4)) '_Contours_meanIntensity.TIF'], 'WriteMode', 'append');
            close figure 5
    
        % Option to display the normal vectors to the cell outline
%             plot(cx_0(:,1), cy_0(:,1));
%             hold on
%             plot([cx_0(:,1) cx_0(:,1)-MARGIN*Normals(:,1)]', [cy_0(:,1) cy_0(:,1)-MARGIN*Normals(:,2)]');
    
        % Figure 6 - print curvature along contour
            figure(6)
            clims = [min(Curv_signal) max(Curv_signal)];
            imagesc(meanCurvplot, clims);
            axis equal
            axis off
            colormap jet
            hold off
            imwrite(uint16((meanCurvplot+1)*1000), [PathName3 'SingleCell-Data\J_Cell_' num2str(i) '_' num2str(FileName2(1:end-4)) '_Contours_meanCurvature.TIF']);
            imwrite(uint16((meanCurvplot+1)*1000), [PathName3 'Ensemble-Data\J_AllCells_' num2str(FileName2(1:end-4)) '_Contours_meanCurvature.TIF'], 'WriteMode', 'append');
            close figure 6
            
        % Measure the correlation coefficient between mean_I_signal and curvature
            [CorrCoefs, CorrPs, CorrLoCIs, CorrHiCIs] = corrcoef(mean_I_signal(isfinite(mean_I_signal)), Curv_signal(isfinite(mean_I_signal)));
            CorrCoef = CorrCoefs(1,2);
            CorrP = CorrPs(1,2);
            CorrLoCI = CorrLoCIs(1,2);
            CorrHiCI = CorrHiCIs(1,2);
    
        % Figure 7 - print anisotropy along contour
            figure(7)
            clims = [min(Anisotropy_signal) max(Anisotropy_signal)];
            imagesc(meanAnisotropyplot, clims);
            axis equal
            axis off
            colormap jet
            hold off
            imwrite(uint16((meanAnisotropyplot+1)*1000), [PathName3 'SingleCell-Data\K_Cell_' num2str(i) '_' num2str(FileName2(1:end-4)) '_Contours_meanAnisotropy.TIF']);
            imwrite(uint16((meanAnisotropyplot+1)*1000), [PathName3 'Ensemble-Data\K_AllCells_' num2str(FileName2(1:end-4)) '_Contours_meanAnisotropy.TIF'], 'WriteMode', 'append');
            close figure 7
    
        % Figure 8 - print orientation along contour
            figure(8)
            clims = [min(Orientation_signal) max(Orientation_signal)];
            imagesc(meanOrientationplot, clims);
            axis equal
            axis off
            colormap jet
            hold off
            imwrite(uint16((meanOrientationplot+1)*1000), [PathName3 'SingleCell-Data\L_Cell_' num2str(i) '_' num2str(FileName2(1:end-4)) '_Contours_meanOrientation.TIF']);
            imwrite(uint16((meanOrientationplot+1)*1000), [PathName3 'Ensemble-Data\L_AllCells_' num2str(FileName2(1:end-4)) '_Contours_meanOrientation.TIF'], 'WriteMode', 'append');
            close figure 8
    
        % Figure 9 - print local LEC along contour
            figure(9)
            clims = [min(localLEC_signal) max(localLEC_signal)];
            imagesc(meanlocalLECplot, clims);
            axis equal
            axis off
            colormap jet
            hold off
            imwrite(uint16((meanlocalLECplot)), [PathName3 'SingleCell-Data\N_Cell_' num2str(i) '_' num2str(FileName2(1:end-4)) '_Contours_localLEC.TIF']);
            imwrite(uint16((meanlocalLECplot)), [PathName3 'Ensemble-Data\N_AllCells_' num2str(FileName2(1:end-4)) '_Contours_localLEC.TIF'], 'WriteMode', 'append');
            close figure 9
            
         % Figure 10 - print normal component of anisotropy tensor
            figure(10)
            clims = [min(Aniso_Norm_signal) max(Aniso_Norm_signal)];
            imagesc(Aniso_Normplot, clims);
            axis equal
            axis off
            colormap jet
            hold off
            imwrite(uint16((Aniso_Normplot+1)*1000), [PathName3 'SingleCell-Data\O_Cell_' num2str(i) '_' num2str(FileName2(1:end-4)) '_Contours_Aniso_Norm.TIF']);
            imwrite(uint16((Aniso_Normplot+1)*1000), [PathName3 'Ensemble-Data\O_AllCells_' num2str(FileName2(1:end-4)) '_Contours_Aniso_Norm.TIF'], 'WriteMode', 'append');
            close figure 10  
         
        % Figure showing all individual cells (histograms and I-vs-Curvature)
            figure(11)
            subplot(2,2,1)
                bins = min(Curv_signal):0.02:max(Curv_signal);
                H = hist(Curv_signal, bins);
                hist(Curv_signal, bins)
                axis([-0.4 0.6 0 2*max(H)])
                xlabel('Curvature [1/px]')
                ylabel('frequency')
                hold off
            subplot(2,2,3)
                bins = min(mean_I_signal):(max(mean_I_signal)-min(mean_I_signal))/30:max(mean_I_signal);
                H = hist(mean_I_signal, bins);
                hist(mean_I_signal, bins)
                axis([0 max(mean_I_signal) 0 2*max(H)])
                xlabel('Mean intensity [a.u.]')
                ylabel('frequency')
                hold off
            subplot(2,2,[2 4])
                binsize = 0.02;
                bins = -0.4:binsize:0.6;

        % Fit regression line to I-vs-Curvature plot
            data_binning_Mean = [];
            data_binning_Std = [];
            data_binning_N = [];
            mean_I_signal_norm = (mean_I_signal-min(mean_I_signal))/(max(mean_I_signal)-min(mean_I_signal));
            for ii=1:length(bins)
                curv_bin_1 = Curv_signal <= bins(1,ii)+binsize/2;
                curv_bin_2 = Curv_signal > bins(1,ii)-binsize/2;
                curv_bin = logical(curv_bin_1.*curv_bin_2);

                Signal_bin_mean = mean(mean_I_signal_norm(curv_bin));
                Signal_bin_std = std(mean_I_signal_norm(curv_bin));
                Signal_bin_n = sum(curv_bin);

                data_binning_Mean = [data_binning_Mean Signal_bin_mean];
                data_binning_Std = [data_binning_Std Signal_bin_std];
                data_binning_N = [data_binning_N Signal_bin_n];
            end
            Ignore_data = ~isnan(data_binning_Mean);
            [Fit3, S3] = polyfit(bins(Ignore_data), data_binning_Mean(Ignore_data), 1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
            [f3, delta3] = polyval(Fit3, bins, S3);
            errorbar(bins, data_binning_Mean, data_binning_Std./sqrt(data_binning_N), 'bo');
            hold on;
            errorbar(bins, f3, delta3, 'r');
            hold on
            xlabel('curvature [1/px]');
            ylabel('norm. intensity [a.u.]');
            axis([-0.4 0.6 0 1]);
            hold off
            % Determine 95% confidence intervals for fit
            x_regress = [ones(size(Curv2)) Curv2];
            [beta, beta_int] = regress(mean_I_signal_norm', x_regress);
            saveas(gcf,[PathName3 'SingleCell-Data\M_' num2str(i) '_' num2str(FileName2(1:end-4)) '_Curv_vs_Intensity_EachCell'], 'jpg');
            close figure 11
    
        %% Write single-cell data sets into spreadsheets to be analyzed later
            data_singlecell{i} = [i cell_Area cell_Perimeter cell_Lobyness beta(2,1) (beta_int(2,2)-beta_int(2,1))/4 mean(Curv2) std(Curv2) mean(Anisotropy2) std(Anisotropy2) mean(Orientation2) std(Orientation2) max(localLEC2) median(localLEC2) std(localLEC2) mean(Aniso_Norm2) std(Aniso_Norm2) CorrCoef CorrP CorrLoCI CorrHiCI];
            data_multi{i} = [cx_0(STOMATA_INDEX) cy_0(STOMATA_INDEX) mean_I_signal(STOMATA_INDEX)' Curv_signal(STOMATA_INDEX) Anisotropy_signal(STOMATA_INDEX) Orientation_signal(STOMATA_INDEX) ];

        %% Print single-cell measurements to excel file
            t_str ='SingleCell_Measurements';
            tab = char(9);
            newline = char(10);
            file_id = fopen([PathName3 'SingleCell-Data\SingleCell_MEASUREMENTS_' num2str(length(cell_outlines)) ' cells -- ' FileName1(8:end-9)  '.txt'],'w');
            str2 = ['Cell number' tab 'Cell area [px^2]' tab 'Cell perimeter [px]' tab 'Lobyness' tab 'Curv-vs-Intensity Correlation (Slope)' tab 'Curv-vs-Intensity [+/- fit error (std)]' tab 'Mean Curvature [1/px]' tab 'Std Curvature' tab 'Mean Anisotropy' tab 'Std Anisotropy' tab 'Mean Orientation (1 = || to Normal)' tab 'Std Orientation' tab 'LEC' tab 'Median local-LEC' tab 'Std local-LEC' tab 'Mean Anisotropy normal to wall' tab 'Std Anisotropy normal' tab 'Correlation Coefficient' tab 'Correlation p-value' tab 'Correlation lower .95 CI' tab 'Correlation upper .95 CI' newline];
            str = [t_str newline str2];
            fprintf(file_id,str);
            for p = 1:length(data_singlecell)
                str = [num2str(data_singlecell{1,p}(1,1),5) tab num2str(data_singlecell{1,p}(1,2),5) tab num2str(data_singlecell{1,p}(1,3),5) tab num2str(data_singlecell{1,p}(1,4),5) tab num2str(data_singlecell{1,p}(1,5),5) tab num2str(data_singlecell{1,p}(1,6),5) tab num2str(data_singlecell{1,p}(1,7),5) tab num2str(data_singlecell{1,p}(1,8),5) tab num2str(data_singlecell{1,p}(1,9),5) tab num2str(data_singlecell{1,p}(1,10),5) tab num2str(data_singlecell{1,p}(1,11),5) tab num2str(data_singlecell{1,p}(1,12),5) tab num2str(data_singlecell{1,p}(1,13),5) tab num2str(data_singlecell{1,p}(1,14),5) tab num2str(data_singlecell{1,p}(1,15),5) tab num2str(data_singlecell{1,p}(1,16),5) tab num2str(data_singlecell{1,p}(1,17),5) tab num2str(data_singlecell{1,p}(1,18),5) tab num2str(data_singlecell{1,p}(1,19),5) tab num2str(data_singlecell{1,p}(1,20),5) tab num2str(data_singlecell{1,p}(1,21),5) newline];
                fprintf(file_id,str);
            end
            fclose(file_id);
end

    %% Global analysis (all cells pooled together)
        mean_I_signal_multi = [];
        mean_Curv_data_multi = [];
        mean_Aniso_Signal_multi = [];
        mean_Orient_Signal_multi = [];
    % Scatter plot
        figure(12)
        for tt = 1:length(data_multi)
            MT_Signal = data_multi{1,tt}(:,3);
            Anisotropy_Signal = data_multi{1,tt}(:,5);
            Orientation_Signal = data_multi{1,tt}(:,6);
            MT_Signal_norm = (MT_Signal-min(MT_Signal))/(max(MT_Signal)-min(MT_Signal));
            Curv_Signal = data_multi{1,tt}(:,4);
            mean_I_signal_multi = [mean_I_signal_multi' MT_Signal_norm']';
            mean_Curv_data_multi = [mean_Curv_data_multi' Curv_Signal']';
            mean_Aniso_Signal_multi = [mean_Aniso_Signal_multi' Anisotropy_Signal']';
            mean_Orient_Signal_multi = [mean_Orient_Signal_multi' Orientation_Signal']';
            plot(Curv_Signal, MT_Signal_norm, 'o');
            hold on
        end
        max_Curv = max(mean_Curv_data_multi);
        min_Curv = min(mean_Curv_data_multi);
        mean_Curv = mean(mean_Curv_data_multi);
        std_Curv = std(mean_Curv_data_multi);
        max_Signal = max(mean_I_signal_multi);
        min_Signal = min(mean_I_signal_multi);
        mean_Signal = mean(mean_I_signal_multi);
        std_Signal = std(mean_I_signal_multi);
        plot(mean_Curv, mean_Signal, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r')
        hold on

        Fit = polyfit(mean_Curv_data_multi, mean_I_signal_multi, 1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
        f = polyval(Fit,mean_Curv_data_multi);
        plot(mean_Curv_data_multi, f, 'r')
        xlim([-0.4 1])
        ylim([0 1])
        %legend('data', 'center point', 'linear fit')
        text(0.9*min(mean_Curv_data_multi), 0.9*max(mean_I_signal_multi), ['y = mx + n; m = ' num2str(round(Fit(1,1),2)) '; n = ' num2str(round(Fit(1,2),2))])
        %T = table(mean_Curv_data_multi, mean_I_signal_multi, Fit, mean_I_signal_multi-Fit,'VariableNames',{'Curvature','MT intensity','Fit','FitError'})
        saveas(gcf,[PathName3 'N_ALL_' num2str(FileName2(1:end-4)) '_Curv_Intensity_AllCells'], 'tiff');
        hold off
        close figure 12

        % Denisty plot
        densityplot(mean_Curv_data_multi, mean_I_signal_multi)
        hold on
        plot(mean_Curv, mean_Signal, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r')
        colorbar
        title('Heatmap - Curvature vs. normalized MT intensity')
        xlabel('curvature [necks (concave) --> lobes (convex)]')
        ylabel('norm. intensity [a.u.]')
        saveas(gcf,[PathName3 'L_ALL_' num2str(FileName2(1:end-4)) '_C-I-Heatmap'], 'tiff');
        close
        binsize = 0.02;
        bins = -0.4:binsize:1;
        data_binning_Mean = [];
        data_binning_Std = [];
        data_binning_N = [];

        for ii=1:length(bins)
            curv_bin_1 = mean_Curv_data_multi <= bins(1,ii)+binsize/2;
            curv_bin_2 = mean_Curv_data_multi > bins(1,ii)-binsize/2;

            curv_bin = logical(curv_bin_1.*curv_bin_2);

            Signal_bin_mean = mean(mean_I_signal_multi(curv_bin));
            Signal_bin_std = std(mean_I_signal_multi(curv_bin));
            Signal_bin_n = sum(curv_bin);

            data_binning_Mean = [data_binning_Mean Signal_bin_mean];
            data_binning_Std = [data_binning_Std Signal_bin_std];
            data_binning_N = [data_binning_N Signal_bin_n];
        end

    %% Print data to excel file
        t_str ='CellContourAnalysis';
        tab = char(9);
        newline = char(10);
        file_id2 = fopen([PathName3 'Cell_Contour_ANALYSIS_' num2str(length(cell_outlines)) ' cells -- ' FileName1(8:end-8)  '.txt'],'w');
        str2 = ['Cell outline curvature' tab 'Mean signal intensity [a.u.]' tab 'Standard deviation of signal intensity [a.u.]' tab 'Number of data points' tab 'Standard error of the mean (S.E.M.)' newline];
        str = [t_str newline str2];
        fprintf(file_id2,str);
        for p = 1:length(data_binning_Mean)
            str = [num2str(bins(1,p),3) tab num2str(data_binning_Mean(1,p),4) tab num2str(data_binning_Std(1,p),4) tab num2str(data_binning_N(1,p),4) tab num2str(data_binning_Std(1,p)/sqrt(data_binning_N(1,p)),4) newline];
            fprintf(file_id2,str);
        end
        fclose(file_id2);

        %% Save Workspace ([...]_Analysis.mat)
        save([PathName3 '\workspace_Analysis.mat']);