% function [avg_rgb_flat, avg_lab_flat, exlab_flat] = BiologPlateReader(image_name, numclus, name)
%%%% Plate Reader for Biolog Assays for C112L Lab Project%%%%%
% willsharpless@berkeley.edu

% Inputs:
% image_name - the string name of plate image file 
% numclus - an integer for the number of clusters to make
% name - name of the plate being analyzed for figure titles

% Outputs:
% avg_rgb_flat - the average rgb within each well that imfindcircle() finds
% avg_lab_flat - avg_rgb_flat converted to CIE-l*a*b*
% exlab_flat - the experiment labels associated with avg_rgb_flat and avg_lab_flat

% 1. You will be prompted to crop the image to the size of the plate
% 2. You will need to select the four corners of the corner wells in order
% to adjust the perspective distortion
% 3. That image will segmented such that each well is partitioned
% individually
% 4. imfindcircles() will locate the wells and then the color will be
% averaged
% 5. Finally the colors will be clustered by similarity and then ordered by
% darkness
name = 'Isolate 9';
image_name = 'iso9.jpg';
numclus = 8;


I = imread(image_name);
% I = rgb2lab(I);
% figure, imshow(I), title('Wild Type Plate')
K = imcrop(I);


% Perspective Fix  %%%%%%%%%%%%%%%%%%%%%%

% The following workflow for adjusting perspective-skew is from
% Michael Chan's "Perspective Control/Correction" file at 
% https://www.mathworks.com/matlabcentral/fileexchange/35531-perspective-control-correction
targetImageData = K;

% 5.1 Undoing Perspective Distortion of Planar Surface
% Note corners must be chosen clockwise or cc!

% a)
imshow(targetImageData);

% b)
fprintf('Corner selection must be clockwise or anti-clockwise.\n');
[X Y] = ginput(4);

%X = uint8(X);
%Y = uint8(Y);
[X Y] = sortPolyFromClockwiseStartingFromTopLeft( X, Y );

[sy sx z] = size(K);

x=[1;sx-100;sx-100;1];
y=[1;1;sy-100;sy-100];

% c)
A=zeros(8,8);
A(1,:)=[X(1),Y(1),1,0,0,0,-1*X(1)*x(1),-1*Y(1)*x(1)];
A(2,:)=[0,0,0,X(1),Y(1),1,-1*X(1)*y(1),-1*Y(1)*y(1)];

A(3,:)=[X(2),Y(2),1,0,0,0,-1*X(2)*x(2),-1*Y(2)*x(2)];
A(4,:)=[0,0,0,X(2),Y(2),1,-1*X(2)*y(2),-1*Y(2)*y(2)];

A(5,:)=[X(3),Y(3),1,0,0,0,-1*X(3)*x(3),-1*Y(3)*x(3)];
A(6,:)=[0,0,0,X(3),Y(3),1,-1*X(3)*y(3),-1*Y(3)*y(3)];

A(7,:)=[X(4),Y(4),1,0,0,0,-1*X(4)*x(4),-1*Y(4)*x(4)];
A(8,:)=[0,0,0,X(4),Y(4),1,-1*X(4)*y(4),-1*Y(4)*y(4)];

v=[x(1);y(1);x(2);y(2);x(3);y(3);x(4);y(4)];

u=A\v;
%which reshape

U=reshape([u;1],3,3)';

w=U*[X';Y';ones(1,4)];
w=w./(ones(3,1)*w(3,:));

% d)
%which maketform
T=maketform('projective',U');

%which imtransform
J=imtransform(targetImageData,T,'XData',[1 sx-100],'YData',[1 sy-100]);

% figure,
% imshow(J)
% title('Perspective Fixed')

% Partitioning Image  %%%%%%%%%%%%%%%%%%%%%%

[rows columns numberOfColorBands] = size(J);

close all

blockSizeR = floor(rows/8); % Rows in image/8 (partition for 96 well plate)
blockSizeC = floor(columns/12); % Columns in image/12 (partition for 96 well plate).

% Figure out the size of each block in rows.
% Most will be blockSizeR but there may be a remainder amount of less than that.
wholeBlockRows = floor(rows / blockSizeR);
blockVectorR = [blockSizeR * ones(1, wholeBlockRows), rem(rows, blockSizeR)];
% Figure out the size of each block in columns.
wholeBlockCols = floor(columns / blockSizeC);
blockVectorC = [blockSizeC * ones(1, wholeBlockCols), rem(columns, blockSizeC)];

% Each cell (except for the remainder cells at the end of the image)
% in the array contains a blockSizeR by blockSizeC by 3 color array.
% This line is where the image is actually divided up into blocks.
if numberOfColorBands > 1
    % It's a color image.
    ca = mat2cell(J, blockVectorR, blockVectorC, numberOfColorBands);
else
    ca = mat2cell(J, blockVectorR, blockVectorC);
end

% Display All the Chopped Images (For Debugging) %%%%%%%%%%%%%%%%%%%%%%

%%Uncomment to display segmented images
% plotIndex = 1;
% numPlotsR = size(ca, 1);
% numPlotsC = size(ca, 2);
% figure
% for r = 1 : numPlotsR
%     for c = 1 : numPlotsC
%         fprintf('plotindex = %d,   c=%d, r=%d\n', plotIndex, c, r);
%         % Specify the location for display of the image.
%         subplot(numPlotsR, numPlotsC, plotIndex);
%         % Extract the numerical array out of the cell
%         rgbBlock = ca{r,c};
%         imshow(rgbBlock); % Could call imshow(ca{r,c}) if you wanted to.
%         [rowsB columnsB numberOfColorBandsB] = size(rgbBlock);
%         % Make the caption the block number.
%         caption = sprintf('Block #%d of %d\n%d rows by %d columns', ...
%             plotIndex, numPlotsR*numPlotsC, rowsB, columnsB);
%         title(caption);
%         drawnow;
%         % Increment the subplot to the next location.
%         plotIndex = plotIndex + 1;
%     end
% end

% Store Well Color in Each Image  %%%%%%%%%%%%%%%%%%%%%%
ca_real = ca(1:8,1:12); %only the wells
numPlotsR = size(ca_real, 1);
numPlotsC = size(ca_real, 2);

avg_rgb = zeros(numPlotsR,numPlotsC,3);
plotIndex = 1;

lastcenter = zeros(1,2);
lastradius = zeros(1,1);
first_col_center = zeros(1,2);
first_col_radius = zeros(1,1);

figure('Position', [10 10 900 600])
set(gcf,'color','w')
for r = 1 : numPlotsR
    for c = 1 : numPlotsC
        subplot(numPlotsR, numPlotsC, plotIndex);
        rgbBlock = ca_real{r,c};
        imshow(rgbBlock);
        
         % Try to find the circle
         % Note these sensitivities have been tailored to work for my
         % images but different lighting or better images might require
         % less sensitivity.
        [centers, radii] = imfindcircles(rgbBlock,[23 27],'ObjectPolarity','dark','Sensitivity',0.9);
        
        % Search and Dial Up Sensitivity
        if isempty(centers)
            [centers, radii] = imfindcircles(rgbBlock,[23 27],'ObjectPolarity','dark','Sensitivity',0.92);
        end
        
        if isempty(centers)
            [centers, radii] = imfindcircles(rgbBlock,[23 27],'ObjectPolarity','dark','Sensitivity',0.94);
        end
        
        if isempty(centers)
            [centers, radii] = imfindcircles(rgbBlock,[23 27],'ObjectPolarity','dark','Sensitivity',0.96);
        end
        
        if isempty(centers)
            [centers, radii] = imfindcircles(rgbBlock,[23 27],'ObjectPolarity','dark','Sensitivity',0.98);
        end
        
        if isempty(centers) %cannot find circle (probably empty well)
            %first well is manually set to second well position
            if lastcenter == zeros(1,2)
                lastcenter = [38.1483 34.7987];
                lastradius = [25];
            end
            
            if c == 1 && r > 1
                lastcenter = first_col_center;
                lastradius = first_col_radius;
            end
            
            %we make well = position of last well       
            centers = lastcenter(1,:);
            radii = lastradius(1,:);
            disp('Warning: Cannot find circle, using the location of previously found well')
        end
        
%         if isempty(centers) (for Debugging)
%             %no circle found is average color of rgbBlock
%             disp('no circles')
%             red = rgbBlock(:,:,1); green = rgbBlock(:,:,2); blue = rgbBlock(:,:,3);
%             avg_red = mean(red(:)); avg_green = mean(green(:)); avg_blue = mean(blue(:));
%             avg_rgb(r,c,:) = [avg_red avg_green avg_blue];
%             title(sprintf('[%0.1f %0.1f %0.1f]',avg_rgb(r,c,1),avg_rgb(r,c,2),avg_rgb(r,c,3)))
%             
%             lastcenter = zeros(1,2);
            
        if numel(centers) == 2 %only 1 circle
            disp('found one circle')
            h = viscircles(centers,radii,'Color','g');
            
            %identifying pixels in circle
            [x,y] = meshgrid(1:size(rgbBlock,2),1:size(rgbBlock,1));
            mask = hypot(x - centers(1), y - centers(2)) < radii(1)/2;
            
            %isoalting color averages
            red = rgbBlock(:,:,1); green = rgbBlock(:,:,2); blue = rgbBlock(:,:,3);
            avg_red = mean(red(mask)); avg_green = mean(green(mask)); avg_blue = mean(blue(mask));
            avg_rgb(r,c,:) = [avg_red avg_green avg_blue];
%             title(sprintf('r=%0.1f s=%0.2f',radii(1)))
            
            lastcenter = centers(1,:);
            lastradius = radii(1);
            
            if c == 1
                first_col_center = centers(1,:);
                first_col_radius = radii(1);
            end
        else
            disp('found more than one circle')
%             h = viscircles(centers,radii,'Color','r');
            
            [x,y] = meshgrid(1:size(rgbBlock,2),1:size(rgbBlock,1));
            
            %use last well (assume dx,dy are neglible between wells in my segmentation)
            centers = lastcenter(1,:);
            radii = lastradius(1,:);
            h = viscircles(centers,radii,'Color','g'); %circle taken is green
            
            %identifying pixels in circle
            [x,y] = meshgrid(1:size(rgbBlock,2),1:size(rgbBlock,1));
            mask = hypot(x - centers(1), y - centers(1)) < radii(1)/2;
            
            %isoalting color averages
            red = rgbBlock(:,:,1); green = rgbBlock(:,:,2); blue = rgbBlock(:,:,3);
            avg_red = mean(red(mask)); avg_green = mean(green(mask)); avg_blue = mean(blue(mask));
            avg_rgb(r,c,:) = [avg_red avg_green avg_blue];
            title(sprintf('r=%0.1f s=%0.2f',radii(1)))
            
            lastcenter = centers(1,:);
            lastradius = radii(1);
            
            if c == 1
                first_col_center = centers(1,:);
                first_col_radius = radii(1);
            end
        end
        
         title(sprintf('[%0.1f %0.1f %0.1f]',avg_rgb(r,c,1),avg_rgb(r,c,2),avg_rgb(r,c,3)))
        
%        title([label_ca{r,c} sprintf(': %d', ps_idx_l(r,c))])
        drawnow;
        plotIndex = plotIndex + 1;
    end
end

sgtitle('\qquad  \textbf{\textit{Finding the Wells with imfindcircles() and Averaging Color}}', 'FontSize', 20,'Interpreter', 'Latex')

% Creating Labels and Lists %%%%%%%%%%%%%%%%%%%%%%

label_ca = cell(8,12);
well_letters = {'A','B','C','D','E','F','G','H'};
well_numbers = 12:-1:1; %backwards since I took plate photos from the underside for better resolution

for r = 1:8
    for c=1:12
        label_ca{r,c} = strcat(well_letters{r},string(well_numbers(c)));
    end
end

avg_rgb_flat = reshape(avg_rgb, [96, 3]);
label_ca_flat = reshape(label_ca, [96, 1]);

%%Uncomment to Plot RGB Data
% figure('Position', [10 10 900 600])
% set(gcf,'color','w')
% for i=1:96
%     plot3(avg_rgb_flat(i,1),avg_rgb_flat(i,2),avg_rgb_flat(i,3),'bo'); hold on
%     text(avg_rgb_flat(i,1)+1,avg_rgb_flat(i,2)+1,avg_rgb_flat(i,3)+1,label_ca_flat(i));
% end
% grid on
% title('RGB distribution of Wells')
% hold off

avg_lab_flat = rgb2lab(avg_rgb_flat); 
%%Uncomment to Plot CIE LAB Data
% convert RGB to CIE-L*a*b* because 
% "rgb isn't a good model for color distances"
% 
% figure('Position', [10 10 900 600])
% set(gcf,'color','w')
% for i=1:96
%     plot3(avg_lab_flat(i,1),avg_lab_flat(i,2),avg_lab_flat(i,3),'bo'); hold on
%     text(avg_lab_flat(i,1)+1,avg_lab_flat(i,2)+1,avg_lab_flat(i,3)+1,label_ca_flat(i));
% end
% grid on
% title('CIE-L*A*B* distribution of Wells')
% hold off

% Cluster the Wells by Color %%%%%%%%%%%%%%%%%%%%%%

[idx_r, Cr] = kmeans(avg_rgb_flat,numclus); %clustering by rgb (better if only presence of one color)
[idx_l, Cl] = kmeans(avg_lab_flat,numclus); %clustering by cie-lab (better because minor differences due to light/color behind plate)
colors=rand(numclus,3);

%%Uncomment to plot clustering
% figure('Position', [10 10 900 600])
% set(gcf,'color','w')
% for a=1:numclus
%     %disp(a) %progress report
%     plot3(avg_rgb_flat(idx_r==a,1),avg_rgb_flat(idx_r==a,2),avg_rgb_flat(idx_r==a,3),'o','color',colors(a,:),"LineWidth",2)
%     hold on
% end
% for i=1:96
%     text(avg_rgb_flat(i,1)+1,avg_rgb_flat(i,2)+1,avg_rgb_flat(i,3)+1,label_ca_flat(i));
% end
% grid on
% title(sprintf('RGB Clustering of Wells with %d Clusters',numclus));
% hold off
% 
figure('Position', [10 10 900 600])
set(gcf,'color','w')
for a=1:numclus
    %disp(a) %progress report
    plot3(avg_lab_flat(idx_l==a,1),avg_lab_flat(idx_l==a,2),avg_lab_flat(idx_l==a,3),'o','color',colors(a,:),"LineWidth",2)
    hold on
end
for i=1:96
    text(avg_lab_flat(i,1)+1,avg_lab_flat(i,2)+1,avg_lab_flat(i,3)+1,label_ca_flat(i));
end
grid on
title(sprintf('CIE LAB Clustering of Wells with %d Clusters',numclus));
hold off

ps_idx_r = reshape(idx_r, [8, 12]); %plate shape of clusters
ps_idx_l = reshape(idx_l, [8, 12]);

% Gap Statistic %%%%%%%%%%%%%%%%%%%%%%
% % Tibshirani, R., G. Walther, and T. Hastie. ?Estimating the number of clusters in a data set via the gap statistic.? Journal of the Royal Statistical Society: Series B. Vol. 63, Part 2, 2001, pp. 411?423.
% eva = evalclusters(avg_rgb_flat,'kmeans','gap','KList',[1:15]);
% plot(eva)

% Display All the Chopped Images w Cluster %%%%%%%%%%%%%%%%%%%%%%

% plotIndex = 1;
% 
% figure('Position', [10 10 900 600])
% set(gcf,'color','w')
% for r = 1 : 8
%     for c = 1 : 12
%         fprintf('plotindex = %d,   c=%d, r=%d\n', plotIndex, c, r);
%         % Specify the location for display of the image.
%         subplot(8, 12, plotIndex);
%         % Extract the numerical array out of the cell
%         rgbBlock = ca{r,c};
%         imshow(rgbBlock); % Could call imshow(ca{r,c}) if you wanted to.
%         [rowsB columnsB numberOfColorBandsB] = size(rgbBlock);
%         % Make the caption the block number.
%         caption = [label_ca{r,c} sprintf(': %d', ps_idx_r(r,c))];
%         title(caption);
%         drawnow;
%         % Increment the subplot to the next location.
%         plotIndex = plotIndex + 1;
%     end
% end

% Associating Wells with GEN3 Biolog Experiments %%%%%%%%%%%%%%%%%%%%%%
extract = readcell("/Users/willsharpless/Documents/MATLAB/c112/gen3microplate.csv");
%arranged A1 to 12 then B1 to 12 etc (flipped from our plates)
exlab_wrongorder = extract(2:2:end,1);

i = 1;
ps_exlab = cell(8,12);
for r=1:8
    for c=12:-1:1
       ps_exlab(r,c) = exlab_wrongorder(i);
       i=i+1;
    end
end

exlab_flat = reshape(ps_exlab,[96,1]);

% Display By Cluster (Final Figure) %%%%%%%%%%%%%%%%%%%%%%
ca_flat = reshape(ca_real,[96,1]);

figure('Position', [10 10 2000 2000])
set(gcf,'color','w')


%order them by darkness = distance from darkest point
[~, idxm] = min(vecnorm(Cl'));

[~, idx_v] = sort(vecnorm((Cl-Cl(idxm,:))'));

pi=1;
for a=idx_v
    subplot(numclus,1, pi)
    tbc = ca_flat(idx_l==a);
    clus_labs = exlab_flat(idx_l==a);
    images = [];
    ylabs = [];
    for aa=1:numel(tbc)
        images = [images imresize(tbc{aa},1.5)];
        if aa==numel(tbc)
            ylabs = [ylabs clus_labs{aa}];
        else
            ylabs = [ylabs clus_labs{aa} ', '];
        end
    end
    imshow(images);
%     title(sprintf('Cluster = %d',a));
    ylabel(ylabs, 'FontSize', 20);
    yl = get(gca,'ylabel');
    ylp = get(yl,'Position');
    offset = size(images,2)/2;
    set(get(gca,'YLabel'),'Rotation',0,'Position', ylp+[offset 167 0])
    drawnow;
    pi=1+pi;
end

sgtitle(['\qquad  \textbf{' name '}'], 'FontSize', 30,'Interpreter', 'Latex')

% end

