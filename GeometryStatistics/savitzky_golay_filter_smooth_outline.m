% Demo to smooth a 2D outline (boundary) of a blob using Savitzky-Golay filters.
function savitky_golay_filter_smooth_outline()
clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 25;

% Read in a standard MATLAB gray scale demo image.
button = menu('Use which demo image?', 'CameraMan', 'Moon', 'MRI', 'Glass', 'Pout');
thresholdValue = 125;
if button == 1
	baseFileName = 'cameraman.tif';
	thresholdValue = 75;
elseif button == 2
	baseFileName = 'moon.tif';
	thresholdValue = 25;
elseif button == 3
	baseFileName = 'mri.tif';
	thresholdValue = 33;
elseif button == 4
	baseFileName = 'glass.png';
	thresholdValue = 125;
else
	baseFileName = 'pout.tif';
	thresholdValue = 125;
end

%===============================================================================
% Read in a standard MATLAB gray scale demo image.
folder = fileparts(which('cameraman.tif')); % Determine where demo folder is (works with all versions).
% baseFileName = 'moon.tif';
% Get the full filename, with path prepended.
fullFileName = fullfile(folder, baseFileName);
% Check if file exists.
if ~exist(fullFileName, 'file')
	% File doesn't exist -- didn't find it there.  Check the search path for it.
	fullFileNameOnSearchPath = baseFileName; % No path this time.
	if ~exist(fullFileNameOnSearchPath, 'file')
		% Still didn't find it.  Alert user.
		errorMessage = sprintf('Error: %s does not exist in the search path folders.', fullFileName);
		uiwait(warndlg(errorMessage));
		return;
	end
end
grayImage = imread(fullFileName);
% Get the dimensions of the image.  
% numberOfColorBands should be = 1.
[rows, columns, numberOfColorBands] = size(grayImage);
if numberOfColorBands > 1
	% It's not really gray scale like we expected - it's color.
	% Convert it to gray scale by taking only the green channel.
	grayImage = grayImage(:, :, 2); % Take green channel.
end
% Display the original gray scale image.
subplot(2, 3, 1);
imshow(grayImage, []);
axis on;
title('Original Grayscale Image', 'FontSize', fontSize, 'Interpreter', 'None');
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Give a name to the title bar.
set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off') 

% Let's compute and display the histogram.
[pixelCount, grayLevels] = imhist(grayImage);
subplot(2, 3, 2); 
bar(grayLevels, pixelCount); % Plot it as a bar chart.
grid on;
title('Histogram of original image', 'FontSize', fontSize, 'Interpreter', 'None');
xlabel('Gray Level', 'FontSize', fontSize);
ylabel('Pixel Count', 'FontSize', fontSize);
xlim([0 grayLevels(end)]); % Scale x axis manually.

% Threshold the image
binaryImage = grayImage > thresholdValue;
% Get rid of holes in the blobs.
binaryImage = imfill(binaryImage, 'holes');
%---------------------------------------------------------------------------
% Extract the largest area using our custom function ExtractNLargestBlobs().
numberToExtract = 1;
biggestBlob = ExtractNLargestBlobs(binaryImage, numberToExtract);
%---------------------------------------------------------------------------
% Display the binary image.
subplot(2, 3, 3);
imshow(biggestBlob, []);
title('Binary Image', 'FontSize', fontSize, 'Interpreter', 'None');

% Get the boundaries.
% First, display the original gray scale image.
subplot(2, 3, 4);
imshow(grayImage, []);
axis image; % Make sure image is not artificially stretched because of screen's aspect ratio.
hold on;
% Now get the boundaries.
boundaries = bwboundaries(biggestBlob);
numberOfBoundaries = size(boundaries, 1);
for k = 1 : numberOfBoundaries
	thisBoundary = boundaries{k};
	plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 2);
end
hold off;
caption = sprintf('%d Outlines, from bwboundaries()', numberOfBoundaries);
title(caption, 'FontSize', fontSize); 
axis on;

firstBoundary = boundaries{1};
% Get the x and y coordinates.
x = firstBoundary(:, 2);
y = firstBoundary(:, 1);

% Now smooth with a Savitzky-Golay sliding polynomial filter
windowWidth = 45
polynomialOrder = 2
smoothX = sgolayfilt(x, polynomialOrder, windowWidth);
smoothY = sgolayfilt(y, polynomialOrder, windowWidth);
% First, display the original gray scale image.
subplot(2, 3, 5);
imshow(grayImage, []);
axis image; % Make sure image is not artificially stretched because of screen's aspect ratio.
hold on;
plot(smoothX, smoothY, 'r-', 'LineWidth', 2);
grid on;
title('Smoothed Boundary', 'FontSize', fontSize);

%==============================================================================================
% Function to return the specified number of largest or smallest blobs in a binary image.
% If numberToExtract > 0 it returns the numberToExtract largest blobs.
% If numberToExtract < 0 it returns the numberToExtract smallest blobs.
% Example: return a binary image with only the largest blob:
%   binaryImage = ExtractNLargestBlobs(binaryImage, 1);
% Example: return a binary image with the 3 smallest blobs:
%   binaryImage = ExtractNLargestBlobs(binaryImage, -3);
function binaryImage = ExtractNLargestBlobs(binaryImage, numberToExtract)
try
	% Get all the blob properties.  Can only pass in originalImage in version R2008a and later.
	[labeledImage, numberOfBlobs] = bwlabel(binaryImage);
	blobMeasurements = regionprops(labeledImage, 'area');
	% Get all the areas
	allAreas = [blobMeasurements.Area];
	if numberToExtract > length(allAreas);
		% Limit the number they can get to the number that are there/available.
		numberToExtract = length(allAreas);
	end
	if numberToExtract > 0
		% For positive numbers, sort in order of largest to smallest.
		% Sort them.
		[sortedAreas, sortIndexes] = sort(allAreas, 'descend');
	elseif numberToExtract < 0
		% For negative numbers, sort in order of smallest to largest.
		% Sort them.
		[sortedAreas, sortIndexes] = sort(allAreas, 'ascend');
		% Need to negate numberToExtract so we can use it in sortIndexes later.
		numberToExtract = -numberToExtract;
	else
		% numberToExtract = 0.  Shouldn't happen.  Return no blobs.
		binaryImage = false(size(binaryImage));
		return;
	end
	% Extract the "numberToExtract" largest blob(a)s using ismember().
	biggestBlob = ismember(labeledImage, sortIndexes(1:numberToExtract));
	% Convert from integer labeled image into binary (logical) image.
	binaryImage = biggestBlob > 0;
catch ME
	errorMessage = sprintf('Error in function ExtractNLargestBlobs().\n\nError Message:\n%s', ME.message);
	fprintf(1, '%s\n', errorMessage);
	uiwait(warndlg(errorMessage));
end
