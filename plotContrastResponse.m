%%  plotContrastResponse.m
%%
%%       usage: interrogator function
%%
%%          by: akshay jagadeesh
%%        date: 03/03/2017
%%     purpose: gets voxel RFs parameters + stimulus image, uses them to 
%%              figure out when the bars cross through each voxel, and plots 
%%              each voxel's peak response amplitude as a function of contrast.
%%
%%       input: 
function plotContrastResponse(v, overlayNum, scanNum, x,y,z,roi)

mglType = 'doubleBars';
%mglType = 'singleBar';

% Load analysis and get computed parameters
analysis = viewGet(v, 'analysis');
if ~isfield(analysis,'type') || ~strcmp(analysis.type,'pRFAnal')
  disp(sprintf('(pRFPlot) pRF analysis has not been run on this scan'));
  return
end
d = viewGet(v, 'd', scanNum);
scanDims = viewGet(v, 'scanDims');
stimfile = viewGet(v, 'stimfile');
rfParams = d.params;

% Define colors for plotting each contrast level
c0625 = [0 .5 0]; c125 = [0 .5 .5];
c675 = [0 0 .5]; c1 = [1 0 .5];
c25 = [0.2 0 0.4];
c75 = [1 .5 0];

% get Stim Image
stimuli = d.stim;
stimY = d.stimY;
stimX = d.stimX;

% Compute stimuli
if strcmp(mglType,'doubleBars')
  conditions = stimfile{1}.stimulus.conditions; 
  seed = 4084807311;
  rng('default'); rng(seed);
  conditionOrder = randsample(size(conditions,1), size(conditions,1), false);
  orderedConds = conditions(conditionOrder(:),:);
  xBar1 = stimfile{1}.stimulus.xBar1;
  yBar1 = stimfile{1}.stimulus.yBar1;
  numSegs = 24;
else
  barCenter = stimfile{1}{1}.stimulus.barCenter;
  barWidth = stimfile{1}{1}.stimulus.barWidth;
  barHeight = stimfile{1}{1}.stimulus.barHeight;
  xBar1 = barCenter(:,1); yBar1 = barCenter(:,2);
  xBar1 = [xBar1-barWidth/2 xBar1-barWidth/2 xBar1+barWidth/2 xBar1+barWidth/2]';
  yBar1 = [yBar1-barHeight/2 yBar1+barHeight/2 yBar1+barHeight/2 yBar1-barHeight/2]';
  angles = [-1; 270; 0; -1; 225; 315;-1; 135; 90; -1; 180; 45];
  orderedConds = [repmat(1,12,1) repmat(1,12,1) angles repmat(-1,12,1)];
  numSegs = size(xBar1,2);
end

% Load time series for this one voxel
disp('Loading voxel time series');
tSeries = loadTSeries(v, scanNum, z, [], x, y); 
tSeries = tSeries(:);
whichVoxel = find(d.linearCoords == sub2ind(scanDims, x, y, z));

% Get RF params for this voxel
thisX = rfParams(1, whichVoxel);
thisY = rfParams(2, whichVoxel);
thisRfWidth = rfParams(3, whichVoxel);
disp(sprintf('Voxel [%d,%d,%d] Receptive Field: (%0.3g, %0.3g) with width %0.3g', x,y,z, thisX, thisY, thisRfWidth));

% Get intersection between each bar and voxel
if strcmp(mglType, 'doubleBars')
  sweepLen = 96;
  numSweeps = length(d.stim)*10;
else
  sweepLen = 48;
  numSweeps = size(orderedConds,1);
end
blankIdx = 0;

%Initialize arrays
intersectLength = 25;
mean0625 = zeros(intersectLength+6,1); count0625 = 0;
mean1 = zeros(intersectLength+6,1); count1 =0;
mean125 = zeros(intersectLength+6,1); count125 = 0;
mean675 = zeros(intersectLength+6,1); count675=0;
mean25 = zeros(intersectLength+6,1); count25=0;
mean75 = zeros(intersectLength+6,1); count75=0;
max0625 = []; max125 = []; max675 = []; max1 = []; max25 = []; max75=[];

% Open a new figure;
figure;

for i = 1:numSweeps
  if ~strcmp(mglType,'doubleBars') % half sweep
    if orderedConds(i,3) == -1
      blankIdx = blankIdx + 1;
      continue
    end
  end
  numNonBlank = i - blankIdx;
  tSweepStart = sweepLen*(numNonBlank-1)+blankIdx*(sweepLen/2); 
  bar1Angle = orderedConds(i,3); bar2Angle = orderedConds(i,4);
  bar1Contrast = orderedConds(i,1); bar2Contrast = orderedConds(i,2);
  cosine = cos(pi*bar1Angle/180); sine = sin(pi*bar1Angle/180);
  barRotMatrix = [cosine sine; -sine cosine];
  initial = 1;

%  if bar2Angle == -1 && bar1Angle ~= -1
    for j = 1:numSegs %iterate through segments
      xc = xBar1(:,j); yc = yBar1(:,j);
      bar1Coords = barRotMatrix * [xc(:,1)'; yc(:,1)'];
      [m1,b1] = calculateLine(bar1Coords(:,1)', bar1Coords(:,2)');
      [m2,b2] = calculateLine(bar1Coords(:,3)', bar1Coords(:,4)');
      distToBar = min(getShortestPath([thisX,thisY], [m1,b1]), getShortestPath([thisX,thisY], [m2,b2]));

      entryTime = tSweepStart+4*(j-1);
%      disp(sprintf('Bar enters receptive field at time %d', entryTime));
      if entryTime +intersectLength <= size(tSeries,1)
        intersectRange = entryTime:entryTime+intersectLength+5;
      else
        continue
      end

      if bar2Angle ~= -1
        barRotMatrix2 = [cos(pi*bar2Angle/180) sin(pi*bar2Angle/180); -sin(pi*bar2Angle/180) cos(pi*bar2Angle/180)];
        bar2Coords = barRotMatrix2 * [xc(:,1)'; yc(:,1)'];
        [m1,b1] = calculateLine(bar2Coords(:,1)', bar2Coords(:,2)');
        [m2,b2] = calculateLine(bar2Coords(:,3)', bar2Coords(:,4)');
        distToBar2 = min(getShortestPath([thisX,thisY], [m1,b1]), getShortestPath([thisX,thisY],[m2,b2]));

        if distToBar2 <= thisRfWidth && distToBar > thisRfWidth*2 && initial == 1
          initial = 0;
          subplot(2,3,5);
          if bar2Contrast == .25
            p25 = plot(tSeries(intersectRange), 'Color', c25); hold on;
            max25(end+1) = max(tSeries(intersectRange));
            mean25 = mean25 + tSeries(intersectRange);
            count25 = count25+1;
          elseif bar2Contrast == .75
            p75 = plot(tSeries(intersectRange), 'Color', c75); hold on;
            max75(end+1) = max(tSeries(intersectRange));
            mean75 = mean75 + tSeries(intersectRange);
            count75 = count75+1;
          end
        elseif distToBar<=thisRfWidth && distToBar2 > thisRfWidth*2 && initial == 1
          initial = 0;
          subplot(2,3,5);
          if bar1Contrast == .25
            p25 = plot(tSeries(intersectRange), 'Color', c25); hold on;
            max25(end+1) = max(tSeries(intersectRange));
            mean25 = mean25 + tSeries(intersectRange);
            count25 = count25+1;
          elseif bar1Contrast == .75
            p75 = plot(tSeries(intersectRange), 'Color', c75); hold on;
            max75(end+1) = max(tSeries(intersectRange));
            mean75 = mean75 + tSeries(intersectRange);
            count75 = count75+1;
          end
        end
      elseif distToBar <= thisRfWidth && initial == 1

        initial = 0;

        if bar1Contrast == .0625
          subplot(2,3,1); p1 = plot(tSeries(intersectRange), 'Color', c0625); hold on;
          max0625(end+1) = max(tSeries(intersectRange));
          mean0625 = mean0625 + tSeries(intersectRange);
          count0625 = count0625 + 1;
        elseif bar1Contrast == .125
          subplot(2,3,1); p2 = plot(tSeries(intersectRange), 'Color', c125); hold on;
          max125(end+1) = max(tSeries(intersectRange));
          mean125 = mean125 + tSeries(intersectRange);
          count125 = count125 + 1;
        elseif bar1Contrast == .675
          subplot(2,3,1); p3 = plot(tSeries(intersectRange), 'Color', c675); hold on;
          max675 = max(tSeries(intersectRange));
          mean675 = mean675 + tSeries(intersectRange);
          count675 = count675 + 1;
        elseif bar1Contrast == 1
          subplot(2,3,1); p4 = plot(tSeries(intersectRange), 'Color', c1); hold on;
          max1 = max(tSeries(intersectRange));
          mean1 = mean1 + tSeries(intersectRange);
          count1 = count1 + 1;
        end 
        continue
      end
    end % End for loop through sweeps
%  end
end
%keyboard
subplot(2,3,1); title('Voxel time courses to all single bar sweeps');
%legend([p1,p2,p3,p4],'6.25%', '12.5%', '67.5%', '100%');

% Average by number of times this voxel saw each contrast
mean0625 = mean0625./count0625;
mean1 = mean1./count1;
mean125 = mean125./count125;
mean675 = mean675./count675;
mean25 = mean25./count25;
mean75 = mean75./count75;

% Plot mean response time courses to each contrast
subplot(2,3,2);
plot(mean0625, 'Color', c0625, 'LineWidth', 0.75); hold on; plot(mean125, 'Color', c125); hold on;
plot(mean675, 'Color', c675, 'LineWidth', 0.75); hold on; plot(mean1, 'Color', c1);
plot(mean25, 'Color', c25, 'LineWidth', 0.75); hold on; plot(mean75, 'Color', c75);
legend('6.25%', '12.5%', '67.5%', '100%', '25%', '75%');
title('Mean response time course to each contrast');

% Plot peak amplitude of each mean time course
subplot(2,3,3);
meanMax = [mean(max0625), mean(max125),mean(max25), mean(max675), mean(max75),mean(max1)];
bGraph = bar(1:6, diag(meanMax)); 
set(bGraph(1), 'facecolor', c0625); set(bGraph(2), 'facecolor', c125); set(bGraph(3), 'facecolor', c25); set(bGraph(4), 'facecolor', c675); set(bGraph(5), 'facecolor', c75); set(bGraph(6), 'facecolor', c1);
legend('6.25%', '12.5%', '25%','67.5%','75%', '100%');
title('Mean peak response amplitude by contrast');
xlabel('Bar Contrast'); ylabel('Mean peak response across voxels'); ylim([min(meanMax)-.05 max(meanMax)+.05]);

% Plot voxel RF
subplot(2,3,4);
cPlot = circle(thisX, thisY, thisRfWidth); hold on;
title('Voxel Receptive Field Location');
xlim([stimX(1), stimX(end)]); ylim([stimY(1), stimY(end)]);
hline(0,':'); vline(0,':');


keyboard

%%%%%% < end of main > %%%%%%%%%%

%----------------------------------------------
% calculateLine
%       calculates the line between two points
%   
%   Arguments:
%       p1: [x1 y1]
%       p2: [x2 y2]
%----------------------------------------------
function [slope, intercept] = calculateLine(p1, p2)
slope = (p2(2)-p1(2)) / (p2(1) - p1(1));
if slope == Inf % return the x-value instead
  intercept = p2(1);
else
  intercept = p2(2) - slope*p2(1);
end

%----------------------------------------------
% getShortestPath
%     gets the distance of the shortest path from a point to a line
%
%  Arguments:
%      point: [x, y]
%      line: [slope, intercept]
%----------------------------------------------
function distance = getShortestPath(point, line)
x1 = point(1); y1 = point(2);
m = line(1); b = line(2);
if m == Inf
  distance = abs(x1 - b);
else
  distance = abs(-m*x1 + y1 + (-b))/(sqrt(1+(-m)^2));
end


%----------------------------------------------
% circle
%    plots a circle, given center and radius
%
%----------------------------------------------
function cPlot = circle(x,y, r)
angle = 0:0.01:2*pi;
xp = r*cos(angle);
yp = r*sin(angle);
cPlot = plot(x+xp, y+yp, 'g');

