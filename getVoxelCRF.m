%%  getVoxelCRF.m
%%
%%       usage: 
%%
%%          by: akshay jagadeesh
%%        date: 03/03/2017
%%     purpose: gets voxel RFs parameters + stimulus image, uses them to 
%%              figure out when the bars cross through each voxel, and plots 
%%              each voxel's peak response amplitude as a function of contrast.
%%
%%       input: 
%


function getVoxelCRF(pRFAnal, stimImage)

%%%%% Hardcoded parameters
%analysis = 'pRF_gaussprefit10.mat';
analysis = 'pRF_gaussprefit.mat';
scanNum = 3; %% 3=concat2 -- change to 2 for concat10

% Create new view and load analysis
v = newView;
v = viewSet(v, 'curGroup', 'Concatenation');
v = loadAnalysis(v, ['pRFAnal/' analysis]);
d = viewGet(v, 'd', scanNum);

% Get model RF parameters
rfParams = d.params;

% get Stim Image
stimuli = d.stim;
im1 = stimuli{1}.im;

% best voxel
[a,b] = max(d.r(:,1));
bestParams = params(:,b); x = bestParams(1); y = bestParams(2); rfWidth = bestParams(3);

stimY = d.stimY;
stimX = d.stimX;

% Step 1: Convert stimulus image to screen coordinates

% Step 2: Figure out when stimulus image at the rf is nonzero

keyboard
% Iterate through each voxel
%for i = 1:size(rfParams,2)

  % Iterate through time points and find the stim time when the stimulus intersects with this voxel
%  for i = 1:size(im1, 3)     
    
    

%  end
  % Separate by different contrast bars and average across same contrasts
  % Grab the time series response when it is exposed to bar of each contrast


%end

%% Start by doing this for best voxel


keyboard
