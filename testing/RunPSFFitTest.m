%RunPSFFitTest


%% Set parameters
optionsSim.ampMean = 300; optionsSim.ampSigma = 100; % peak amplitude above background in photons
optionsSim.background = 0; % background count in photons
optionsSim.PSFsigma = 1.25; % PSF stdev in pixels
optionsSim.halfwMask = 5; % half window size of PSF mask
optionsSim.sigmaRange = [optionsSim.PSFsigma-0.6,optionsSim.PSFsigma+0.6];
optionsSim.thetaRange = 20; % rotation of Gaussian function. Set to 0 to simulate perfect axis alignment
optionsSim.imWidth = 64;

% fit parameters
optionsFit.usePixelIntegratedFit = false; %use pixel integrated fit, not possible when fitting angles
optionsFit.useMLE = false; % use MLE
optionsFit.halfw = 5; %half window size of fit mask
optionsFit.varsToFit = [1,1,1,1,1,1,1]; %fit parameters [x,y,A,BG,sigma_x,sigma_y,angle] true/false

plotFits = true;

% Create candidate positions
nr_candidates = 10;
candidatePos = randi([optionsSim.halfwMask,optionsSim.imWidth-optionsSim.halfwMask],nr_candidates,2);

%% Create simulated image with predetermined positions
[candidatePos,img,img_truth,ground_truth] = createPSFFitTest(candidatePos, optionsSim,optionsFit,plotFits);


%% Test performance
[params] = psfFit_Image( img, candidatePos.',optionsFit.varsToFit,optionsFit.usePixelIntegratedFit,optionsFit.useMLE,optionsFit.halfw,optionsSim.PSFsigma);
fitData = {params(:,params(end,:)==1).'};

figure; mim(img); hold on; plot(fitData{1}(:,1),fitData{1}(:,2),'bo'); hold off;

tic;
for i=1:1e2
[params] = psfFit_Image( img, candidatePos.',optionsFit.varsToFit,optionsFit.usePixelIntegratedFit,optionsFit.useMLE,optionsFit.halfw);
end
toc
