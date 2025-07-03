% 1. input: table from generates_barcodes and psf from SMAP

%% Generate/import the codebook
% generate code book
numGenes = 100;
numBits = 16;
numOnes = 4;
% codebook = generate_barcodes(numGenes, numBits, numOnes, 1);

% import the codebook
codebook = read_codebook('codebook\C1E1_codebook.csv');

%% Load the Point Spread Function from SMAP
load('psf\1006PSF1000_3dcal');
psf = SXY.PSF{1};
psf = psf/sum(psf(:));
[numRows, numCols] = size(psf(:, :, 1));


% Grid dimension (1608 x 1608 x 401 layers)
xDim = 300 + size(psf, 1) - 1;
yDim = 300 + size(psf, 2) - 1;
zDim = size(psf, 3);

% Create a circular binary mask
radius = 14./2;
maskX = -(floor(numRows./2)):(floor(numRows./2));
maskY = -(floor(numCols./2)):(floor(numCols./2));
[X, Y] = meshgrid(maskX, maskY);
circlemask = (X).^2+(Y).^2 >= radius.^2;

for i = 1:zDim
    zslice = psf(:,:,i);
    meanedge = sum(zslice(circlemask))/sum(circlemask(:));
    zslice = zslice - meanedge;
    zslice(zslice<0) = 0;
    psf(:,:,i) = zslice;
end

%% Define Simulation Parameters

% Mean number of emitters (poisson)
N = 200;

% Camera Parameters
QE = [0.87, 0.95, 0.89, 0.71]; % Quantum efficiency - 473, 561, 647, 750
gain = 1.33; % e/ADU
bias = 100; % ADU
dark = 0; % e/s
readNoise = 1.6; % electrons standard dev


%% Simulate flourescent probes

% Generate the number of emitters
volume = xDim*yDim*zDim;
lambda=N/volume; % density (volume)
npoints = poissrnd(lambda*volume);

% Randomly distribute the points within the volume
pproc = rand(npoints,3);
% emitterPos = pproc.*[xDim, yDim, zDim];
emitterPos = zeros(size(pproc));
emitterPos(:,1) = pproc(:,1)*xDim;
emitterPos(:,2) = pproc(:,2)*yDim;
emitterPos(:,3) = pproc(:,3)*zDim;

% r contain the number of lines in the table for the genes
r = randi([1 numGenes], 1, npoints);

% Creating the ground truth table
varTypes = ["int32","int32","int32","string","string"];
varNames = ["x","y","z","genes","barcodes"];
groundTruthTable = table('Size',[npoints,5],'VariableTypes',varTypes,'VariableNames',varNames);
emitterBarcodes = false(npoints, numBits); % logical matrix that contains emitter barcodes


% Create a folder to store the images
t = datetime('now','Format','dd-MMM-yyyy HH-mm-ss');
S = char(t);
dirPath = fullfile('generated_data', [S, ' merfish_sim']);
mkdir(dirPath);

for i = 1:npoints
    % grab gene and barcode from the codebook
    gene = codebook{r(i), 1};
    barcode = codebook{r(i), 2};
    
    groundTruthTable(i,:) = {emitterPos(i,2) - floor(size(psf, 2)/2), emitterPos(i,1) - floor(size(psf, 1)/2), emitterPos(i,3), gene, barcode};
    emitterBarcodes(i,:) = logical(str2num(barcode));  %#ok<ST2NM>
end


%Main loop
for i = 1:numBits

    % Get emitter positions that will light up this round
    onEmitterPos = emitterPos(emitterBarcodes(:,i), :);
    
    % Create an empty matrix to represent emitter positions in 3D
%     meanBackground = 20;
%     stddevBackground = 5;
    meanEmitter = 50000000;
    stddevEmitter = 200;
    backgroundMatrix = zeros(xDim, yDim, zDim);
    
    for j = 1:size(onEmitterPos,1)
        pos = ceil(onEmitterPos(j,:));
        backgroundMatrix(pos(1), pos(2), pos(3)) = normrnd(meanEmitter, stddevEmitter);
    end
    
    %conv with psf
    photonIm = convn(backgroundMatrix, psf, 'valid');

    % Add photon shot noise
    photonIm = arrayfun(@(x) poissrnd(x), photonIm);
    
    % Quantum efficiency (Electron conversion - photodetector)
    electronIm = photonIm*QE(2);
    
    % Adding in Camera Noise
    [row, col] = size(electronIm);
    readNoiseIm = normrnd(0, readNoise, row, col);
    electronIm2 = electronIm + readNoiseIm;
    
    % Converting electrons to ADUs
    digitalIm = electronIm2 ./gain + bias;
    fprintf('Round %i completed \n', i);

    % Save the image
    imwrite(uint16(digitalIm), fullfile(cd,dirPath,"bit" + num2str(i) + ".Tiff"));
end

%% Save relevant files

fprintf("Done! \n");
writetable(groundTruthTable, fullfile(cd,dirPath,'ground_truth.txt'));
        