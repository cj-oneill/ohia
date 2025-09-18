function scene = create_independent_scene()
%CREATE_INDEPENDENT_SCENE Creates a fresh, identical radar scenario for a parallel worker.

% This function contains all the necessary setup code copied from the main script.

%% Copied from "Generate Simulated Terrain"
rng(2021);
xLimits         = [900 1200];
yLimits         = [-200 200];
roughnessFactor = 1.75;
initialHgt      = 0;
initialPerturb  = 200;
numIter         = 8;
[x,y,A] = helperRandomTerrainGenerator(roughnessFactor,initialHgt, ....
    initialPerturb,xLimits(1),xLimits(2), ...
    yLimits(1),yLimits(2),numIter);
A(A < 0) = 0;
xvec = x(1,:); 
yvec = y(:,1);

%% Copied from "Specify the SAR System and Scenario"
freq = 1e9;
[lambda,~] = freq2wavelen(freq);
bw = 30e6;
fs = 60e6;
tpd = 3e-6;
rngRes = bw2rangeres(bw);
apertureLength = 6;
v = 100;
dur = 1;
rdrhgt = 1000;
rdrpos1 = [0 0 rdrhgt];
rdrvel = [0 v 0];
len = sarlen(v,dur);
targetpos = [1000,len/2,0;1020,len/2,0;1040,len/2,0];
tgthgts = 110*ones(1,3);
for it = 1:3
    [~,idxX] = min(abs(targetpos(it,1) - xvec)); 
    [~,idxY] = min(abs(targetpos(it,2) - yvec)); 
    tgthgts(it) = tgthgts(it) + A(idxX,idxY); 
    targetpos(it,3) = tgthgts(it); 
end
rc = sqrt((rdrhgt - mean(tgthgts))^2 + (mean(targetpos(:,1)))^2);
depang = depressionang(rdrhgt,rc,'Flat','TargetHeight',mean(tgthgts));
prf = 500;
scene = radarScenario('UpdateRate',prf,'IsEarthCentered',false,'StopTime',dur);
rdrplat = platform(scene,'Trajectory',kinematicTrajectory('Position',rdrpos1,'Velocity',[0 v 0]));
rcs = rcsSignature('Pattern',5); 
for it = 1:3
    platform(scene,'Position',targetpos(it,:),'Signatures',{rcs});
end

%% Copied from "Define the Land Surface Reflectivity"
grazTable = 20:0.1:60;
freqTable = [1e9 10e9]; 
numSurfaces = 2;
reflectivityLayers = zeros(numel(grazTable),numel(freqTable),numSurfaces);
reflectivityLayers(:,:,1) = landreflectivity('Woods', grazTable,freqTable);
reflectivityLayers(:,:,2) = landreflectivity('WoodedHills', grazTable,freqTable);
reflectivityType = ones(size(A)); 
reflectivityType(A > 100) = 2; 
reflectivityMap = surfaceReflectivity('Custom','Frequency',freqTable, ...
    'GrazingAngle',grazTable,'Reflectivity',reflectivityLayers, ...
    'Speckle','Rayleigh');
landSurface(scene,'Terrain',A,'Boundary',[xLimits;yLimits], ...
    'RadarReflectivity',reflectivityMap, ...
    'ReflectivityMap',reflectivityType);

%% Copied from "Configure the Radar Transceiver"
maxRange = 2500;
mountAngles = [0 depang 0];
rdr = radarTransceiver('MountingAngles',mountAngles,'NumRepetitions',1, ...
    'RangeLimits',[0 maxRange]);
rdr.Transmitter.PeakPower = 50e3; 
rdr.Receiver.SampleRate = fs;
rdr.Receiver.NoiseFigure = 30; 
antbw = ap2beamwidth(apertureLength,lambda); 
ant = phased.SincAntennaElement('FrequencyRange',[1e9 10e9],'Beamwidth',antbw);
rdr.TransmitAntenna.Sensor = ant;
rdr.TransmitAntenna.OperatingFrequency = freq;
rdr.ReceiveAntenna.Sensor = ant;
rdr.ReceiveAntenna.OperatingFrequency = freq;
antennaGain = aperture2gain(apertureLength^2,lambda); 
rdr.Transmitter.Gain = antennaGain;
rdr.Receiver.Gain = antennaGain;
rdr.Waveform = phased.LinearFMWaveform('SampleRate',fs,'PulseWidth',tpd, ...
    'PRF',prf,'SweepBandwidth',bw); 
rdrplat.Sensors = rdr;

%% Copied from "Generate the Datacube"
clutterGenerator(scene,rdr,'Resolution',rngRes/2, 'RangeLimit',maxRange);

end