tic;

import matlab.io.*
dataFolder = 'D:\s1d\19mar2016\';
files = dir(fullfile(dataFolder,'*s1d_A.fits'));

%files = files(1:20:end);
N = length(files);
s1dDate = zeros(1,N);

for Run = 1:N
    
    dataFileA = files(Run).name;
    
    %Extract date/time info
    s1dDate(Run) = datenum([dataFileA(7:16),', ',dataFileA(18:29)],'yyyy-mm-dd, HH-MM-SS');
end

% Days = unique(floor(s1dDate));
% Days = Days + 13/24; %Find date closest to 1 PM each day
% M = length(Days);
% s1dDateSel = zeros(1,M);

M = N;

Efficiencies = zeros(1,M);
AirMass = zeros(1,M);
Centroids = zeros(1,M);
MeanCounts =zeros(1,M);
BERV = zeros(1,M);
BERVMX = zeros(1,M);

for Run2 = 1:M
        Run = Run2;
%     [~,Run] = min(abs(CloudCutDates - Days(Run2)));
%     [~,Run] = min(abs(CloudCutDates(Run) - s1dDate));
%     s1dDateSel(Run2) = s1dDate(Run);
    dataFileA = files(Run).name;
    
    fptr = fits.openFile([dataFolder, dataFileA]);
    AirMass(Run2) = str2double(fits.readKey(fptr,'AIRMASS'));
    Centroids(Run2) =str2double(fits.readKey(fptr,'TNG EXP_METER_A EXP CENTROID'));
    Efficiencies(Run2) = str2double(fits.readKey(fptr,'TNG EXP_METER_A EFFICENCY MEAN'));
    MeanCounts(Run2) = str2double(fits.readKey(fptr,'TNG EXP_METER_A COUNTS MEAN'));
    BERV(Run2) = str2double(fits.readKey(fptr,'TNG DRS BERV'));
    BERVMX(Run2) = str2double(fits.readKey(fptr,'TNG DRS BERVMX'));
    %Fix HARPS-N Barycentric correction. Remove HARPS-N DRS correction, then
    %subtract off JPL calculated Sun-La Plama relative velocity
    
    v_JPL=1e3*interp1(DaTJPL,RV_JPL,s1dDate(Run),'linear','extrap');
    v_DRS=-1e3*BERV(Run2);
    %v_DRS = 1e3*interp1(DaTHARPSN,vrad,s1dDate(Run),'linear','extrap');
    
    %Load spectra/comb
    %s1d files have different numbers of pixels. UGH. For now, let's use
    %the first spectrum as a template, and interpolate onto that one's
    %wavelength solution.
    if Run2 == 1
        j=0;
        firstspectra = fitsread([dataFolder, dataFileA]);
        spectra = nan(M,length(firstspectra));
        wavelengthsSave = nan(M,length(firstspectra));
        spectra(Run2,:) = firstspectra;
        spectrainfo = fitsinfo([dataFolder, dataFileA]);
        spectraRef = str2double(fits.readKey(fptr,'CRVAL1'));
        spectraPixRef = str2double(fits.readKey(fptr,'CRPIX1'));
        spectraInc = str2double(fits.readKey(fptr,'CDELT1'));
        spectraStart = spectraRef - spectraInc*(spectraPixRef-1);
        spectraWavelength = (spectraStart:spectraInc:spectraStart + spectraInc*(length(firstspectra)-1))/10;
        %spectraWavelength = spectraWavelength*(1-(v_JPL-v_DRS)/3e8);
        wavelengthsSave(Run2,:) = spectraWavelength;
        clear firstspectra
    else
        %Determine wavelength ranges
        spectrainfo = fitsinfo([dataFolder, dataFileA]);
        spectraRef = str2double(fits.readKey(fptr,'CRVAL1'));
        spectraPixRef = str2double(fits.readKey(fptr,'CRPIX1'));
        spectraInc = str2double(fits.readKey(fptr,'CDELT1'));
        spectraStart = spectraRef - spectraInc*(spectraPixRef-1);
        templatespectra = fitsread([dataFolder, dataFileA]);
        templateWavelength = (spectraStart:spectraInc:spectraStart + spectraInc*(length(templatespectra)-1))/10;
        %templateWavelength = templateWavelength*(1-(v_JPL-v_DRS)/3e8);
        if length(templateWavelength) > size(spectra,2)
            spectra = [spectra,nan(M,length(templateWavelength)-size(spectra,2))];
            wavelengthsSave = [wavelengthsSave,nan(M,length(templateWavelength)-size(wavelengthsSave,2))];
            spectra(Run2,:) = templatespectra;
            wavelengthsSave(Run2,:) = templateWavelength;
            j = j+1;
        else
            spectra(Run2,1:length(templatespectra)) = templatespectra;
            wavelengthsSave(Run2,1:length(templateWavelength)) = templateWavelength;
        end
%        spectra(Run2,:) = interp1(templateWavelength,templatespectra,spectraWavelength,'linear','extrap');
    end
    fits.closeFile(fptr);

    
end

wl_min = max(min(wavelengthsSave'));
wl_max = min(max(wavelengthsSave'));
%wl_max = 683.2486; %Start of order 69
wavelengthgrid = wl_min:spectraInc/10:wl_max;
spectragrid = zeros(M,length(wavelengthgrid));
for Run2 = 1:M
    templatespectra = spectra(Run2,:);
    templatespectra = templatespectra(isnan(templatespectra)==0);
    templateWavelength = wavelengthsSave(Run2,:);
    templateWavelength = templateWavelength(isnan(templateWavelength)==0);
    spectragrid(Run2,:) = interp1(templateWavelength,templatespectra,wavelengthgrid,'linear','extrap');
end
toc
%%
clearvars -except wavelengthgrid spectragrid s1dDateSel

%% Divide out overall changes in shape
spectragrid2 = spectragrid;
for Run3 = 1:M
   
    spectragrid2(Run3,:) = spectragrid(Run3,:)./polyval(pmat(Run3,:),fitgrid);
    
end

%% Remove waterlines from HARPS-N mask
spectragrid3 = spectragrid2;

for Run4 = 1:length(Bandpass)
    blockinds = find(abs(wavelengthgrid - Bandpass(Run4))<0.015);
    spectragrid3(:,blockinds) = 0;
end

%% Do PCA
[coeff,~,var] = pca(spectragrid);
figure;scatter(wavelengthgrid,spectragrid(1,:),4,coeff(:,1));shg
colorbar; colormap('jet')

%% Find lines of interest
sigma = sqrt(mean(coeff(:,3).^2));

locs = find(abs(coeff(:,3)) > 5*sigma);
minlocs = zeros(size(locs));
r = 10;
for kk = 1:length(locs)
    
    spectratest = spectragrid3(1,locs(kk)-r:locs(kk)+r);
    [~, minInd] = min(spectratest);
    minlocs(kk) = (locs(kk)-r+minInd)-1;
    
    %         spectratest = spectragrid(1,locs(kk)mi-r:end);
    %         [~, minInd] = min(spectratest);
    %         minlocs(kk) = (locs(kk)-r+minInd);
    
end

minlocs = unique(minlocs);

%%
sigma = sqrt(mean(coeff(:,3).^2));

locs = find(abs(coeff(:,3)) > 5*sigma);
minlocs = zeros(size(locs));

r = 40;
specdiff = diff(spectragrid(1,:));

for kk = 1:length(locs)
    if specdiff(locs(kk)) > 0
        [~,minlocs(kk)] = findpeaks(-spectragrid(1,locs(kk)-r:locs(kk)+5),'NPeaks',1);
        minlocs(kk) = minlocs(kk) + locs(kk)-r;
    else
        [~,minlocs(kk)] = findpeaks(-spectragrid(1,locs(kk)-5:locs(kk)+r),'NPeaks',1);
        minlocs(kk) = minlocs(kk) + locs(kk)-5;
    end
end
minlocs = unique(minlocs);

%%
derivmat = zeros(length(locs),length(s1dDateSel));

for day = 1:length(s1dDateSel)
    for pix = 1:length(locs)
        
        derivmat(pix,day) = (spectragrid(day,locs(pix)+1) - spectragrid(day,locs(pix)-1)) / (wavelengthgrid(locs(pix) + 1) - wavelengthgrid(locs(pix)-1));
        
    end
end