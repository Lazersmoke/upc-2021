h = 6.63e-34; % J s
c = 2.99e8; % m/s

%% Load in absorbance coefficients
XSdat = textscan(fopen('CO2_Cross_Sectional_Absorbance_at_current_Mars.txt'),'%f%f%f','headerlines',6);
lambdasC = 100./XSdat{1};
XSC = XSdat{2};
newLs = [];
newXS = [];
for i = 1:size(lambdasC)
    if XSC(i) > 0
        newLs(end+1) = lambdasC(i);
        newXS(end+1) = XSC(i);
    end
end
lambdasC = newLs; % m
XSC = newXS; % cm^2 / molecule
XSC = XSC * (6.022e23/44.01*1000) * 10^-4; % m^2/kg


XSdat = textscan(fopen('Watervaporcrosssection.txt'),'%f%f%f','headerlines',6);
lambdasW = 100./XSdat{1};
XSW = XSdat{2};
newLs = [];
newXS = [];
for i = 1:size(lambdasW)
    if XSW(i) > 0
        newLs(end+1) = lambdasW(i);
        newXS(end+1) = XSW(i);
    end
end
lambdasW = newLs; % m
XSW = newXS; % cm^2 / molecule
massAvgWater = 0.01801528;
XSW = XSW * (6.022e23/massAvgWater) * 10^-4; % m^2/kg

%% Preparing the initial conditions from the viking data
numLayers = 7*5;
zPerLayer = 3; %km
zmax = numLayers*zPerLayer;

vdat = textscan(fopen('viking.txt'),'%d%f%f%f%f%f%f','headerlines',1);

zs = flip(vdat{3});
rho = flip(vdat{4}); 
pres = flip(vdat{5});
temp = flip(vdat{6});


pressuresC = zeros(numLayers,1); % kg / m^3
rhosC = zeros(numLayers,1); % mb
temps = zeros(numLayers,1); % K
ind = 1;
for i = 1:numLayers
    maxZ = zPerLayer * i;
    tempAve = 0;
    rhoAve = 0;
    pAve = 0;
    counter = 0;
    while zs(ind) < maxZ
        tempAve = tempAve + temp(ind);
        rhoAve = rhoAve + rho(ind);
        pAve = pAve + pres(ind);
        ind = ind + 1;
        counter = counter + 1;
    end
    if counter > 0 
        tempAve = tempAve / counter;
        rhoAve = rhoAve / counter;
        pAve = pAve / counter;
    end
    temps(i) = tempAve;
    rhosC(i) = rhoAve;
    pressuresC(i) = pAve;
end

zs = linspace(0,zmax .* 1000,numLayers); % km

diffTemp = 333.15 - temps(1);
for i = 1:numLayers
    temps(i) = temps(i) + diffTemp*10^(-zs(i)/10000); %adjusts base temp to be
    %60, assumes added temp falls by factor of 1/10 every 10 km
end


%% Time stepping


% Day night cylce (heat and cool) and over the net it cools

atmosTemps = temps;
surfaceTemp = atmosTemps(1);
surfaceTemp = 273 + 60;
atmosTemps = atmosTemps;

massAvgWater = 0.01801528;
p0Initial = 0.2 * 101325; % Pa

presZero = computePressures(150,0.2,zs,massAvgWater);
densZero = computeDensities(150,0.2,zs,massAvgWater);
intens = computeIntensities(densZero,XSW,lambdasW,zs,surfaceTemp); % dE/dt function of z*

dz = diff(zs);
dz(end+1) = dz(end);

cpwater = 4186; % J/(kg C)

omegaDayNight = 2 * pi * 0.1;
    
dt = 0.1;
close all
figure(1)
for t = 1:dt:1000
    phiDayNight = omegaDayNight * t;
    dayNightCycle = max(0,cos(phiDayNight));
    
%     surfaceTemp = atmosTemps(1); % Equilibrium with the first layer
    dens = computeDensities(mean(atmosTemps),p0Initial,zs,massAvgWater);
    
    % Radiate surface to atmosphere
    intens = computeIntensities(dens,XSW,lambdasW,zs,surfaceTemp);
    
    delta = 1 ./ dens' ./ cpwater .* intens' ./ dz';
    
    
    % Sun heats surface
    
    h = 6.63e-34; % J s
    c = 2.99e8; % m/s
    kb = 1.38e-23;
    sigmaSF = 2 * pi^5 * kb^4/(15 * c^2 * h^3);
    rSun = 7e8;
    rMars = 3.3895e6;
    dSunMars = 2.395e11;
    sunTemp = 5778;
    jStar = sigmaSF * sunTemp^4;
    Lsun = jStar * 4 * pi * rSun^2;
    
    
%     cpMartianSurface = cpwater;
%     cpMartianSurface = 0.16e26 / (4 * pi * rMars^2);
    cpMartianSurface = 741;
%     741 J/(kg K)
    
    % Amount of sun per meter hitting mars
    hittingMars = dayNightCycle .* Lsun./(4 * pi * dSunMars^2);
    jStarMars = sigmaSF .* surfaceTemp^4;
    % Neglect conduction
    

    backRad = zeros(size(intens));
    heatLossAtm = zeros(size(intens));
    % Gas radiates to space
    for zi = 1:size(zs,2)
%         idxDown = flip(1:zi);
%         computeIntensities(dens(idxDown),)
        backRadAmt = getBBSpec(atmosTemps(zi),lambdasW) .* 0.4 .* dens(zi) .* dz(zi);
        backRad(zi) = sum(backRadAmt .* 0.1); % dLambda = 0.1
        heatLossAtm(zi) = sigmaSF .* atmosTemps(zi)^4;
    end
    hittingMars;
    jStarMars;
    backRadiation = sum(backRad) ./ cpMartianSurface;
%     delta = delta - backRad' ./ cpwater; % Atmosphere loses heat to surface
    delta = delta - heatLossAtm' ./ cpwater; % Atmosphere loses heat to surface
    
    % Heat eqn
    kappa = 1;
    delta = delta + kappa .* [diff([surfaceTemp; atmosTemps],2);0];
    
    dta = hittingMars - jStarMars + backRadiation;
    surfaceTemp
    backRadiation
    jStarMars
    hittingMars
    dta
    hold off
    plot([surfaceTemp; atmosTemps])
    hold on
    plot([surfaceTemp + dt .* dta; atmosTemps + dt .* delta])
    surfaceTemp = surfaceTemp + dt .* dta;
    atmosTemps = atmosTemps + dt .* delta;
    
    drawnow
end
;



%%
function pres = computePressures(Teff,p0,zs,massAvg)
    gravConst = 6.674e-11;
    marsMass = 6.39e23;
    gasConst = 8.314;
    radMars = 3.3895e6;
    kConst = gravConst .* massAvg .* marsMass ./ gasConst ./ Teff;
    pres = p0 .* exp(kConst .* (1 ./(radMars + zs) - 1./radMars));
end

function dens = computeDensities(Teff,p0,zs,massAvg)
    gasConst = 8.314;
    dens = (massAvg ./ gasConst ./ Teff) .* computePressures(Teff,p0,zs,massAvg);
end

function Intensity = computeIntensities(rho,XS,lambdas,zs,Teff)
    Intensity = zeros(size(zs)); % dE/dt function of z
    dz = diff(zs);
    dz(end+1) = dz(end);
%     dLambda = lambdas.^2;
    dLambda = 0.1;
    numLayers = size(zs,2);
    
    for zi = 1:numLayers % Index of z*
        int1 = sum(rho(1:zi).*dz(1:zi)); % Integral of rho up to z*
        expTerm = exp(-int1 .* XS);
        
        h = 6.63e-34; % J s
        c = 2.99e8; % m/s
        
        % Radiation from surface
        frontFactor = getBBSpec(Teff,lambdas) .* (pi * c);
        intensityAtAllLambda = frontFactor .* (rho(zi) .* XS) .* expTerm;
%         intensityAtAllLambda = log(frontFactor) + log(rho(zi) .* XS) + log(expTerm);
        Intensity(zi) = sum(intensityAtAllLambda .* dLambda);
        
        % Radiation from sun
        idxSunIntegral = flip(zi:size(zs,2));
        optDensIntFromSun = trapz(rho(idxSunIntegral) .* dz(idxSunIntegral));
        expTermSun = exp(-optDensIntFromSun .* XS);
        
        Asun_Amars = 8.43*10^-6; % (radius of sun/distance to mars)^2
        sunTemp = 5778; %Kelvin
        frontFactorSun = getBBSpec(sunTemp,lambdas) .* (pi * c) .* Asun_Amars;
        intensitySun = frontFactorSun .* (rho(zi) .* XS) .* expTermSun;
        
        Intensity(zi) = Intensity(zi) + sum(intensitySun .* dLambda);
        
%         IntensityPerLambdadL = zeros(size(lambdas,2),1);
%         for i = 1:(size(lambdas,2)-1)
%             IntensityPerLambdadL(i) = blackBody(i).Photons * h*c/(lambdasC(i))*N(i)*(XS(i)*rho(zi))*((lambdas(i)-lambdas(i+1))*10000);
%         end
%         Intensity(zi) = sum(IntensityPerLambdadL); % Integral over wavelength
    end

%     %copy for water vapor and liquid water as well
%     for zi = numLayers:1
%         int1 = sum(rho(numLayers:zi)*zPerLayer);
%         N = exp(-XSC*int1); %change?
%         IntensityPerLambdadL = zeros(size(lambdasC,2),1);
%         %daytimeFactor = max(sin(omega*t),0) % less light at different times of day
%         for i = 1:(size(lambdasC,2)-1)
%             IntensityPerLambdadL(i) = blackBodyCSun(i).Photons * h*c/(lambdasC(i))*N(i)*(XSC(i)*rhosC(zi))*((lambdasC(i)-lambdasC(i+1))*10000);
%         end
%         Intensity(zi) = Intensity(zi) + sum(IntensityPerLambdadL) * Asun_Amars;
%         %IntensityC(index) = IntensityC(index) + sum(IntensityPerLambdadL) *
%         %dayTimeFactor * Asun_Amars; % uncomment when doing actual simulation, takes into
%         %account daytime factor
%     end
    
end

function blackBody = getBBSpec(temp,lambdas)

    h = 6.63e-34; % J s
    c = 2.99e8; % m/s
    kb = 1.38e-23;
%     blackBody = zeros(size(lambdas));
%     for i = 1:size(lambdas,2)
%          blackBody(i) = BlackBody(temp,lambdas(i) .* 1e6).Photons;
%         blackBody(i) = h * c^3/(lambdas.^5)
         % cm to micrometer
%     end
    blackBody = pi * 2 * h * c^2 ./ (lambdas.^5) ./ (exp(h * c ./ (kb * temp .* lambdas)) - 1);
    blackBody(isnan(blackBody)) = 0;
end

% figure(4)
% plot(zs,IntensityC)
% title('Intensity')

