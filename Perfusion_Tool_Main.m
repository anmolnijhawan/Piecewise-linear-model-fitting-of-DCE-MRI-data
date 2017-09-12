%% PERFUSION TOOL MAIN FILE
%%% RAKSHIT MARCH 2016

%% INPUT PARAMETERS
global TR TE TR_Dynamic TE_Dynamic ang R1 R2 timeAIF ntpoints nuslices step_fun H W time
totalSlices = 12;  % total no of slices in data
nuslices = 12;     % remaining slices after exclusion of extreme and lower slices
ntpoints = 32;     % no of time points in data
% W = 256;           % width of image
% H = 256;           % width of image
Dynamic_images = 384; % total no of images in T1_Dynamic data folder
time = 0:0.065:2.015; % Time
step_fun = time(2)-time(1);  % Step function (in mins)  (3.9sec)

theta = 10;  % flip angle in degree
ang = cosd(theta);
R1 = 6.3/1000;  % T1_relax/1000;
R2 = 17.5/1000; % T2_relax/1000; GD-BOPAT 3T
timeAIF = 0:0.065:1.82; %% shift 3 time points to match minBAT = BAT_Ct for selective voxels(arteries)
globalAIF = 3;  % global AIF peak value for PVE correction

%% DATA READING in dicom format
% We are reading dicom images to get acquisition parameters by using dicominfo
% We can get acquisition parameters info from nifti header also, in that
% case there is no need to read dicom images
InputDir = uigetdir('Select outer data folder')
cd(InputDir)
mkdir(InputDir,'\Perfusion_Tool_Matlab_OUTPUT');  % create output folder
outpath_perfusion = strcat(InputDir,'\Perfusion_Tool_Matlab_OUTPUT\');

addpath('D:\dicom_tools'); % add path of dicom_header

[filenamePD, pathnamePD] = uigetfile('multiselect', 'on', '*.*', 'select PD w images'); % uigetfile selects multiple files from folder
[filenameT2, pathnameT2] = uigetfile('multiselect', 'on', '*.*','select T2 w images');
[filenameT1, pathnameT1] = uigetfile('multiselect', 'on','*.*', 'select T1 w images');
[filenameDynamic, pathnameDynamic] = uigetfile('multiselect', 'on','*.*', 'select DYNAMIC FSE images');

InfoFirstPD = dicom_header(fullfile(pathnamePD,filenamePD{1}));  % Info file privides information about data acquisition protocol
W = InfoFirstPD.Height; H = InfoFirstPD.Width;
TR1 = InfoFirstPD.RepetitionTime; TE1 = InfoFirstPD.EchoTime;
InfoFirstT2 = dicom_header(fullfile(pathnameT2,filenameT2{1}));
TR2 = InfoFirstT2.RepetitionTime; TE2 = InfoFirstT2.EchoTime;
InfoFirstT1 = dicom_header(fullfile(pathnameT1,filenameT1{1}));
TR3 = InfoFirstT1.RepetitionTime; TE3 = InfoFirstT1.EchoTime;
InfoFirstDynamic = dicom_header(fullfile(pathnameDynamic,filenameDynamic{1}));
TR_Dynamic = InfoFirstDynamic.RepetitionTime; TE_Dynamic = InfoFirstDynamic.EchoTime;  % Relaxation and Echo time of Dynamic data
theta = InfoFirstDynamic.FlipAngle;
TR = [TR1, TR2, TR3];  % Relaxation Times
TE = [TE1, TE2, TE3];   % Ech0 Times

SI_PD = zeros(W,H,nuslices);
SI_T2 = zeros(W,H,nuslices);
SI_T1 = zeros(W,H,nuslices);
DYNAMIC = zeros(W,H,nuslices,ntpoints);

%% Reading nifti images
% After applying some preprocessing steps (Brain extraction and Motion correction)
OuterDir = uigetdir('Select outer data folder')
cd(OuterDir)

addpath('D:\Tools for NIfTI and ANALYZE image');

[filenamePD,pathnamePD] = uigetfile('multiselect','off', '*', 'select PD data files');
Data2Read = fullfile(pathnamePD,filenamePD);
HeaderInfo = spm_vol(Data2Read);
Z1 = spm_read_vols(HeaderInfo);
for n =1:nuslices
    SI_PD_ori(:,:,n) = Z1(:,:,n)';
end
% figure, imagesc(SI_PD_ori(:,:,6)), colormap(jet), colorbar, axis image, axis off;

[filenameT2,pathnameT2] = uigetfile('multiselect','off', '*', 'select T2 data files');
Data2Read = fullfile(pathnameT2,filenameT2);
HeaderInfo = spm_vol(Data2Read);
Z2 = spm_read_vols(HeaderInfo);
for n =1:nuslices
    SI_T2_ori(:,:,n) = Z2(:,:,n)';
end
% figure, imagesc(SI_T2_ori(:,:,6)), colormap(gray), colorbar, axis image, axis off;

[filenameT1,pathnameT1] = uigetfile('multiselect','off', '*', 'select T1 data files');
Data2Read = fullfile(pathnameT1,filenameT1);
HeaderInfo = spm_vol(Data2Read);
Z3 = spm_read_vols(HeaderInfo);
for n =1:nuslices
    SI_T1_ori(:,:,n) = Z3(:,:,n)';
end
% figure, imagesc(SI_T1_ori(:,:,6)), colormap(jet), colorbar, axis image, axis off;

% Reading T1-W images binary mask created by FSL BET tool
[filenamemask,pathnamemask] = uigetfile('multiselect','off', '*', 'select T1 mask data files');
Data2Read = fullfile(pathnamemask,filenamemask);
HeaderInfo = spm_vol(Data2Read);
Z3 = spm_read_vols(HeaderInfo);
for n =1:nuslices
    SI_mask(:,:,n) = Z3(:,:,n)';
end
% figure, imagesc(SI_mask(:,:,1)), colormap(jet), colorbar, axis image, axis off;


[filenameDynamic,pathnameDynamic] = uigetfile('multiselect','on', '*', 'select Dynamic data files');
Data2Read = fullfile(pathnameDynamic,filenameDynamic);
HeaderInfo = spm_vol(Data2Read);
Z4 = spm_read_vols(HeaderInfo);
for n =1:nuslices*ntpoints
    SI_Dynamic_ori(:,:,n) = Z4(:,:,n)';
end
% figure, imagesc(SI_Dynamic_ori(:,:,1)), colormap(jet), colorbar, axis image, axis off;

% converting Dynamic data into 4-D form
SI_Dynamic = zeros(W,H,nuslices,ntpoints);
for n= 1:nuslices
    for t=1:ntpoints
        SI_Dynamic(:,:,n,t) = SI_Dynamic_ori(:,:,(n-1)*ntpoints+t);
    end
end
%  figure, imagesc(SI_Dynamic(:,:,12,32)), colormap(gray), colorbar, axis image, axis off;

% applying mask on structural images
for n=1:nuslices
    SI_T1mask(:,:,n) = SI_T1_ori(:,:,n).*SI_mask(:,:,n);
    SI_T2mask(:,:,n) = SI_T2_ori(:,:,n).*SI_mask(:,:,n);
    SI_PDmask(:,:,n) = SI_PD_ori(:,:,n).*SI_mask(:,:,n);
end
% figure, imagesc(SI_T1mask(:,:,10)), colormap(gray), axis image, axis off;

%%% Apply mask on T1 Dynamic images
for n=1:nuslices
    for t=1:ntpoints
        SI_DYNAMICmask(:,:,n,t) = (SI_Dynamic(:,:,n,t)).*SI_mask(:,:,n);
    end
end
% figure, imagesc([SI_PDmask(:,:,12), SI_DYNAMICmask(:,:,12,32)]), colormap(gray), colorbar, axis image, axis off;

%% smoothing using gaussian kernal filter
h = fspecial('gaussian', [3 3], 0.5); % higher s.d. and kernal size will make the image more blur

SI_PD = imfilter(SI_PDmask, h);
% figure, imagesc(SI_PD(:,:,6)), colormap(gray), colorbar, axis image, axis off;
SI_T2 = imfilter(SI_T2mask, h);
% figure, imagesc(SI_T2(:,:,6)), colormap(jet), colorbar, axis image, axis off;
SI_T1 = imfilter(SI_T1mask, h);
% figure, imagesc(SI_T1(:,:,6)), colormap(jet), colorbar, axis image, axis off;
DYNAMIC = imfilter(SI_DYNAMICmask, h);
% figure, imagesc(SI_DYNAMICmask(:,:,6,12)), colormap(gray), colorbar, axis image, axis off;
% figure, imagesc(DYNAMIC(:,:,6,12)), colormap(gray), colorbar, axis image, axis off;

%% Temporal Smoothing
% For temporal smoothing we were using smoothspline fitting but it is not
% so well, so this part has to improve further
%
% figure, plot(squeeze(DYNAMIC(145,150,6,:)));
% Dyfft = fft(squeeze(DYNAMIC(145,150,6,:)));
% figure, plot(time, fftshift(abs(Dyfft)),'-og');


%% T1 ESTIMATION
%%% lsqnonlin fitting for T1 calculation

% cd('E:\data_rakshit\PERFUSION TOOL MATLAB   1-Apr-2016') %change directory

parT10 = 600.0; parT1LB = 100.0; parT1UB = 5000.0; % initial guess for T1 in ms
options = optimset('Algorithm', 'levenberg-marquardt', 'MaxIter', 100, 'MaxFunEvals', 100, 'display','off');
T1map = zeros(W, H, nuslices);
warning('off');
tic
for n = 1:nuslices
    s1 =  SI_PD(:,:,n); s2 =  SI_T2(:,:,n); s3 =  SI_T1(:,:,n);
    T1_val = zeros(W,H); resu = zeros(W,H);
    
    for i=1:W
        for j =1:H
            if (s1(i,j)>0)
                %if((s1(i,j)>0) && (s2(i,j)>0) && (s3(i,j)>0))
                s(1) = s1(i,j); s(2) = s2(i,j); s(3) = s3(i,j);
                [parFit,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@FT1_calculations, parT10, parT1LB, parT1UB, options,s);
                
                T1_val(i,j) = parFit;
                resu(i,j) = residual;
            end
        end
    end
    T1_val(resu(:,:) > 0.5) = 0;
    T1map(:,:,n) = T1_val; % calculated T10 maps
    n
end
% %%% Total Time taken in T1 calculation
T10_time_min = toc/60.0  % time in minute

%  figure, imagesc(T1map(:,:,10)), colormap(gray), axis image, axis off;
% figure, imagesc(T1map(:,:,9), [100, 1500]), colormap(jet), colorbar, axis image, axis off;

%% SIGNAL INTENSITY TO CONCENTRATION TIME CURVE (SI-CTC)CONVERSION
Con0 = zeros(W,H, nuslices,ntpoints);
% Con1 = zeros(W,H, nuslices,ntpoints);
options = optimset('Algorithm', 'levenberg-marquardt', 'MaxIter', 100, 'MaxFunEvals', 30, 'display','off');
DYNAMIC_baseAv = mean(DYNAMIC(:,:,:,1:4), 4); % average of 4 pre contrast Dynamic images
% figure, imagesc(squeeze(DYNAMIC_baseAv(:,:,4,1))), colormap(gray), colorbar, axis image, axis off;
% figure, imagesc(squeeze(DYNAMIC(:,:,4,1))), colormap(gray), colorbar, axis image, axis off;

for n=1:nuslices
    for t=1:ntpoints
        for i = 1:W
            for j = 1:H
                if (DYNAMIC_baseAv(i,j,n,1)>0)
                    T10 = T1map(i,j,n);   % T10 map of slices
                    s(1) = DYNAMIC_baseAv(i,j,n,1); % mean signal intensity of 4 pre-contrast images
                    s(2) = DYNAMIC(i,j,n,t); %signal intensity of post-contrast image
                    
                    Con0(i,j,n,t) = Concentration(T10,s); %% initial guess for concentration
                    %                     C0 = squeeze(Con0(i,j,n,t));
                    %                     Con1(i,j,n,t) = fminsearch(@CTC_fmin,C0,options,T10,s); % finding a minimum of function
                end
            end
        end
    end
    n
end
Con0(Con0(:)<0)=0; Con0(Con0(:)>10)=10;  % Eliminate negative values of concentration
 figure, imagesc(Con0(:,:,8,25),[0 0.5]), colormap(jet), axis image, axis off;
%
% figure, plot(time,squeeze(Con0(208,110,10,:)),'-*r');

%% PIECE-WISE LINEAR MODEL

%%%Continious Piece Wise Linear fitting model (not generalized)
%%% Piece-wise Linear Model (PLM) fitting to calculate BAT,Beta, Slope1 and Slope2 for each voxel


%clear BAT Beta Slope1 Slope2s = zeros(5,1,1);

BAT = zeros(W,H,nuslices);
Beta = zeros(W,H,nuslices);
Slope1 = zeros(W,H,nuslices);
Slope2 = zeros(W,H,nuslices);
opts = optimset('Display','off');
warning('off')
tic
mex lm_PL_fitting.c;
for n= 10
    for i = 125
        for j = 125
            if(SI_mask(i,j,n)>0)
                batCTC = squeeze(Con0(i,j,n,:));
                batCTC = batCTC(:)';
%                 ss = lm_PL_fitting(batCTC);
%                 ss(1)= ss(1)*0.065;
%                 ss(2)= ss(2)*0.065;
%                 if(ss(6)== 1)
%                     disp(i);
%                     disp(j);
%                     disp(n);
%                 end
%                 %disp(ss(1));
                %disp(ss(2));
                    
                
                %OnsetTavg = mean(batCTC(1:6));    % Taking average of initial Onset Time points
                
                %    Initial Guess, Lower Bound, and Upper Bound for fitting parameters
                %    [ALPHA(BAT), beta(TTP), C(OnsetTime), Slop-1, slop-2]
                p0_PL = [4*step_fun; 7*step_fun; OnsetTavg; 0.5; 0];
                LB_PL = [4*step_fun; 7*step_fun; OnsetTavg-(0.06*OnsetTavg); 0; -1];
                UB_PL = [9*step_fun; 12*step_fun; OnsetTavg+(0.06*OnsetTavg); 1; 1];
                
                %if OnsetTavg>0  % remove all voxels whose initial onset time points average is zero
                    
                    [par_PL,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@PLMmain,p0_PL,time,batCTC,LB_PL,UB_PL,options);
                    
                     PL_fit = PLMmain(par_PL,time);
                     figure, plot(time,PL_fit,'-o',time,batCTC,'-*');
                    %%% PL model fitting parameters
                    %disp(ss(2));
                     BAT(i,j,n) = par_PL(1);
                    Beta(i,j,n) = par_PL(2);
                    Slope1(i,j,n) = par_PL(4);
                    Slope2(i,j,n) = par_PL(5);
                    
                    
                    
                
            end
        end
    end
    n
end
PLfitting_time_min = toc/60.0  % time in minute

figure, imagesc(BAT(:,:,12),[0,0.5]), colormap(jet), colorbar, axis image, axis off;
figure, imagesc(Beta(:,:,12),[0,0.8]), colormap(jet), colorbar, axis image, axis off;
figure, imagesc(Slope1(:,:,12),[0, 0.03]), colormap(jet), colorbar, axis image, axis off;
figure, imagesc(Slope2(:,:,12),[-0.005, 0.005]), colormap(jet), colorbar, axis image, axis off;
figure, imagesc(Slope2(:,:,12)>0, [0, 2]), colormap(jet), colorbar, axis image, axis off;
% figure, imagesc(Slope2(:,:,10))
figure, plot(squeeze(Con0(195,106,10,:)),'-o');


figure, imagesc(BAT(:,:,9)), colormap(jet), colorbar, axis image, axis off;
figure, imagesc(Beta(:,:,9)), colormap(jet), colorbar, axis image, axis off;
 figure, imagesc(Slope1(:,:,9),[0, 0.03]), colormap(jet), colorbar, axis image, axis off;
figure, imagesc(Slope2(:,:,9),[-0.05, 0.05]), colormap(jet), colorbar, axis image, axis off;
figure, imagesc(Slope2(:,:,10)>0, [0, 2]), colormap(jet), colorbar, axis image, axis off;
%figure, plot(T1, avAIF_fnl, '-Ob');
%clear BAT Beta Slope1 Slope2


CTC_PL = load('best_fit.txt');
figure, plot(time, squeeze(Con0(195,106,10,:)),'-o',time, CTC_PL,'-or');

%%% Java calculated BAT BETA comparisons
% Reading BAT
[filename, pathname] = uigetfile('multiselect','on','*.*','Select JAVA calculated BAT file');
outputpath = fullfile(pathname,filename);
fid = fopen(outputpath,'r','b');% 'b' means big endians\
nuslices = 12;
BAT_JAVA = zeros(256,256,12);
Z = fread(fid,'float32','b');
B = reshape(Z, [256, 256, 12]);
for i = 1:nuslices
    
    BAT_JAVA(:,:,i) =  (B(:,:,i).*0.065)';
end
figure, imagesc(BAT_JAVA(:,:,10),[0,0.5]), colormap(jet), colorbar, axis image, axis off;

% Reading Beta
[filename, pathname] = uigetfile('multiselect','on','*.*','Select JAVA calculated BAT file');
outputpath = fullfile(pathname,filename);
fid = fopen(outputpath,'r','b');% 'b' means big endians\
nuslices = 12;
Beta_JAVA = zeros(256,256,12);
Z = fread(fid,'float32','b');
B = reshape(Z, [256, 256, 12]);
for i = 1:nuslices
    Beta_JAVA(:,:,i) =  B(:,:,i)';
end
figure, imagesc(Beta_JAVA(:,:,10),[0, 10]), colormap(jet), colorbar, axis image, axis off;


%% AUTOMATIC AIF CALCULATION
%%% Automatic AIF calculation by calling avAIFcal function file

aifCTC = zeros(W,H,ntpoints);
BATaif = zeros(W,H);
Betaaif = zeros(W,H);
avAIF = zeros(nuslices,length(timeAIF));

for n =1:nuslices
    aifCTC = squeeze(Con0(:,:,n,:));
    aifCTC = double(aifCTC);
    %   figure, imagesc(aifCTC(:,:,32),[0 0.5]), colormap(jet), colorbar, axis image, axis off;
    BATaif = BAT(:,:,n);
    %   figure, imagesc(BAT(:,:,n)), colormap(jet), colorbar, axis image, axis off;
    Betaaif = Beta(:,:,n);
    %   figure, imagesc(Beta(:,:,n)), colormap(jet), colorbar, axis image, axis off;
    avAIF(n,:) = avAIFcal(aifCTC,BATaif,Betaaif); % avAIFcal function will estimate avg AIF for all slices
end

%% GTKM MODEL FITTING and First Pass Analysis
%%% Calling average AIF after automatic AIF calculation

selAIF = avAIF(6,:);  %select good avg AIF of any slice
% selAIF = load('seema_AIF.txt');
figure, plot(selAIF,'-*');
CTCfpass = double(Con0);
OnsetTavg = mean(selAIF(1:6));

%%% PL model fitting for avg AIF
% [ALPHA (BAT), beta (TTP), C, Slop-1, slop-2]
PLp0 = [4*step_fun; 7*step_fun; OnsetTavg; 0.5; 0];
PLLB = [4*step_fun; 7*step_fun; OnsetTavg-(0.06*OnsetTavg); 0; -1];
PLUB = [9*step_fun; 12*step_fun; OnsetTavg+(0.06*OnsetTavg); 1; 1];

[par] = lsqcurvefit(@PLMmain,PLp0,timeAIF,selAIF,PLLB,PLUB);
BAT_avAIF = round(par(1)/step_fun);  % BAT in form of time points
[a1, a2] = max(selAIF);
TPeak = a2;   % Time to Peak of avg AIF

FP1 = 2*(TPeak-BAT_avAIF);  %taking full first cycle into account
FP = BAT_avAIF:FP1;

rho = 1.04;% units of density of brain tissue in g/mL
Hem = 0.75;% Hemtocrit = (1-Hlv)/(1-Hsv) = (1.0-0.45)/(1.0-0.025)

%%% Initial guess, Lower Bound, and Upper Bound for GTKM tracer kinetic model fitting
%    [Ktrans, Kep, Vp]
p0 = [0.5, 0.6, 0.5];
LB = [0, 0, 0];
UB = [2, 4, 1];

Ktrans = zeros(W,H,nuslices);
Kep = zeros(W,H,nuslices);
Vp = zeros(W,H,nuslices);
CBV = zeros(W,H,nuslices);
CBF_R = zeros(W,H,nuslices,length(FP));
CBF = zeros(W,H,nuslices);
CBV_Corr = zeros(W,H,nuslices);
Ve = zeros(W,H,nuslices);

for n=1:nuslices
    for i = 1:W
        for j =1:H
            if(SI_mask(i,j,n)>0)
                opts = optimset('Display','off');
                CTCfpass_PL = squeeze(Con0(i,j,n,1:length(timeAIF)));
                BAT_Ct = round(BAT(i,j,n)/step_fun);
                % good avAIF selected from middle slice
                Cp = selAIF;
                
                if (BAT_Ct>0)
                    % GTKM tracer kinetic model fitting for estimation of kinetic parameters (Ktrans, Kep and Vp)
                    
                    [parm, Resnorm] = lsqcurvefit(@CtFit_FP,p0,Cp,CTCfpass_PL,LB,UB,opts,BAT_avAIF,BAT_Ct);
                    
                    Ct_fit1 = CtFit_FP(parm,Cp,BAT_avAIF,BAT_Ct);
                    % figure, plot(timeAIF,Ct_fit1,'-*',timeAIF,CTCfpass_PL,'-o');
                    
                    % Kinetic Parameters
                    Ktrans(i,j,n) = parm(1);
                    Kep(i,j,n) = parm(2);
                    Vp(i,j,n) = parm(3);
                    
                    %%% We are using simulated AIF (first gaussian of AIF given by parker) for first pass analysis
                    Cp_sim = load('simAIF.txt');
                    %figure, plot(Cp_sim/3,'-*');
                    Ct_fit = parm(3)*Cp_sim;     % Ct = Vp*Cp
                    %figure, plot(timeAIF, Ct_fit,'-*', timeAIF, Cp_sim,'-O');
                    Ve(i,j,n) = Ktrans(i,j,n)./Kep(i,j,n);  % Kep = Ktrans/Ve
                    
                    %%%%%%%%% First Pass Analysis
                    %%% CBV and CBV Corrected
                    
                    CBV_aif = sum(Cp_sim(1:length(FP))); % integral of Cp within FirstPass Int
                    CBV_Ct = sum(Ct_fit(1:length(FP))); % integralof Ct within FirstPass Int
                    % figure, plot(FP, Cp_sim(1:length(FP)),'-*', FP, Ct_fit(1:length(FP)),'-O');
                    CBV(i,j,n) = (Hem/rho)*(100*(CBV_Ct./CBV_aif));    % CBV(mL/100gm)
                    CBV_Corr(i,j,n) = CBV(i,j,n)-(Ve(i,j,n).*CBV(i,j,n)); % Corrected CBV after removal of contribution of leakage space
                    
                    %%% CBF Calculation
                    fftCT = fft(Ct_fit(1:length(FP))); %fourier transform within FirstPass
                    fftAIF = fft(Cp_sim(1:length(FP)));
                    % figure, plot(FP, fftshift(abs(fftCT)),'-*r', FP, fftshift(abs(fftAIF)),'-O');
                    CBF_R(i,j,n,1:length(FP)) = (Hem/rho)*(ifft((fftCT)./(fftAIF)));
                    % CBF_R(i,j,n,1:length(FP)) = (Hem/rho)*(ifft((fftCT)./(fftAIF)));
                    
                    CBF1(i,j,n) = max(CBF_R(i,j,n,:));
                    CBF(i,j,n) = 100*(60/3.9)*CBF1(i,j,n);  %CBF(mL/100gm/min)
                    
                end
            end
        end
        i
    end
    n
end

% figure, imagesc(Ktrans(:,:,1),[0, 0.3]),colormap(jet),colorbar, axis off, axis image;
% figure, imagesc(Kep(:,:,10)),colormap(jet), colorbar, axis off, axis image;
% figure, imagesc(Vp(:,:,10),[0, 0.08]),colormap(jet), axis off, axis image;
% figure, imagesc(Ve(:,:,10),[0, 0.5]),colormap(jet), axis off, axis image;
%
% figure, imagesc(CBF(:,:,10),[0, 70]),colormap(jet), axis off, axis image;
% figure, imagesc(CBV(:,:,10),[0, 5]),colormap(jet),colorbar, axis off, axis image;
% figure, imagesc(CBV_Corr(:,:,1), [0, 15]),colormap(gray), axis off, axis image;
% figure, plot(squeeze(CBF_R(194,106,10,:)));

%% Writing perfusion parameters maps
%OutputDir = uigetdir('select Output folder')
OutputDir = outpath_perfusion;
%cd(OutputDir)

%%% writing T10(pre-contrast T1 relaxation) maps
outputpath = fullfile(OutputDir,'T10.raw');
fid = fopen(outputpath,'w','b');% 'b' means big endians\
for i = 1:nuslices
    B = T1map(:,:,i)';
    count = fwrite(fid,B,'float32','b');
end
%%% writing Concentration(CTC) maps
mkdir(OutputDir,'\Gd_Concentration\'); savepath = strcat(OutputDir,'\Gd_Concentration\');
for n = 1:nuslices
    outputpath = fullfile(savepath,['Con',num2str(n),'.raw']);
    fid = fopen(outputpath,'w','b');
    for i = 1:ntpoints
        B = Con0(:,:,n,i)';
        count = fwrite(fid,B,'float32','b');
    end
end

%%% writing PL model fitting parameters
outputpath = fullfile(OutputDir,'BAT.raw');
fid = fopen(outputpath,'w','b');% 'b' means big endians\
for i = 1:nuslices
    B = BAT(:,:,i)';     %matrix transpose
    count = fwrite(fid,B,'float32','b');
end
outputpath = fullfile(OutputDir,'Beta.raw');
fid = fopen(outputpath,'w','b');% 'b' means big endians\
for i = 1:nuslices
    B = Beta(:,:,i)';     %matrix transpose
    count = fwrite(fid,B,'float32','b');
end
outputpath = fullfile(OutputDir,'Slope1.raw');
fid = fopen(outputpath,'w','b');% 'b' means big endians\
for i = 1:nuslices
    B = Slope1(:,:,i)';     %matrix transpose
    count = fwrite(fid,B,'float32','b');
end
outputpath = fullfile(OutputDir,'Slope2.raw');
fid = fopen(outputpath,'w','b');% 'b' means big endians\
for i = 1:nuslices
    B = Slope2(:,:,i)';     %matrix transpose
    count = fwrite(fid,B,'float32','b');
end

%%% Writing kinetic parameters
outputpath = fullfile(OutputDir,'Ktrans.raw');
fid = fopen(outputpath,'w','b');% 'b' means big endians\
for i = 1:nuslices
    B = Ktrans(:,:,i)';     %matrix transpose
    count = fwrite(fid,B,'float32','b');
end
outputpath = fullfile(OutputDir,'Kep.raw');
fid = fopen(outputpath,'w','b');
for i = 1:nuslices
    B = Kep(:,:,i)';
    count = fwrite(fid,B,'float32','b');
end
outputpath = fullfile(OutputDir,'Vp.raw');
fid = fopen(outputpath,'w','b');
for i = 1:nuslices
    B = Vp(:,:,i)';
    count = fwrite(fid,B,'float32','b');
end
outputpath = fullfile(OutputDir,'Ve.raw');
fid = fopen(outputpath,'w','b');
for i = 1:nuslices
    B = Ve(:,:,i)';
    count = fwrite(fid,B,'float32','b');
end

%%% Writing hemodynamic parameters
outputpath = fullfile(OutputDir,'CBF.raw');
fid = fopen(outputpath,'w','b');
for i = 1:nuslices
    B = CBF(:,:,i)';
    count = fwrite(fid,B,'float32','b');
end
outputpath = fullfile(OutputDir,'CBV.raw');
fid = fopen(outputpath,'w','b');
for i = 1:nuslices
    B = CBV(:,:,i)';
    count = fwrite(fid,B,'float32','b');
end
outputpath = fullfile(OutputDir,'CBV_Corr.raw');
fid = fopen(outputpath,'w','b');
for i = 1:nuslices
    B = CBV_Corr(:,:,i)';
    count = fwrite(fid,B,'float32','b');
end