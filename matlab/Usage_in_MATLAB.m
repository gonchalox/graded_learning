%% This Script demonstrates how to correctly load the HDR-Images from the HdM-HDR-2014 data set.

%% Init 
% Adjust this Path to the Path where you downloaded the HdM-HDR-2014 data set
ImgPath = '/media/gluzardo/Data/Stuttgart/'; 

ImgName = {'beerfest_lightshow_01','beerfest_lightshow_02','beerfest_lightshow_03',...
    'beerfest_lightshow_04','beerfest_lightshow_05','beerfest_lightshow_06',...
    'beerfest_lightshow_07','bistro_01','bistro_02','bistro_03','carousel_fireworks_01',...
    'carousel_fireworks_02','carousel_fireworks_03','carousel_fireworks_04',...
    'carousel_fireworks_05','carousel_fireworks_06','carousel_fireworks_07',...
    'carousel_fireworks_08','carousel_fireworks_09','cars_closeshot','cars_fullshot',...
    'cars_longshot','fireplace_01','fireplace_02','fishing_closeshot','fishing_longshot',...
    'hdr_testimage','poker_fullshot','poker_travelling_slowmotion','showgirl_01',...
    'showgirl_02','smith_hammering','smith_welding'};

StartFrame = [1591; 1225; 1702; 3467; 4953; 6101; 6951; 295; 237;  858; 1187;  936; 348; 115;...
    756; 5034; 128;   1; 3715;  599; 350;  258;  660; 319; 354;  420; 1020; 370; 1034; 183;...
    348; 448; 413];

EndFrame =   [1684; 1412; 1791; 4327; 5048; 6417; 7142; 445; 848; 1027; 1625; 1193; 533; 322;...
    929; 5561; 272; 339; 4242; 1012; 791; 1077; 1150; 779; 723; 1253; 1499; 969; 2980; 958;...
    688; 914; 1514];

AlexaWideGamut2sRGB = [1.617523436306807  -0.070572740897816  -0.021101728042793;...
  -0.537286622188294   1.334613062330328  -0.226953875218266;...
  -0.080236814118512  -0.264040321432512   1.248055603261060];

sRGBDeLinearize = @(x)((x>0.0031308).*(1.055.*x.^(1/2.4)-0.055)+(x<=0.0031308).*12.92.*x);

%% Load a Frame:
% Select an Image from the last sequence
ImgSequence = 20;
% Select Frame 556 from the range of 413-1514
fcc = 599;     

% Read OpenEXR-File with exrread for MATLAB, downloaded from: http://www.mit.edu/~kimo/software/matlabexr/
Img = exrread([ImgPath,cell2mat(ImgName(ImgSequence)),'/',cell2mat(ImgName(ImgSequence)),'_',...
    num2str(fcc,'%06d'),'.exr']);

%% Display - Play around with exposure values between -10 and 5:
Exposure = -10;

linearsRGBImg = reshape(reshape((Img.*2^Exposure),[],3)*AlexaWideGamut2sRGB,1080,1920,3);
imshow(sRGBDeLinearize(linearsRGBImg));

%% Gratulation, you have sucessfully loaded the HdM-HDR-2014 footage. 
% Good luck finding new Algorithms in HDR-processing and display. 
% If our images helped you, please cite our paper.
%
%% Title: 
% Creating cinematic wide gamut HDR-video for the evaluation of tone mapping operators and HDR-displays
%% Authors: 
% Jan Froehlich, Stefan Grandinetti, Bernd Eberhardt, Simon Walter, Andreas Schilling, Harald Brendel
%% Presented at the:
% SPIE/IS&T Electronic Imaging Conference, San Francisco, February the 5th 2014
%
%% Abstract: 
% High quality video sequences are required for the evaluation of tone mapping operators and high 
% dynamic range (HDR) displays. We provide scenic and documentary scenes with a dynamic range of up 
% to 18 stops. The scenes are staged using professional film lighting, make-up and set design to 
% enable the evaluation of image and material appearance. To address challenges for HDR-displays and
% temporal tone mapping operators, the sequences include highlights entering and leaving the image,
% brightness changing over time, high contrast skin tones, specular highlights and bright, saturated
% colors. HDR-capture is carried out using two cameras mounted on a mirror-rig. To achieve a 
% cinematic depth of field, digital motion picture cameras with Super-35mm size sensors are used. 
% We provide HDR-video sequences to serve as a common ground for the evaluation of temporal tone 
% mapping operators and HDR-displays. They are available to the scientific community for further 
% research.
%% Keywords:
% High Dynamic Range, HDR-Video, Wide Gamut, Tone Mapping.