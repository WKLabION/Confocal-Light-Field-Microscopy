clear all;
close all;
addpath("./subfunctions/");
addpath("./example_data");

%%
FilePath=['.\']; %the path of the zip file
FileNamePrefix=['example_data\example'];
ROIFileNamePrefix=[FilePath 'roi_example'];

%%
Param.SampleRate=400;           % image frame rate (unit Hz)
Param.SpikeTemplateLength=8;    % time window for template generation (time points), time window of [-TemplateLength:TemplateLength] is used for template generation
Param.SpikeTemplateN=100;        % number of peaks used for template generation    
Param.CutOffFreq=1;             % the cutoff frequency to highpass filter the raw trace;
Param.SpikePolarity=-1;         % -1 for negtively polarized spikes
Param.MinSpikeTemplateN=[8 4 2 1];      % minimum number of spikes to calcuate template, if the number of this case is small, it is not counted as spiking neuron.
Param.SNRList=[3.5 4 4.5 5];
Param.HeadTailSize=200;
Param.CellEnvSize=15;           % the size extending from soma to estimate local background to remove background noise and interferences


Tmp=load([FilePath  'pca.mat'],'NoisePCA');
NoisePCA=Tmp.NoisePCA;

    %%
load([FilePath FileNamePrefix '.mat']);
ImStack=noise_pca_crt(ImStackCrt,NoisePCA);
CrtImg=ones(size(ImStackCrt)).*mean(ImStackCrt,3)+double(ImStack);
imstackwrite(uint16(CrtImg),[FilePath FileNamePrefix  '_crt.tif']);
    %%
ImPreview=std(single(ImStack(:,:,1:min(1000,size(ImStack,3)))),0,3);
ImPreview=ImPreview/max(ImPreview(:));
ImPreview=repmat(ImPreview,1,1,3);
    
load([ROIFileNamePrefix '.mat']);

ROITot=ImPreview(:,:,1)*0;
for ii=1:length(ROIs)
     ROITot(ROIs{ii})=1;
end
figure(1);imagesc(ROITot);axis image;
    
ROIMask=ImPreview(:,:,1)*0;
    
Neuron=[];
for ii=1:length(ROIs)
     Neuron(ii).ROI=ROIs{ii};
     Neuron(ii).SpikeInfo=[];
     if ~isempty(Neuron(ii).ROI)
        disp(['roi # : ' num2str(ii)]);
    
        ROIMask=ROIMask*0;
        ROIMask(Neuron(ii).ROI)=1;
    
        ImPreview=display_roi(ImPreview,ROIMask);
        figure(1);imagesc(ImPreview);axis image;title(['roi # ' num2str(ii)]);
    
        [DataROI SomaIdx BkgIdx]=extract_roi(Param,ImStack,ROIMask,ROITot);
        SpikeInfo=spike_extract(Param,DataROI,SomaIdx,BkgIdx);
        Neuron(ii).SpikeInfo=SpikeInfo;
    
%         pause(0.1);
            pause;
    
    end
end
    
NumberSpikeNeuron=count_spike_neuron(Neuron)
save([FilePath  'spikeinfo.mat'],'Neuron','Param');

    
 %%   
    
function NumberSpikeNeuron=count_spike_neuron(Neuron)
Count=0;
for ii=1:size(Neuron,2)
    if ~isempty(Neuron(ii).SpikeInfo)
        Count=Count+1;
    end
end
NumberSpikeNeuron=Count;
end

function ImgCrt=noise_pca_crt(ImgRaw,NoisePCA)
NoiseCrtN=6;
Trace=NoisePCA.Trace(1:NoiseCrtN,:);

ImgSize=size(ImgRaw);
ImgRaw=reshape(single(ImgRaw),ImgSize(1)*ImgSize(2),ImgSize(3));
Coef=inv(Trace*Trace')*(Trace*ImgRaw');
ImgCrt=ImgRaw-Coef'*Trace;
ImgCrt=uint16(ImgCrt-min(ImgCrt(:)));
ImgCrt=reshape(ImgCrt,ImgSize(1),ImgSize(2),ImgSize(3));
end

function Img=display_roi(Img,ROIMask)
ROIMask=repmat(ROIMask,1,1,3);
ROIMask=ROIMask.*rand(1,1,3)/6;
Img=Img+ROIMask;
end

function [DataROI1 SomaIdx BkgIdx]=extract_roi(Param,ImStack1,ROIMask,ROITot)
SpikePolarity=Param.SpikePolarity;
CellEnvSize=Param.CellEnvSize;

Tmp=find(max(ROIMask,[],1));
xmin=min(Tmp);
xmax=max(Tmp);
Tmp=find(max(ROIMask,[],2));
ymin=min(Tmp);
ymax=max(Tmp);

xmin=max(1,xmin-CellEnvSize);
xmax=min(size(ROIMask,2),xmax+CellEnvSize);
ymin=max(1,ymin-CellEnvSize);
ymax=min(size(ROIMask,1),ymax+CellEnvSize);

DataROI1=double(ImStack1(ymin:ymax,xmin:xmax,:));

SomaMask=ROIMask(ymin:ymax,xmin:xmax);
BkgMask=1-ROITot(ymin:ymax,xmin:xmax);
% BkgMask=imgaussfilt(double(SomaMask),3);
% BkgMask=1-(BkgMask>(max(BkgMask(:))/20));
DataROI1=DataROI1*SpikePolarity;
SomaIdx=find(SomaMask(:));
BkgIdx=find(BkgMask(:));
figure(5000);subplot(1,2,1);imagesc(SomaMask);axis image;
figure(5000);subplot(1,2,2);imagesc(BkgMask);axis image;

end



function SpikeInfo=spike_extract(Param,DataROI1,SomaIdx,BkgIdx)

SomaMask=DataROI1(:,:,1)*0;
SomaMask(SomaIdx)=1;

DataROIFilt=highpass_video_filt(Param,DataROI1);
DataROIFilt=reshape(DataROIFilt,[],size(DataROIFilt,3));
RawTrace=mean(DataROIFilt(SomaIdx,:),1);
RawTrace=RawTrace(:);
Bkg=mean(DataROIFilt(BkgIdx,:),1);
Bkg=Bkg(:);
Trace=RawTrace-Bkg;
Trace1=Trace-mean(Trace(:));

Mask=(Trace1<0);
NoiseAmp1=sqrt(sum(Trace1.^2.*Mask)/sum(Mask));

Trace=Trace1;

SpikeInfo=spike_denoise(Param,Trace);

figure(1000);ax1=subplot(3,1,1);plot(Trace);
figure(1001);subplot(2,2,1);imagesc(std(DataROI1,0,3));axis image;

if ~isempty(SpikeInfo)
    figure(1000);
    ax1=subplot(3,1,1);
    hold on;plot(SpikeInfo.SpikeIdx,ones(length(SpikeInfo.SpikeIdx),1)*max(Trace(:))*1.1,'go');
    plot(SpikeInfo.SpikeTemplateIdx,ones(length(SpikeInfo.SpikeTemplateIdx),1)*max(Trace(:))*1.15,'ro');hold off;
    ax2=subplot(3,1,2);
    plot(SpikeInfo.RawTrace);
    ax3=subplot(3,1,3);
    plot(SpikeInfo.FiltTrace);
    linkaxes([ax1 ax2 ax3],'x');

    figure(1001);
    subplot(2,2,3);imagesc(SomaMask);axis image;
    subplot(2,2,4);plot(SpikeInfo.SpikeTemplate);
else
    figure(1000);
    ax2=subplot(3,1,2);
    plot(Trace*0);
    ax3=subplot(3,1,3);
    plot(Trace*0);
    linkaxes([ax1 ax2 ax3],'x');

    figure(1001);
    subplot(2,2,3);imagesc(SomaMask);axis image;
    subplot(2,2,4);plot([0 0]);
end

end

function VideoFilt=highpass_video_filt(Param,Video)
CutOffFreq=Param.CutOffFreq;
SampleRate=Param.SampleRate;

NormFreq=CutOffFreq/(SampleRate/2);
[bb,aa]=butter(3,NormFreq, 'high');
Video=permute(Video,[3 1 2]);
VideoFilt=filtfilt(bb,aa,Video); 
VideoFilt=permute(VideoFilt,[2 3 1]);
end