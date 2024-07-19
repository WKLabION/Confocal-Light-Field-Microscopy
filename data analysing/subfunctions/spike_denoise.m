function SpikeInfo=spike_denoise(Param,RawTrace)

SampleRate=Param.SampleRate;
SpikeTemplateLength=Param.SpikeTemplateLength;
SpikeTemplateN=Param.SpikeTemplateN;
CutOffFreq=Param.CutOffFreq;
HeadTailSize=Param.HeadTailSize;
SNRList=Param.SNRList;
MinSpikeTemplateN=Param.MinSpikeTemplateN;

%%
DataHP=highpass_filt(RawTrace,CutOffFreq,SampleRate);
[SpikeTemplateIdx SNRThd PosNegPeakCnt]=template_threshold(DataHP,HeadTailSize,SNRList,MinSpikeTemplateN);
if ~isempty(SpikeTemplateIdx)
    [SpikeIdx RawTraceNoiseAmp]=simple_threshold(DataHP,SNRThd);

    SpikeTemplate=spike_template_gen(DataHP,SpikeTemplateIdx,SpikeTemplateLength,SpikeTemplateN);
    DataFilt=whitened_matched_filter(DataHP,SpikeTemplateIdx, SpikeTemplateLength,SpikeTemplateN);
    [SpikeIdx FiltTraceNoiseAmp]=simple_threshold(DataFilt,SNRThd);
    SpikeRecon=DataFilt*0;
    SpikeRecon(SpikeIdx)=1;
    SpikeRecon=conv(SpikeRecon,SpikeTemplate,'same');
    Coef=sum(SpikeRecon.*RawTrace)/sum(SpikeRecon.^2);
    SpikeSub=RawTrace-Coef*SpikeRecon;

    SpikeInfo.RawTrace=RawTrace;
    SpikeInfo.FiltTrace=DataFilt;

    SpikeInfo.SpikeIdx=SpikeIdx;
    SpikeInfo.SpikeTemplate=SpikeTemplate;
    SpikeInfo.SpikeRecon=SpikeRecon;
    SpikeInfo.SpikeSub=SpikeSub;
    SpikeInfo.SpikeTemplateIdx=SpikeTemplateIdx;
    SpikeInfo.SNRThd=SNRThd;
    SpikeInfo.PosNegPeakCnt=PosNegPeakCnt;

    figure(2000);
    ax1=subplot(3,1,1);plot(RawTrace);
    ax2=subplot(3,1,2);plot(DataFilt);
    ax3=subplot(3,1,3);plot(SpikeRecon);
    linkaxes([ax1 ax2 ax3],'x');
else
    SpikeInfo=[];
end



end

%%
function DataOut = whitened_matched_filter(DataIn,SpikeIdx, SpikeTemplateLength,SpikeTemplateN)
SpikeRecon=DataIn*0;
SpikeRecon(SpikeIdx)=1;
SpikeMask=conv(SpikeRecon,ones(2*SpikeTemplateLength+1,1),'same');
NoiseMask=SpikeMask<0.5;

Noise=DataIn(find(NoiseMask));
NoisePSD=pwelch(Noise,1000,[],2*length(DataIn)-1);
NoisePSD=[NoisePSD(:);flipud(NoisePSD(1:end-1))];
DataWhiten=real(ifft(fft(DataIn,length(NoisePSD))./sqrt(NoisePSD)));
DataWhiten=DataWhiten(1:length(DataIn));
SpikeTemplate=spike_template_gen(DataWhiten,SpikeIdx,SpikeTemplateLength,SpikeTemplateN);
DataOut=conv(DataWhiten,flipud(SpikeTemplate),'same');
end

function SpikeTemplate=spike_template_gen(Data,SpikeIdx,SpikeTemplateLength,SpikeTemplateN)
SpikeIdxMask=(SpikeIdx>SpikeTemplateLength).*(SpikeIdx<(length(Data)-SpikeTemplateLength));
SpikeIdx=SpikeIdx(find(SpikeIdxMask));

if length(SpikeIdx)>SpikeTemplateN
    [~,TmpIdx]=sort(Data(SpikeIdx),'descend');
    SpikeIdx=SpikeIdx(TmpIdx(1:SpikeTemplateN));
end

SpikeTemplate=[];
for ii=1:length(SpikeIdx)
    SpikeTemplate(:,ii)=Data(SpikeIdx(ii)-SpikeTemplateLength:SpikeIdx(ii)+SpikeTemplateLength);
end

% figure(100);
% plot(SpikeTemplate,'color', [0.5 0.5 0.5]);hold on;
SpikeTemplate=mean(SpikeTemplate,2);
% plot(SpikeTemplate,'k','linewidth',2);hold off;
end


function [PeakIdx NoiseAmp]=simple_threshold(Data,SNRThd)
Mask=Data<0;
NoiseAmp=sqrt(sum(Mask.*Data.^2)/sum(Mask));
Data=Data.*(Data>(SNRThd*NoiseAmp));
[~,PeakIdx]=findpeaks(Data);
end

function [PeakTemplateIdx SNRThd PosNegPeakCnt]=template_threshold(Data,HeadTailSize,SNRList,MinSpikeTemplateN)
%% spike has to be a single peak with 2 preceding and following data points smaller than the peak

Data=Data(HeadTailSize+1:end-HeadTailSize);
Mask=Data<0;
NoiseAmp=sqrt(sum(Mask.*Data.^2)/sum(Mask));
Data=Data/NoiseAmp;
PosPeaks=(Data>circshift(Data,1)).*(Data>circshift(Data,-1));
NegPeaks=(Data<circshift(Data,1)).*(Data<circshift(Data,-1));
PosPeakIdx=find(PosPeaks);
NegPeakIdx=find(NegPeaks);

SNRThdTheory=1/2*(1+erf(-SNRList/sqrt(2)))*length(Data(:));

SinglePeaks=(Data>circshift(Data,1)).*(Data>circshift(Data,2)).*(Data>circshift(Data,-1)).*(Data>circshift(Data,-2));
SinglePeakIdx=find(SinglePeaks);

PosNegPeakCnt=[];
PeakTemplateIdx={};
for ii=1:length(SNRList)
    PosNegPeakCnt(ii,1)=sum(Data(PosPeakIdx)>SNRList(ii));
    PosNegPeakCnt(ii,2)=sum(Data(NegPeakIdx)<-SNRList(ii));
end
Tmp1=PosNegPeakCnt(:,1)>(max(PosNegPeakCnt(:,2),SNRThdTheory(:))*3);
Tmp2=(PosNegPeakCnt(:,1)-PosNegPeakCnt(:,2))>MinSpikeTemplateN(:);
SNRThd=SNRList(find(Tmp1.*Tmp2,1));
PosNegPeakCnt

PeakTemplateIdx=[];
if ~isempty(SNRThd)
    SinglePeaks=(Data>circshift(Data,1)).*(Data>circshift(Data,2)).*(Data>circshift(Data,-1)).*(Data>circshift(Data,-2));
    SinglePeakIdx=find(SinglePeaks);
    TmpIdx=find(Data(SinglePeakIdx)>SNRThd);
    PeakTemplateIdx=SinglePeakIdx(TmpIdx)+HeadTailSize;
end
end

function DataOut=highpass_filt(DataIn,CutOffFreq,SampleRate)
NormFreq=CutOffFreq/(SampleRate/2);
[bb,aa]=butter(2,NormFreq,'high');
DataOut=filtfilt(bb,aa,double(DataIn));
DataOut=DataOut-median(DataOut(:));
end