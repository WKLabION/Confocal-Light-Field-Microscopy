clear all;
close all;
addpath("./subfunctions/");
addpath("./example_data");

%%
FilePath=['.\']; %the path of the zip file
FileNamePrefix=['example_data\'];
FileNameSuffix=['spikeinfo.mat'];

Layer=load([FilePath  FileNameSuffix]);

Preview=imread([FilePath FileNamePrefix 'STD_example.tif']);


Data.Param=Layer(1).Param;
Data.Layer=Layer;
Data.Preview=Preview;

%%
Data.ActiveIdxList=[];
Data.RawTraceMatrix=[];
Data.FiltTraceMatrix=[];
Data.SpikeReconMatrix=[];
Data.SpikeIdxMatrix={};
Data.Groups=[];

kk=1;
for ii=1:length(Layer)
    NeuronTmp=Layer(ii).Neuron;
    for jj=1:length(NeuronTmp)
        if ~isempty(NeuronTmp(jj).SpikeInfo)
            Data.ActiveIdxList(end+1,:)=[ii jj];
            Data.RawTraceMatrix(:,end+1)=NeuronTmp(jj).SpikeInfo.RawTrace(:);
            Data.FiltTraceMatrix(:,end+1)=NeuronTmp(jj).SpikeInfo.FiltTrace(:);
            Data.SpikeReconMatrix(:,end+1)=NeuronTmp(jj).SpikeInfo.SpikeRecon(:);
            Data.SpikeIdxMatrix{kk}=NeuronTmp(jj).SpikeInfo.SpikeIdx(:);
            kk=kk+1;
        end
    end
end

%% clustering into groups
Cutoff=0.8;
Groups=neuron_cluster(Data,Cutoff);
Data.Groups=Groups;
% pause;


%% manual check grouping results

GroupCheck=neuron_cluster_check(Data);
Data.GroupCheck=GroupCheck;

%%
MatrixUpdate=neuron_combine_maxsnr(Data);
Data.ActiveIdxListCm=MatrixUpdate.ActiveIdxListCm;

save([FilePath 'result.mat'],'Data');

%%
function MatrixUpdate=neuron_combine_maxsnr(Data)
GroupCheck=Data.GroupCheck;
Layer=Data.Layer;
Param=Data.Param;

ActiveIdxList=Data.ActiveIdxList;
ActiveIdxListCm=[];

for ii=1:max(GroupCheck)
    ii
    TmpIdx=find(GroupCheck==ii);
    if ~isempty(TmpIdx)
        NeuronIdx=ActiveIdxList(TmpIdx,:);
        RawTraceTotDisp=[];
        SNR=[];

        SpikeN=[];
        for jj=1:size(NeuronIdx,1)
            Neuron=Layer(NeuronIdx(jj,1)).Neuron(NeuronIdx(jj,2));
            SpikeN(jj)=length(Neuron.SpikeInfo.SpikeIdx);
        end

        for jj=1:size(NeuronIdx,1)
            Neuron=Layer(NeuronIdx(jj,1)).Neuron(NeuronIdx(jj,2));
            FiltTrace=Neuron.SpikeInfo.FiltTrace(:);
            SpikeIdx=Neuron.SpikeInfo.SpikeIdx;
            Mask=FiltTrace<0;
            NoiseAmp=sqrt(sum(Mask.*FiltTrace.^2)/sum(Mask));
            SpikeAmp=FiltTrace(SpikeIdx);
            SpikeAmp=sort(SpikeAmp,'descend');
            SNR(jj)=mean(SpikeAmp(1:min(SpikeN)))/NoiseAmp;


            RawTraceTmp=Neuron.SpikeInfo.RawTrace(:);
            RawTraceTotDisp(:,jj)=RawTraceTmp/max(RawTraceTmp(:))+jj*1.5;
        end
        figure(3001);plot(RawTraceTotDisp);

        figure(3000);subplot(1,2,1);plot(SNR);
        figure(3000);subplot(1,2,2);plot(SpikeN);

        [~,IdxSNRMax]=max(SNR);
        IdxSNRMax
        ActiveIdxListCm(end+1,:)=ActiveIdxList(TmpIdx(IdxSNRMax),:);
        pause(1);
    end
end

MatrixUpdate.ActiveIdxListCm=ActiveIdxListCm;

end

%%

function GroupCheck=neuron_cluster_check(Data)
Groups=Data.Groups;
GroupCheck=Groups;
SpikeReconMatrix=Data.SpikeReconMatrix;
RawTraceMatrix=Data.RawTraceMatrix;
FiltTraceMatrix=Data.FiltTraceMatrix;
ActiveIdxList=Data.ActiveIdxList;
Layer=Data.Layer;
Preview=Data.Preview;

for ii=1:max(Groups)
    ii

    TmpIdx=find(Groups==ii);
    SpikeReconDisp=SpikeReconMatrix(:,TmpIdx);
    RawTraceDisp=RawTraceMatrix(:,TmpIdx);
    FiltTraceDisp=FiltTraceMatrix(:,TmpIdx);
    SpikeReconDisp=SpikeReconDisp./max(SpikeReconDisp,[],1)+[1:length(TmpIdx)];
    RawTraceDisp=RawTraceDisp./max(RawTraceDisp,[],1)+[1:length(TmpIdx)];
    FiltTraceDisp=FiltTraceDisp./max(FiltTraceDisp,[],1)+[1:length(TmpIdx)];

    figure(1);plot(RawTraceDisp);title('Raw Trace')
    figure(2);plot(FiltTraceDisp);title('Filted Trace')
    figure(3);plot(SpikeReconDisp);title('Spike Reconstructed Trace')

    NeuronIdx=ActiveIdxList(TmpIdx,:);
    neuron_disp(NeuronIdx,Layer,Preview);

    Reply=input('more than one neuron?:','s');
    if Reply=='y'
        Reply=input('input the neuron groups: ');
        while length(Reply)~=length(TmpIdx)
            Reply=input('wrong input, please input the neuron groups');
        end
        NewGroup=[];
        for jj=1:length(Reply)
            NewGroup(jj)=Reply(jj);
        end
        NewGroup
        TmpIdx1=find(NewGroup==0);
        GroupCheck(TmpIdx(TmpIdx1))=0;
        for jj=2:max(NewGroup)
            TmpIdx1=find(NewGroup==jj);
            GroupCheck(TmpIdx(TmpIdx1))=max(GroupCheck(:))+1;
        end
    end
end
end

%%
function Groups=neuron_cluster(Data,Cutoff)
M=Data.SpikeReconMatrix;
CM=corrcoef(M);
figure(1);
hist(CM(:),50);
Dist=1-CM;
figure(2);imagesc(CM);
Z=linkage(squareform(Dist),'average');
Groups=cluster(Z,'cutoff',Cutoff,'criterion','distance');
max(Groups)

[TmpTmp IndexTmp]=sort(Groups);
OBG(:,[1:length(Groups)])=M(:,IndexTmp);
CMG=corrcoef(OBG);
figure(3);imagesc(CMG);axis image;
end


%%
function neuron_disp(NeuronIdx,Layer,Preview)
PreviewDisp=Preview(:,:,1)*0;
CellEnvSize=15;
ROITmp=PreviewDisp*0;
for ii=1:size(NeuronIdx,1)
    Neuron=Layer(NeuronIdx(ii,1)).Neuron(NeuronIdx(ii,2));
    ROITmp=ROITmp*0;
    ROITmp(Neuron.ROI)=1;
    ROIStack(:,:,ii)=ROITmp;
 end
PreviewDisp=sum(ROIStack,3);
figure(1000); imagesc(PreviewDisp);axis image;
for ii=1:size(NeuronIdx,1) 
    Neuron=Layer(NeuronIdx(ii,1)).Neuron(NeuronIdx(ii,2));
    [x,y]=ind2sub(size(PreviewDisp),Neuron.ROI(1));
    text(y,x,num2str(ii));
end


Tmp=find(max(PreviewDisp,[],1));    
xmin=min(Tmp);
xmax=max(Tmp);
Tmp=find(max(PreviewDisp,[],2));
ymin=min(Tmp);
ymax=max(Tmp);

xmin=max(1,xmin-CellEnvSize);
xmax=min(size(PreviewDisp,2),xmax+CellEnvSize);
ymin=max(1,ymin-CellEnvSize);
ymax=min(size(PreviewDisp,1),ymax+CellEnvSize);

ImgDisp=[];
SNRDisp=[];
for ii=1:size(NeuronIdx,1)
    Neuron=Layer(NeuronIdx(ii,1)).Neuron(NeuronIdx(ii,2));
    ImgDisp(:,:,ii)=Preview(ymin:ymax,xmin:xmax,NeuronIdx(ii,1));
    ImgDisp(:,:,ii)=ImgDisp(:,:,ii)+double(ROIStack(ymin:ymax,xmin:xmax,ii))*(max(max(ImgDisp(:,:,ii)))/10);
    Peaks=Neuron.SpikeInfo.FiltTrace(Neuron.SpikeInfo.SpikeIdx);
    [Tmp1,Tmp2]=histcounts(Peaks,[0:0.5:10]);
    SNRDisp(ii,:)=Tmp1;
end
ImgDisp=reshape(ImgDisp,size(ImgDisp,1),[]);
figure(2000);imagesc(ImgDisp);axis image;title(num2str(NeuronIdx(:,1)'));
% figure(3000);plot(SNRDisp');legend;
end
