clear all;
close all;
addpath("./subfunctions/");
addpath("./example_data");

FilePath=['.\']; %the path of the zip file
FileNamePrefix=['example_data\example'];
FileNameSuffix=['.tif'];
BNxy=4;

%% generate motion correction template and estimate motion shifts
ImStack=imstackread([FilePath FileNamePrefix FileNameSuffix]);
Template=mean(single(ImStack),3);

%%
ROIPos=[1 1 size(ImStack)]; %[y0 x0 y_width x_width]
Shift=motion_est(Template,ROIPos,ImStack);
save([FilePath  'motion_estimation.mat'],'Template','ROIPos','Shift');

%%

ImStackCrt=motion_crt(ImStack,Shift);
ImStackCrt=ImStackCrt(:,:,2:end)+ImStackCrt(:,:,1:end-1);
save([FilePath FileNamePrefix '.mat'],'ImStackCrt','-v7.3');

ImPCA=[];
Tmp=ImStackCrt(1:floor(size(ImStackCrt,1)/BNxy)*BNxy,1:floor(size(ImStackCrt,2)/BNxy)*BNxy,:);
Tmp=reshape(Tmp,BNxy,size(Tmp,1)/BNxy,BNxy,size(Tmp,2)/BNxy,[]);
ImPCA(:,:,:)=squeeze(mean(mean(single(Tmp),1),3));
%%
NoisePCA=noise_pca_extract(ImPCA);
save([FilePath  'pca.mat'],'NoisePCA','BNxy','ImPCA');

%%
function NoisePCA=noise_pca_extract(Img)
CmpN=50;
Img=permute(Img,[1 2 4 3]);
Img=reshape(Img,size(Img,1),[],size(Img,4));
ImgSize=size(Img)
Img=reshape(Img,ImgSize(1)*ImgSize(2),[]);
%%
tic
[U S V]=svds(Img,CmpN);
toc

%%
Patterns=reshape(U,ImgSize(1),ImgSize(2),[]);
Trace=S*V';
clear S V

%%
for ii=1:size(Patterns,3)
    figure(1);subplot(1,2,1);imagesc(Patterns(:,:,ii));axis image;colorbar;
    figure(1);subplot(1,2,2);plot(Trace(ii,:));
    title(num2str(ii));
%     pause;
end
NoisePCA.Patterns=Patterns;
NoisePCA.Trace=Trace;
end


function Shift=motion_est(Template,ROIPos,ImStack)
%%
MaxShift=4;

ROI=ImStack(ROIPos(1):ROIPos(1)+ROIPos(3)-1,ROIPos(2):ROIPos(2)+ROIPos(4)-1,:);
ROI=double(ROI);
TemplateROI=Template(ROIPos(1):ROIPos(1)+ROIPos(3)-1,ROIPos(2):ROIPos(2)+ROIPos(4)-1,:);
Mask=ROI(:,:,1)*0;
Mask(MaxShift+1:end-MaxShift,MaxShift+1:end-MaxShift)=1;

ROI=(ROI-circshift(ROI,[1 0 0]))+j*(ROI-circshift(ROI,[0 1 0]));
TemplateROI=(TemplateROI-circshift(TemplateROI,[1 0]))+j*(TemplateROI-circshift(TemplateROI,[0 1]));


figure(1);imagesc(Mask);axis image;
figure(2);imagesc(abs(TemplateROI));axis image;


%%
ROIRaw=ROI.*Mask;
TemplateROI=TemplateROI.*Mask;

XCorr=abs(fftshift(fftshift(ifft2(fft2(TemplateROI).*conj(fft2(ROIRaw))),1),2));
Ref=sqrt(abs(fftshift(fftshift(ifft2(fft2(Mask).*conj(fft2(abs(ROIRaw).^2))),1),2)));
Ref=Ref.*sqrt(abs(fftshift(fftshift(ifft2(fft2(abs(TemplateROI).^2).*conj(fft2(Mask))),1),2)));
XCorr=XCorr./(Ref+eps);
Center=floor(size(XCorr)/2+1);
XCorrROI=XCorr(Center(1)-MaxShift:Center(1)+MaxShift,Center(2)-MaxShift:Center(2)+MaxShift,:);
Tmp=reshape(XCorrROI,[],size(XCorrROI,3));
[XCorrMax Idx]=max(Tmp,[],1);
[YIdx XIdx]=ind2sub([size(XCorrROI,1) size(XCorrROI,2)],Idx);

PeakFitSize=2;
[xx yy]=meshgrid([-PeakFitSize:1:PeakFitSize],[-PeakFitSize:1:PeakFitSize]);
for ii=1:size(XCorr,3)
    Tmp=XCorrROI(YIdx(ii)-PeakFitSize:YIdx(ii)+PeakFitSize,XIdx(ii)-PeakFitSize:XIdx(ii)+PeakFitSize,ii);
    Tmp=Tmp-min(Tmp(:));
    YIdxFine(ii)=sum(sum(yy.*Tmp))/sum(Tmp(:))+YIdx(ii);
    XIdxFine(ii)=sum(sum(xx.*Tmp))/sum(Tmp(:))+XIdx(ii);
end

YShift=YIdx-MaxShift-1;
XShift=XIdx-MaxShift-1;

YShiftFine=YIdxFine-MaxShift-1;
XShiftFine=XIdxFine-MaxShift-1;


figure(100);plot([YShift(:) XShift(:) YShiftFine(:) XShiftFine(:)]);
figure(101);plot(XCorrMax);

Shift.YShift=YShiftFine;
Shift.XShift=XShiftFine;
end


%%
function ImgCrt=motion_crt(Img,Shift)
MaxShift=4;

XShiftFine=Shift.XShift;
YShiftFine=Shift.YShift;
TN=length(XShiftFine);

ImgCrt=double(Img(:,:,end-TN+1:end));


[x0 y0]=meshgrid([1:size(Img,2)],[1:size(Img,1)]);

for ii=1:length(XShiftFine)
    ImgCrt(:,:,ii)=interp2(x0,y0,ImgCrt(:,:,ii),x0-XShiftFine(ii),y0-YShiftFine(ii),'linear',0);
end
Mask=ImgCrt(:,:,1)*0;
Mask(MaxShift+1:end-MaxShift,MaxShift+1:end-MaxShift)=1;
ImgCrt=uint16(ImgCrt.*Mask);
end
