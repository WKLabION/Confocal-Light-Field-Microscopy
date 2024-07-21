clear all; 
close all;
addpath("./subfunctions/");
addpath("./example_data");
FilePath=['.\']; %the path of the zip file

FilePrefix=['example_data\ss_single_']; 
FileIdx=[101:115]; %index of the imaging data "ss_single_XXX.tiff"
N=37; % number of microlenses
Bkg=22;
CellPosition=[3 13 23 33 43 53 63 73 83];  %the depth of cells in imaging field of view(FOV), the number could be changed from 1 to the thickness of PSFIn
%% reconstruct the imaging FOV and indexing by cell depth
PSFIn=imstackread([FilePath '20230522_REDpsf_spa2_160UM.tif'],[30:135]);
for kk=1:length(CellPosition)
SaveFileName=['ReconExample_' num2str(CellPosition(kk)) '.tif']; 
Template=padarray(single(PSFIn(:,:,CellPosition(kk))),[100,100])*1000;
Template=(circshift(Template,[1 0])-circshift(Template,[-1,0]))+j*(circshift(Template,[0,1])-circshift(Template,[0,-1]));
figure(1);imagesc(abs(Template));axis image;

%% make a template from the center sub-image
Center=[1275 1433];
ROI=[200 200];
Ref=Template(Center(1)-ROI(1)+1:Center(1)+ROI(1),Center(2)-ROI(2)+1:Center(2)+ROI(2));
figure(2);imagesc(abs(Ref));axis image;

%% calculate the shift distance of each sub-image from corresponding microlens in particular FOV depth
XCorr=fftshift(ifft2(fft2(Template).*conj(fft2(padarray(Ref,[size(Template)-size(Ref)]/2,0,'both')))));
figure(3);imagesc(abs(XCorr));axis image;
Tmp=abs(XCorr);
Peaks=[];
for ii=1:N
    [~,Idx]=max(Tmp(:));
    [yn xn]=ind2sub(size(Tmp),Idx);
    Peaks(ii,:)=[yn xn];
    Tmp(yn-50:yn+50,xn-50:xn+50)=0;
    figure(4);imagesc(Tmp);axis image;
    text(xn,yn,num2str(ii),'color','w');
end

%% Superimpose the sub-images based on the shift distance above 
Recon=zeros(ROI(1)*2,ROI(2)*2,length(FileIdx));
for ii=1:length(FileIdx)
    Data=double(imread([FilePath FilePrefix num2str(FileIdx(ii)) '.tiff']))-Bkg;  
    Img=padarray(Data,[100,100]);
    for jj=1:N
        Recon(:,:,ii)=Recon(:,:,ii)+Img(Peaks(jj,1)-ROI(1)+1:Peaks(jj,1)+ROI(1),Peaks(jj,2)-ROI(2)+1:Peaks(jj,2)+ROI(2));
    end
end
imstackwrite(uint16(Recon),[FilePath SaveFileName]);
end
