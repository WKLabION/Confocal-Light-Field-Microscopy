function imstack=imstackread(varargin)
if nargin==2
    filename=varargin{1};num=varargin{2};
    info=imfinfo(filename);
    for ii=1:length(num)
        imstack(:,:,ii)=imread(filename,num(ii));
    end
    
else
    filename=varargin{1};
    info=imfinfo(filename);
    for ii=1:size(info)
        imstack(:,:,ii)=imread(filename,'Info',info(ii));
    end
end
% disp(['load: ' filename ': finished']);
end