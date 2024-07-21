function imstackread(imstack,filename)

imwrite(imstack(:,:,1),filename);
for ii=2:size(imstack,3)
%     ii
    imwrite(imstack(:,:,ii),filename,'WriteMode','append');
end

end