function output=scaleTo(r,max_val)
if (max(r(:))-min(r(:))==0)
    output=zeros(size(r(:)));
else
    output=(r-min(r(:)))/(max(r(:))-min(r(:)))*max_val;
end
