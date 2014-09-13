% hightlight the a pardicular patch centered @ center with a fixed size

function highlight(edge, center, size)

output(:,:,1)=edge;
output(:,:,2)=edge;
output(:,:,3)=edge;

%drawlines
output(center(1)-size, center(2)-size:center(2)+size, 1)=1;
output(center(1)+size, center(2)-size:center(2)+size, 1)=1;
output(center(1)-size:center(1)+size, center(2)-size, 1)=1;
output(center(1)-size:center(1)+size, center(2)+size, 1)=1;

figure;imshow(uint8(output*255));

