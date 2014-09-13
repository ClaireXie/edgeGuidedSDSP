function output=quadtree(input)

% tappering
input=double(input);
tmp=zeros(256,256);
tmp(1:size(input,1),1:size(input,2))=input;
input=tmp;

input=uint8(input);
S = qtdecomp(input,.05, 8);
blocks = repmat(0,size(S));

for dim = [128 64 32 16 8];    
  numblocks = length(find(S==dim));    
  if (numblocks > 0)        
    values = repmat(1,[dim dim numblocks]);
    values(2:dim,2:dim,:) = 0;
    blocks = qtsetblk(blocks,S,dim,values);
  end
end

blocks(end,1:end) = 1;
blocks(1:end,end) = 1;

output=(blocks)*255;