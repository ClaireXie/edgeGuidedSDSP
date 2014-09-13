clear;clc;

names=cell(4);
names{1}='cones';
names{2}='teddy';
names{3}='tsukuba';
names{4}='venus';
scale=3;

for i=1:4
    edgeGuided=double(imread(['../outputs/v2/s',num2str(scale), '/', names{i},'2_', num2str(scale),'.png']));
    save(['outputs/v2/',names{i},'_eg.mat'],'edgeGuided');
end


