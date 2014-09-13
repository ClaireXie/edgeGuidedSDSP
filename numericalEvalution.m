clc;clear;close all;

addpath('mainCode/');

% parameters to change
names=cell(4);
names{1}='cones';
names{2}='teddy';
names{3}='tsukuba';
names{4}='venus';
%----------------%
%names{1}='cones';
%scaleFact = 4;
%----------------%


global window sigma_d;
window=7;
sigma_d=0.5;
scale=4;
threshold=0.08;
offset=0;
w2=0;

global fid;
fid=fopen('resultPara_7_10_2.txt','w');

for w1=10
    for localSize=1
        
        for i=1:4
           if (i==1 || i==2)
                scaleFact = 4;
            elseif (i==3)
                scaleFact = 16;
            elseif (i==4)
                scaleFact = 8;
           end
           
           %fprintf(['window= ', num2str(window), '; sigma_d= ', num2str(sigma_d),'\n']);
           fprintf(fid, ['w1= ', num2str(w1), '; localSize= ', num2str(localSize),'\n']);
           
           inputFile=names{i};
           %fprintf(['runnning image ',inputFile,'\n']);
           [highres, edges]=mrfLearning2(names, i, w1, w2, localSize, scale, threshold, 0);

           % run evaluation
           border=2*window;
           runEvaluation(inputFile, scale, highres, edges, scaleFact, border, threshold);
           fprintf(fid, '\n');
        end
        
    end
end



