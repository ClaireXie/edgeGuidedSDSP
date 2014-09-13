clc;close all;clear;

points = [113, 114, 115, 116, 117, 118, 119, 120, 121; ...
          63,  63,  63,  62,  62,  61,  60,  60,  60];

plot(points(2,:)*4,points(1,:)*4, ':') 
hold on;
plot(points(2,:)*4,points(1,:)*4, 'r*') 
values = spcrv(points,4)*4; 
plot(values(2,:),values(1,:));