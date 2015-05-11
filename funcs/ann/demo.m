% test knn_mex.cpp
clc; clear;

data = load('data.txt');
query = load('query.txt');
[index, dist] = knn_mex(data, query, 5, 1);
