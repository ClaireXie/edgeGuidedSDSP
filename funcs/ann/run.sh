#compile
g++ knn.cpp -I. -L. -lANN -o knn

#run
dim=882 #dimension of the data is 882
k=5
./knn data.txt query.txt $dim 5