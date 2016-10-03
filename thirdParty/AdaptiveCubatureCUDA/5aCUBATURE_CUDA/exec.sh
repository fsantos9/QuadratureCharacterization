#g++ -c CubatureTest.cpp
nvcc -c -arch=sm_20 -o Adap.o Adap.cu
nvcc -c -arch=sm_20 -o Heap.o Heap.cu
nvcc -c -arch=sm_20 -o CubaCuda.o CubaCuda.cu
nvcc -arch=sm_20 -o CubaCuda  *.o
rm *.o
