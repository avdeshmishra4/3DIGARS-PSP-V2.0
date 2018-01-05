Compile: mpiCC *.cpp -o main.exe -lm
Run: mpirun -np 4 ./main.exe > console.txt