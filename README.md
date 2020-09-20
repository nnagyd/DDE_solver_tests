### Mintakódok a TDK dolgozathoz

# Csomagok és útmutatások
### C++ Hand-Tuned_CPU
* Linux és Windows - GNU Compiler: ```g++ -O3 -std=c++17 -mavx -o output -Ipath/to/vcldir input.cpp``` és elképzelhető, hogy a GNU compilerben az alligned_alloc máshogy van definiálva
* Windows - Visual Studio: VCL könyvtár elérési útvonal megadása és ```/Ot /GL /O2 /arch:AVX /std:c++17```

### C++ Hand-Tuned_GPU
* Linux és Windows - NVidia Cuda Compiler: ```nvcc -O3 --std=c++14 --ptxas-options=-v --gpu-architecture=sm_61 -lineinfo -maxrregcount=128 -w --resource-usage main.cu -Ipath/to/Hand-Tuned_GPU_RK4/common```
* gpu-architecture: elérhető CUDA Compute Capability megadása
* Több információ a CUDA-ról [itt](https://developer.nvidia.com/cuda-downloads)

### Julia DifferentialEquation.jl
* Linux és Windows: ```julia -O3 input.jl```
* Vagy IDE-n keresztül, pl: Juno
* Útmutató a könyvtárhoz [itt](https://diffeq.sciml.ai/stable/)

### C++ ParallelDDE
* Linux és Windows - GNU Compiler: ```g++ -O3 -std=c++17 -mavx -o output -Ipath/to/vcldir -I/path/to/parallelDDE input.cpp```
* Windows Visual Studio: VCL és ParellelDDE könyvtár elérési útvonal megadása és ```/Ot /GL /O2 /arch:AVX /std:c++17```
* A ParallelDDE könyvtár letölthető [innen](https://github.com/nnagyd/Parallel-DDE-Solver)

### C++ Retard
* Linux és Windows - GNU Compiler: ```gcc -O3 -Ipath/to/RETARD path/to/retard.c input.c```
* Windows Visual Studio: RETARD könyvtár elérési útvonal megadása
* A RETARD könyvtár letölthető [innen](http://www.unige.ch/~hairer/software.html)

### Wolfram NDSolve
* Linux és Windows - Wolfram Engine: ```wolframscript -f input.wls```
* Windows: Wolfram Mathematicával megnyitható és futtatható
* A kódok az ingyenes [Wolfram Engine](https://www.wolfram.com/engine/)-al készültek
