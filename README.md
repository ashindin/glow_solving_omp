# glow_solving_omp
Inverse problem solving for 2 images of artificial airglow (bistatic airglow registration) powered by openmp
## Dependencies:
1. wget.
2. python3 with requests, json and os modules.
3. gfortran or ifort compiler.
## To download test data use:
```bash
./python download_test.py
```
## To build with gfortran use:
```bash
./build.sh
```
## To build with ifort use:
```bash
./build_ifort.sh
```
## To build with ifort for Intel MIC device use:
```bash
./build_ifort_mic.sh
```