# Usage

## 安装和编译
```bash
git clone --recursive git@github.com:FR13ndSDP/OpenCFD-AMR.git
cd OpenCFD-AMR

cmake -B build; cd build
make -j <N> # thread number
```

## 运行
编译选项见 `CMakeLists.txt`

```CMake
set(EBAMR_CASE            "7_JET"  CACHE STRING "Case folder")

#
# Physics options
#
option(OPTION_USE_CHEM      "Enable chemical reaction"     ON)

#
# HPC options
#
option(OPTION_MPI    "Enable MPI"    ON)
option(OPTION_OPENMP "Enable OpenMP" OFF)
option(OPTION_CUDA   "Enable CUDA"   OFF)
option(OPTION_HIP    "Enable HIP"    OFF)
option(OPTION_SYCL   "Enable SyCL"   OFF)
```
以此设置为例

```bash
cp ../Exec/7_JET/inputs .
mpirun -n 8 ./EBR.exe inputs
```