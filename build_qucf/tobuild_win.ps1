$CPPFiles = "../../src/QCircuit.cpp;"
$CPPFiles += "../../src/QLib.cpp;"
$CPPFiles += "../../src/QuEST_lib.cpp;"
$CPPFiles += "../../src/QGates.cpp;"
$CPPFiles += "../../src/BaseTool.cpp;"
$CPPFiles += "../../src/QuCF.cpp;"
$CPPFiles += "../../src/launch_QuCF.cpp;"

cmake `
    -S "../submodules/QuEST" `
    -B build `
    -DCMAKE_TOOLCHAIN_FILE="..\..\..\..\vcpkg/scripts/buildsystems/vcpkg.cmake" `
    -DOUTPUT_EXE="QuCF" `
    -DUSER_SOURCE="$CPPFiles" `
    -DGPUACCELERATED=1 `
    -DGPU_COMPUTE_CAPABILITY=89 `
    -DCUDA_LIBRARIES=CUDA::cudart