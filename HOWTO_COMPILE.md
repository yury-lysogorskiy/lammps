# Compilation with `cmake`

`mkdir build; cd build`

`cmake -DCMAKE_BUILD_TYPE=Release -D BUILD_MPI=ON -DPKG_USER-PACE=ON  -DPKG_USER-PACE-AL=ON  -DPKG_MANYBODY=ON ../cmake`

`cmake --build .`

NOTE: If there are compilation error messages regarding C++11 standards, run 
`cmake ../cmake -DPKG_USER-PACE=ON -DPKG_USER-PACE-AL=ON -D CMAKE_CXX_FLAGS="-std=c++11"`
