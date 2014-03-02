rmdir b
mkdir b
cd b
cmake .. -DCGAL_DIR=k:\CGAL-4.3\build
msbuild cs.sln /p:Configuration=Release
cd ..
