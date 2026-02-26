rem cl NCEPbin.cpp /O2 /GL /Oi
g++ NCEPbin.cpp -o NCEPbin.exe -std=c++11
erase *.obj
NCEPbin.exe
pause
