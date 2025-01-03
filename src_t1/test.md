```
cd src_t1
gfortran -c ranlux.f
gfortran test-music.f ranlux.o -o test-music
./test-music
rm -rf ranlux.o test-music music-cross-sections.dat music-eloss.dat
```