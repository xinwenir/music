music.f CALL:  
    1. music *(emu,depth,emu_f,x,y,z,t,theta,phi,dr,path_max,idim,idim1) #Y
    2. coord_transform(x1,y1,z1,x,y,z,thetac,phic) #Y
    3. angle_transform(theta0,phi0,theta,phi) #Y
    4. ranlux(yfl,1)     ------------------------------------------Need to add
    5. vdem(em1,B0,BS) #Y
    6. vem(em1,b0,bs,t2,em3) #Y
    7. multiple *(t2,em1,em3,b0,bs,theta2,phi2,deltax,deltay,deltar,err) #Y
    8. vcs(j2,yfl,ip,v1) #Y
    9. mu_scatt(ip,em,v1,theta1,phi1) #Y
    10. defl(E,jp,v,theta,phi) #Y
    11. corset(vv,cc,2)  ------------------------------------------Need to add
    12. corgen(cc,vcx,2) ------------------------------------------Need to add

music-crosssections.f CALL:
    1. SIGBK(E,EM,Z,A,FB) #Y
    2. SIGP(E,EM,Z,A,FP) #Y
    3. SIGNU_BS(E,EM,A,FN) #Y
    4. sige1(e,em,z,a,fe) #Y
    5. SIGE(E,EM,Z,A,FE) #Y

test-music.f CALL:
    1. mucrsec(minv,zz0,a0,fr0) #Y------music-crosssections.f
    2. mulos(minv,zz0,a0,fr0,par_ion) #Y------music-crosssections.f
    3. initialize_music(minv,rho,rad) #Y------music.f
    4. rluxgo(3,iranlux,0,0) ------------------------------------------Need to add
    5. rmarin(iranlux1,0,0)  ------------------------------------------Need to add
    6. ranlux(yfl,1)         ------------------------------------------Need to add
    7. muon_transport(x,y,z,cx,cy,cz,emu,depth,ttime,idim,idim1) #Y------music.f


    test-music.f--------|
                        |-------ranlux.f
                        |-------music-crosssections.f
                        |-------music.f---------|
                                                |--------ranlux.f
                                                |--------corset.f
                                                |--------corgen.f-------|       
                                                                        |-------rnormal.f-------|
                                                                                                |-------ranmar.f

The bash file: music.sh is used to generate executable file for the program.

In the above bash script:
First, the Fortran compiler FC is defined as gfortran. You can change it according to the actual Fortran compiler you have installed and are using, such as ifort, etc.
All the relevant source files are listed, including test-music.f and the files it depends on. All these files are in the music directory.
Some compilation options are set. Here, we simply add the options to display warning messages -Wall and -Wextra, and specify the output executable file name as test-music.
Finally, the compilation command is executed. If the compilation is successful (the return value is 0), a success message will be output; otherwise, a failure message will be output and you will be prompted to check the error messages.
You can save the above content as a .sh file, for example, compile_test_music.sh. Then, enter the directory where the file is located in the terminal, grant executable permissions to the file (chmod +x compile_test_music.sh), and finally execute the script (./compile_test_music.sh) to try to compile and generate the test-music executable file.


                                                                                            
