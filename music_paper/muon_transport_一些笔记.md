## muon_transport

| muon_transport             |                              |
| -------------------------- | ---------------------------- |
| **x0,y0,z0**               | 缪子初始位置                 |
| **cx0,cy0,cz0**            | 缪子初始方向(cos值)          |
| **Emuin0**,depth0,**tmu0** | 初始能量,传播深度,初始化时间 |
| idim,idim1                 | 3D                           |

```fortran
!初始化定义
        theta=0.                !angle with respect to z-axis, input and output
        phi=0.                  !angle with respect to x-axis, input and output
        dr=0.                   !muon path in the rock, output
        x=0.	y=0.	z=0.  
        t=0.                    !Initial muon path in rock, input
        emu_f=0.                !final muon energy
!根据input赋值
        path_max=1000000.*100.  !maximal path of muon, input
        depth=depth0*n_rho      !depth along z-axis in g/cm^2, input
        emu=emuin0              !initial muon energy, input
!即将带入计算        
        emu_0=emu
        theta0=theta
        phi0=phi
```

CALL music (emu,depth,emu_f,x,y,z,t,theta,phi,dr,path_max,idim,idim1)

```fortran
	   gamma=(Emu_0+Emu_f)/2./0.105658
        beta=dsqrt(gamma*gamma-1.)/gamma
        vmu=beta*29.9792458
        tmu=dr/vmu/n_rho
        tmu0=tmu0+tmu
        
        thetac=dacos(cz0)
        phic=datan(cy0/cx0)		!方向角
```

call coord_transform(x1,y1,z1,x,y,z,thetac,phic)

call angle_transform(theta0,phi0,theta,phi)

```
	   emuin0=emu_f                   !if emu_f is less than 0.106 the muon was stopped
        if(emuin0.le.0.106d0) emuin0=0.d0
```



### music

Emu,zf,Emu_f,x,y,z,t0,theta,phi,dr,t10,ms_flag,d_flag

| music                            |                                        |
| -------------------------------- | -------------------------------------- |
| emu, **emu_f**                   | 入射能量,出射能量                      |
| **x y z, theta phi**, t0, **dr** | 坐标, 方向角, 初始路程, 又走了多少路程 |
| zf, t10                          | 在哪出射(岩石厚度), 一个缪子最多走多远 |
| ms_flag, d_flag                  | 考虑是否进行散射的参数                 |

> IF(em1.LE.mmu) go to 16, else go to 9:
>
> - 16:	IF(em1.LE.mmu) emu_f=mmu -->  结束
> - 9:      next step 开始
>
> if(z.gt.zf) go to 11:
>
> - 准备结束
>
> IF(em1.LE.1.) GO TO 7:
>
> - IF THE MUON ENERGY IS LESS THAN 1 GEV ONLY THE IONISATION LOSSES ARE CALCULATED
>
> IF(z+z1.LE.zf) GO TO 10:
>
> - Calculation of the continuous energy losses of muon up to the next interaction point

CS(100,71,4): 能损(e/emu)划为 j0 份(j0 < 100); log of energy 划为71份; 四种相互作用.

J: 能量对应的序号

b0, bs: continuous ionisation loss, relative energy loss due to process2,3,4 FOR V<10**(-3)



计算连续能损一套流程:

> ​	  CALL vdem(em1,B0,BS)
> ​      emf=em1-(em1*BS+B0)*T
> ​      IF(emf.LE.mmu) emf=mmu
> ​      emi=10.**((dlog10(emf)+dlog10(em1))/2.)
> ​      CALL vdem(emi,B0,BS)
> ​      CALL vem(em1,B0,BS,T,EM2)



CALL multiple(t2,em1,em3,b0,bs,theta2,phi2,deltax,deltay,deltar,err)