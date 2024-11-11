# MUSIC, MUSUN

> Several ==Monte Carlo codes== are able to transport muons through matter with high accuracy. The codes can be split in two categories: 
>
> (i) ==multipurpose particle transport codes==, such as GEANT4 and FLUKA;
>
> (ii) codes developed ==specifically for muon propagation== through **large thickness of material**, such as PROPMU, MUSIC, MUM and MMC.

#### MUSIC和MUSUN的关系

MUSIC (MUon SImulation Code) and MUSUN (MUon Simulations UNderground). 

MUSIC is a **package for muon transport through matter**. It is particularly useful for propagating muons through <u>large thickness of rock or water</u>, for instance from the surface down to underground/underwater laboratory. 

MUSUN is designed to <u>use the results of muon transport through rock/water</u> to **generate muons** in or around underground laboratory taking into account their energy spectrum and angular distribution.



#### PROPMU, MUM and MMC 与 MUSIC 的对比

> MUM is one-dimensional and does not take into account muon deflection in the plane perpendicular to the initial direction. 
>
> MMC considers only muon deviation due to multiple scattering and not many details or results are given in the original paper
>
> PROPMU also treats only multiple scattering
>
> ​	Muon transport codes MUSIC,MUM and MMC were found to agree with each other giving similar muon energy distributions beyond large thickness of rock or water and similar muon survival probabilities. MUSIC and PROPMU are in agreement for muon transport in rock but results obtained with PROPMU in water were found to be different from those obtained with MUSIC and MUM.



### MUSIC

FORTRAN

There are two major versions of MUSIC existing and developed in parallel: 

(i) ‘standard’, dedicated for muon transport through large thickness of matter; 

(ii) ‘thin slab’, developed for muon transport through thin slabs of materials.

> The code takes into account the energy losses of muons due to **four processes**: 
>
> - ionisation (using Bethe–Bloch formula) including knock-on electron production, 
>
>   <img src="/home/liuruiguoguo/.config/Typora/typora-user-images/image-20200406165111615.png" alt="image-20200406165111615" style="zoom:25%;" />
>
> - bremsstrahlung (or braking radiation),
>
> -  electron–positron pair production  
>
> - muon-nucleus in elastic scattering (or photonuclear interactions).

#### standard
> - The standard version of MUSIC considers all interaction processes stochastically if the fraction of energy lost by a muon in the interaction exceeds a predefined value of a parameter, $v_{cut}$. (大概讲了：$v>v_{cut}$时，由随机数从CERN library中抽样得到每次反应的自由程，随机数选择反应类型，随机数控制能量衰减，$v<v_{cut}$时，Bethe–Bloch公式直接计算能损)
> - Molière theory does not provide the lateral displacement that is sometimes more important from the experimental point of view (for instance,for muon bundles underground)
>  - 除了库伦散射，还考虑了其他相互作用造成的散射

#### thin slab

​	This version has been aimed at **providing muon energy, position and direction** at the end of every small segment of muon path in water and at passing the muon energy loss at the segment to another part of software that generated **Cherenkov photons**. Since then this version has also been used to estimate muon deflection in <u>high-A materials</u> (iron, lead and uranium) for possible <u>security applications</u> (searching for hidden high-A materials, like uranium, in cargo).

###### comparing to ‘standard’ version

> The <u>‘standard’ version</u> of MUSIC, in the absence of stochastic interactions on the <u>small segment</u>, always returns the <u>mean value</u> for continuous energy loss at the end of the segment. In the <u>‘thin slab’ version</u> the cut that separates stochastic and continuous parts of the energy loss is reduced to 1 MeV, meaning that practically <u>all muon interactions are stochastic</u>. The ionisation energy loss is calculated using Landau distribution (call to a function from the CERN library).

​	两种版本的MUSIC（标准平板和薄平板）对于物质厚板都能提供一致的结果，但是“薄平板”版本的CPU消耗更多，因为$v_{cut}$的值较低。 薄板版本可以提供比材料薄板更精确的能谱和角度偏差结果。

#### MUSIC应用

MUSIC has been used in the analysis of SNO and MACRO data. The code has also been applied for the calculation of expected background induced by <u>cosmic-ray muons in deep underground experiments</u>, such as KamLAND, Super-Kamiokande, etc.

#### MUSIC与Geant4、FLUKA模拟结果对比



### MUSUN

FORTRAN

##### 应用：

- In a <u>standard version</u> of MUSUN the vertical depth should be <u>more than 500 m w. e.</u> There is no strict upper limit for the vertical depth, but the maximum slant depth should not exceed 15 km w. e.  (w.e.?)

> In a simple version of MUSUN, the flat profile is assumed for the surface above the underground site (the curvature of the Earth is taken into account but other possible fluctuations of the slant depth are ignored).
>
> MUSUN offers the choice of the muon energy spectrum (as described above), the fraction of prompt muons, the vertical depth of the laboratory, the range of zenith and azimuthal angles, and the range of energies. <u>No additional muon propagation is required for different options. Different types of rocks (rock compositions),however, require separate muon transport.</u>（？）

- For practical purposes (for instance, when these muons are used in multipurpose event generators GEANT4 or FLUKA) it is useful to generate muons on the surface of a rectangular parallelepiped or a sphere with predefined dimensions. MUSUN offers a possibility to generates muon positions on the surface of a rectangular parallelepiped with dimensions specified by the user. Muon parameters are written on the disk and can be passed later on to the multipurpose event generators.

- It has been used to study muon-induced neutron background for experiments looking for rare events, such as WIMPs

### Conclusions

The two Monte Carlo codes MUSIC and MUSUN dedicated to muon simulations have been described. MUSIC, a package for muon transport through matter, **can be used for propagating muons through large thickness of rock or water**, for instance from the surface down to underground/underwater laboratory. It can also be implemented in the **event generators** for large underwater/under-ice neutrino telescopes or other neutrino detectors. MUSUN uses the results of muon transport through rock/water to generate muons in or around underground laboratory taking into account their energy spectrum and angular distribution. <u>*Various tests showed good agreement of the codes’ results with experimental data and other packages*</u>.