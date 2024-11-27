# Building Confidence in Deep Generative Protein Design

## Protein Structures Data

```
BLDB:   Beta-lactamases
CYTC:   Cytochrom c
GFP:    GFP
Ras:    Ras GTPase
```

Metadata of protein structures used for training can be found in `/data/*.csv`.

### Generated Structures

```
B:      Backbone only
BQ:     Backbone + Optimal predicted sequences
BQS:    Backbone + Optimal preidcted sequences + homology modeled side-chains
BQSH:   Backbone + Optimal preidcted sequences + homology modeled side-chains + hydrogen atoms added
C:      Conserved residues
Q:      Optimal predicted sequences
```

## Dependencies

### Python Packages

```
TODO
```

### TM-align

Forked from https://zhanggroup.org/TM-align/ [^1].

[^1]: Zhang, Y. TM-align: a protein structure alignment algorithm based on the TM-score. Nucleic Acids Research 33, 2302–2309. issn: 1362-4962. http://dx.doi.org/10.1093/nar/gki524 (Apr. 2005).

```
cd tmalign
g++ -static -O3 -ffast-math -lm -o TMalign TMalign.cpp
chmod +x ./TMalign
```

### GROMACS

https://manual.gromacs.org/current/download.html [^2]

```
cd gromacs
tar xfz gromacs-2024.2.tar.gz
cd gromacs-2024.2
mkdir build
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON
make
make check
sudo make install
source /usr/local/gromacs/bin/GMXRC
```

For cuda users,
```
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_GPU=CUDA -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc
```

`min_sd.mdp`, `nvt_heat.mdp`, `npt_prod.mdp`, `npt_eq.mdp`, `to_origin.tcl`, `del_wat_inside.tcl` are derived from
https://github.com/allison-group/structural-phylogenetics-bootstrap/blob/master/MD/GMX/.

[^2]: Abraham, M. J. et al. GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers. SoftwareX 1–2, 19–25. issn: 2352-7110. http://dx.doi.org/10.1016/j.softx.2015.06.001 (Sept. 2015).

### CHARMM36 Force Field for GROMACS

https://mackerell.umaryland.edu/charmm_ff.shtml#gromacs

charmm36-jul2022.ff.tgz (July 2022, ver.)
https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-jul2022.ff.tgz

```
cd gromacs
tar xfz charmm36-jul2022.ff.tgz
```

### VMD

Download VMD (e.g., `vmd-1.9.4a57.bin.LINUXAMD64-CUDA102-OptiX650-OSPRay185.opengl.tar.gz`) from 
https://www.ks.uiuc.edu/Research/vmd/ to `./vmd`.

```
cd vmd
tar xfz vmd-1.9.4a57.bin.LINUXAMD64-CUDA102-OptiX650-OSPRay185.opengl.tar.gz
cd vmd-1.9.4a57
./configure
cd src
sudo make install
```

### Autodock Vina

https://vina.scripps.edu/

```
conda install conda-forge::vina
```

### ADFR Suite

https://ccsb.scripps.edu/adfr/downloads/

```
conda install hcc::adfr-suite
```

### Reduce

```
conda install bioconda::reduce
```

### ClustalW

```
conda install bioconda::clustalw
```

### Modeller

```
conda install salilab::modeller
```

### Emboss

```
conda install bioconda::emboss
```

conda install bioconda::mgltools
