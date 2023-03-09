**A simple test to uncover signals of CpG hypermutability using posterior predictive simulations**\
Simon Laurin-Lemay, Nicolas Rodrigue,\

---

## 0. Get Started

```
git clone https://github.com/Simonll/CpG-ppred-test.git
cd CpG-ppred-test
```
## 1. Install Docker

https://docs.docker.com/engine/install/

## 2. Install Snakemake

https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

## 3. Build all Docker containers needed

To run the workflow, you need to build the following Docker containers:

### Primary Docker layer
```bash
docker build --build-arg USER_NAME=$(whoami) --build-arg USER_ID=$(id -u ${USER}) --build-arg GROUP_ID=$(id -g ${USER}) -t ubuntu20.04/basic:latest https://github.com/Simonll/docker.git#develop:/dockerfiles/basic --pull
```

### Specific layer for IQTREE
```bash
docker build -t ubuntu20.04/iqtree:latest https://github.com/Simonll/docker.git#develop:/dockerfiles/iqtree
```

### Specific layer for Phylogenetic Simulator
```bash
docker build --build-arg CACHEBUST=$(date +%s) -t ubuntu20.04/lfp:latest https://github.com/Simonll/docker.git#develop:/dockerfiles/LikelihoodFreePhylogenetics
```

### Specific layer for Phylobayes MPI
```bash
docker build -t ubuntu20.04/pbmpi:latest https://github.com/Simonll/docker.git#develop:/dockerfiles/phylobayes-mpi
```

## 4. Replicate experiments

You need to process all Snakemake files in a semi-automatic way. Begin by activating parts of the workflow selectively (comment out certain sections to obtain the necessary results). Using notebook scripts will help you generate the statistics and figures.
