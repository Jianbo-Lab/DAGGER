# DAGGER

Code for replicating the experiments in the paper [DAGGER: A sequential algorithm for FDR control on DAGs](https://arxiv.org/pdf/1709.10250.pdf) by [Aaditya Ramdas](http://people.eecs.berkeley.edu/~aramdas/), [Jianbo Chen](http://www.jianbochen.me), [Martin J. Wainwright](https://people.eecs.berkeley.edu/~wainwrig/), [Michael I. Jordan](https://people.eecs.berkeley.edu/~jordan/). 

## Dependencies
The code for DAGGER runs with Python 2.7. Please `pip install` the following packages:
- `numpy`
- `scipy`  
- `statsmodels`
- `matplotlib`

Or you may run the following and in shell to install the required packages:
```shell
git clone https://github.com/Jianbo-Lab/DAGGER
cd DAGGER
sudo pip install -r requirements.txt
```

To run the experiments in Section 4.1 of the paper, R is required. Please install the following packages in R. 
- `igraph`
- `cherry` 
- `globaltest` in Bioconductor

## Running experiments in Section 4.1
We provide codes for reproducing figures in Section 4.1. Run the following commands in shell:

```shell
###############################################
# Omit if already git cloned.
git clone https://github.com/Jianbo-Lab/DAGGER
cd DAGGER
############################################### 
cd GO
# Generate simulated data (hypothesis and p-values on the GO graph.)
Rscript make_data.R

# Generate all the jobs. 
python generate_jobs.py

# Run an example job.
# set_name: full or cellprolif. which set of nodes on the GO graph is used.
# seed: random seed.
# pvalues_type: simes or ind: whether to use Simes' p-values or independent p-values.
# algorithm: DAGGER,lord,BH,all-goeman,any-goeman,SCR.DAG,BH.DAG,Holm 
# pi0: proportion of nulls. 
python run_simulations.py --set_name full --seed 1 --pvalues_type simes --algorithm DAGGER --pi0 0.5 
# Run all jobs.
sh jobs.sh

# Plot figures.
python plot.py
```

The generated plots can be found in DAGGER/graph_structure/results. See `core/algorithm.py` and `GO/run_simulations.py` for details. 


|*<center>Power comparison for full and subgraph GOs with independent and Simes p-values. </center>*|
|:--:| 
|![alt-text-1](https://github.com/Jianbo-Lab/DAGGER/blob/master/figures/power.png?raw=true "Power comparison for full and subgraph GOs with independent and Simes p-values.")| 

<!-- |*<center>FDR comparison for full and subgraph GOs with independent and Simes p-values. </center>*|
|:--:| 
|![alt-text-2](https://github.com/Jianbo-Lab/DAGGER/blob/master/figures/fdr.png?raw=true)|

|*<center>Time comparison for full and subgraph GOs with independent and Simes p-values. </center>*|
|:--:| 
|![alt-text-3](https://github.com/Jianbo-Lab/DAGGER/blob/master/figures/time.png?raw=true)| --> 

## Running experiments in Section 4.2 and Section 4.3
We provide codes for reproducing figures in Section 4.2 and Section 4.3. Run the following commands in shell:

```shell
###############################################
# Omit if already git cloned.
git clone https://github.com/Jianbo-Lab/DAGGER
cd DAGGER
############################################### 
cd graph_structure
# Running experiments in Section 4.2.
python run_simulations.py --experiment BH_vs_DAGGER --n_replications 10
```

|*<center>BH vs. DAGGER. </center>*|
|:--:| 
|![alt-text-3](https://github.com/Jianbo-Lab/DAGGER/blob/master/figures/bh.png?raw=true)|

```shell
# Running experiments in Section 4.3
# Mountain vs. Valley. 
python run_simulations.py --experiment mountain_vs_valley --n_replications 10 
```

|*<center>Mountain vs. Valley. </center>*|
|:--:| 
|![alt-text-3](https://github.com/Jianbo-Lab/DAGGER/blob/master/figures/mountain_vs_valley.png?raw=true)|

```shell
# Running experiments in Section 4.3
# Shallow vs. Deep.  
python run_simulations.py --experiment shallow_vs_deep --n_replications 10 
```

|*<center>Shallow vs. Deep. </center>*|
|:--:| 
|![alt-text-4](https://github.com/Jianbo-Lab/DAGGER/blob/master/figures/shallow_vs_deep.png?raw=true)|

```shell
# Running experiments in Section 4.3
# Diamond vs. Hourglass.  
python run_simulations.py --experiment diamond_vs_hourglass --n_replications 10
```

|*<center>Diamond vs. Hourglass. </center>*|
|:--:| 
|![alt-text-5](https://github.com/Jianbo-Lab/DAGGER/blob/master/figures/diamond_vs_hourglass.png?raw=true)|

The generated plots can be found in DAGGER/graph_structure/results. See `core/algorithm.py` and `graph_structure/run_simulations.py` for details.

## Citation
If you use this code for your research, please cite our [paper](https://arxiv.org/pdf/1709.10250.pdf):
```
@article{ramdas2017dagger,
  title={DAGGER: A sequential algorithm for FDR control on DAGs},
  author={Ramdas, Aaditya and Chen, Jianbo and Wainwright, Martin J and Jordan, Michael I},
  journal={arXiv preprint arXiv:1709.10250},
  year={2017}
}
```

## Acknowledgements
We thank Jelle Goeman for helping running code from his packages and reproducing results in his papers. We also thank Lihua Lei for sharing code for SCR-DAG and BH-DAG.
