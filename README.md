# DAGGER

Code for replicating the experiments in the paper [DAGGER: A sequential algorithm for FDR control on DAGs](https://arxiv.org/pdf/1709.10250.pdf) by [Aaditya Ramdas](http://people.eecs.berkeley.edu/~aramdas/), [Jianbo Chen](http://www.jianbochen.me), [Martin J. Wainwright](https://people.eecs.berkeley.edu/~wainwrig/), [Michael I. Jordan](https://people.eecs.berkeley.edu/~jordan/). 

## Dependencies
The code for DAGGER runs with Python 2.7. Please `pip install` the following packages:
- `numpy`
- `scipy` 
- `graph_tool`
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

![alt-text-1](https://raw.githubusercontent.com/Jianbo-Lab/DAGGER/master/graph_structure/results/figs/mountain-valley-power.eps) ![alt-text-2](https://github.com/Jianbo-Lab/DAGGER/blob/master/graph_structure/results/figs/mountain-valley-power.eps)
<center>This figure shows the power of various algorithms on the entire GO graph (left) and the subgraph rooted at cell proliferation (right), under the setting of independent p-values (top) and Simes' p-values (bottom). </center>

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

# Running experiments in Section 4.3
python run_simulations.py --experiment mountain_vs_valley --n_replications 10
python run_simulations.py --experiment shallow_vs_deep --n_replications 10
python run_simulations.py --experiment diamond_vs_hourglass --n_replications 10
```

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

## Acknowledgements: 
We thank Jelle Goeman for helping running code from his packages and reproducing results in his papers. We also thank Lihua Lei for sharing code for SCR-DAG and BH-DAG.
