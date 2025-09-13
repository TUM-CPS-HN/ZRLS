## Online Data-Driven Reachability Analysis using Zonotopic Recursive Least Squares
<br/> 
This repository contains the code for the following paper:

<br/><br/>
**Amr Alanwar, Anne Koch, Frank Allgöwer, Karl H. Johansson**  
*"Online Data-Driven Reachability Analysis using Zonotopic Recursive Least Squares"*

## A Short Video About the Idea
(Add link here if available)

## Problem Statement
We present a data-driven reachability analysis framework that computes **over-approximations of reachable sets directly from online state measurements**.  

The method estimates **time-varying unknown models** using an **Exponentially Forgetting Zonotopic Recursive Least Squares (EF-ZRLS)** algorithm, which processes data corrupted by bounded noise.  

Specifically:  
- A time-varying set of models that contains the true system model is estimated recursively.  
- This estimated set is then used to compute the forward reachable sets under process noise and uncertain inputs.  

<br />
The following figure summarizes the idea behind our paper:
<br /> <br />

<p align="center">
<img
src="Figures/diagram.pdf"
raw=true
alt="Framework Diagram"
width=500
/>
</p>
<br />
<br />

## Files Description 
This repository contains the simulation code used in the paper.<br />

## Running 
1- Download [CORA 2025](https://tumcps.github.io/CORA//pages/archive/v2025/index.html).<br />
2- Add CORA and its subfolders to the MATLAB path.  <br />
3- Add this repository and its subfolders to the MATLAB path.  <br />
<br />
<br />

## Basic Reachability (folder: `Examples`)
1- Run **Example1_A1.m** and **Example1_A2.m** for Example 1.1 in the paper (Reachability Analysis of LTI Systems).<br />
2- Run **Example1_B.m** for Example 1.2 in the paper (Reachability Analysis of LTV Systems).<br />
3- Run **Example2.m** for Example 2 in the paper (Reachability Analysis of Lipschitz Nonlinear Systems).<br />
<br />
<br />
<br />
<br />
<br />

## Citation
If you use this repository, please cite our paper:<br />

```
@article{naderi2025online,
  title     = {Online Data-Driven Reachability Analysis using Zonotopic Recursive Least Squares},
  author    = {Naderi Akhormeh, Alireza and Hegazy, Amr and Alanwar, Amr},
  journal   = {To be filled},
  year      = {2025},
  volume    = {},
  number    = {},
  pages     = {},
  publisher = {}
}

```
