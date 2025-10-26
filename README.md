This repository provides a Transformer-based Physics-Informed Neural Network (PINN) framework for simulating the Landau damping effect based on kinetic theory.

Kinetic Simulation
Navigate to the kinetic_simulation directory and run: nvla1d.m. This script performs a kinetic simulation of the Landau damping effect.

Data Processing
After the simulation, run: t_test.m. This script processes the generated data and converts it into a format compatible with the Python training scripts.

Transformer-Based PINN Training
In the transformer_based_pinn directory, execute: torchrun --nproc_per_node=4 ld_tf_ddp_pq.py

Result Visualization
Use: python analysis.py
