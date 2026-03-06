
# Kinetic Proofreading and Evolutionary Dynamics — MATLAB Code Repository

This repository contains the full MATLAB implementation used to investigate how kinetic proofreading regulates the balance between mutational supply and mutational robustness, and how this interaction shapes genome architecture, population size effects, and adaptive responses to environmental perturbations. All scripts and functions reproduce the analytical results, simulations, and figures described in the accompanying manuscript.

---

## 📁 Repository Contents

### 1. `err_fraction.m`
Computes the minimum error rate \( f_{\min} \) and the corresponding optimal driving-rate constant \( m_0 \) for a given dissociation rate constant \( K_x \).

### 2. `error_calculator.m`
Evaluates the replication error rate for any specified driving-rate constant \( m \) and dissociation rate constant \( K_x \).

### 3. `driving_function.m`
Returns the driving-rate constant \( m \) required to achieve a given error rate \( f \), for a specified dissociation rate constant \( K_x \).

### 4. `K_vs_Temp_shift.m`
Simulates how temperature shifts \( \Delta T \) modify the dissociation rate constant, generating the relationship:
\[
\Delta T \rightarrow K_{x,\text{new}} / K_{x,\text{old}}.
\]

### 5. `del_f_vs_temp_shift.m`
Computes the effect of temperature shifts on the replication error rate, producing:
\[
\Delta T \rightarrow f_i / f_{\min}.
\]

### 6. `error_plots_diff_K_vals.m`
Generates error–driving-rate curves for multiple dissociation rate constants \( K_x \), illustrating how proofreading kinetics depend on underlying molecular parameters.

### 7. `KPR_enzyme_evolution_indi_errors.m`
Simulates evolutionary dynamics under temperature perturbation where **each individual experiences a distinct phenotypic shift** in driving-rate constant. Tracks population size, mean trait evolution, and adaptation rate.

### 8. `KPR_based_enzyme_evolution.m`
A population-level simulation assuming **uniform phenotypic change** across all individuals following environmental perturbation.

### 9. `KPR_enzyme_evolution_diff_population_size.m`
Explores how **different initial population sizes** \( N_0 \) influence adaptive potential, collapse probability, and post-perturbation recovery.

### 10. `KPR_enzyme_evolution_diff_coding_region.m`
Simulates evolutionary dynamics across varying **coding-region lengths** \( l \), examining how genome size affects mutation load, evolvability, and robustness.

### 11. `KPR_enzyme_evolution_diff_selection_strength.m`
Evaluates evolutionary trajectories under different values of the steepness factor \( a \), which controls selection strength in the offspring-probability function.

### 12. `selection_strength.m`
Visualizes how the steepness parameter \( a \) shapes the offspring probability \( P_{\text{off}} \), helping determine biologically meaningful selection regimes.

### 13. `KPR_enzyme_evolution_with_1st_viability_2nd_speed.m`
Simulates the evolutionary dynamics of a particular coding region length with a given population size under biased culling. Here, the organisms are chosen based on viability as the primary priority and the speed at which the correct product forms as the secondary priority.

### 14. `correct_product_formation_rate.m`
Provides the rate of correct product formation with given kinetic parameters.

---

## 🔧 Requirements

- MATLAB R2018a or later (earlier versions likely compatible)  
- Parallel Computing Toolbox (required for scripts using `parfor`)

---

## 📘 Citation

If you use this repository or adapt any portion of the code, please cite the associated manuscript:

*“Kinetic proofreading maintains the balance between mutational supply and mutational robustness under environmental perturbations”*  
(Authors, Year).

---

## 📬 Contact

For questions, issues, or collaboration inquiries, please open a GitHub Issue or contact the repository maintainer.

