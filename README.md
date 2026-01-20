# Regime-Switching Heston Model with Onsager–Machlup Functional

This repository contains the codebase for studying regime identification in a Heston stochastic volatility model using both a classical Hidden Markov Model (HMM) approach and a path-based approach built on the Onsager–Machlup (O–M) functional.

The project focuses on comparing standard likelihood-based inference with path-wise optimality criteria under both baseline and challenging experimental settings.

---

## Project Structure & Usage

This section describes the organization of the repository and provides guidance on how to reproduce the experimental results.

### 1. Main Execution Scripts

To replicate the results reported in the study, use the following entry-point scripts:

- **`Base_All`**  
  Executes all baseline experiments, serving as a reference to establish standard performance metrics under well-separated regimes.

- **`rHard_All`**  
  Executes the challenging-case experiments, designed to assess model robustness under adverse conditions such as elevated noise levels or rapid regime switching.

---

### 2. Model Variants

The repository includes two main modeling frameworks:

- **`Heston_Standard`**  
  Implements the classical Hidden Markov Model (HMM) approach applied to the Heston stochastic volatility model.

- **`Heston_OM`**  
  Implements the proposed framework based on the Onsager–Machlup functional, enabling path-based inference for regime identification.

---

### 3. Visualization & Results

- **`visualization/`**  
  Contains scripts for processing logged outputs and generating diagnostic and performance plots.

- **`Figure/`**  
  Stores all generated figures, including regime probability maps and estimated volatility trajectories.




