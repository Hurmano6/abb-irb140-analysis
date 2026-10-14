# abb-irb140-analysis
Kinematic and dynamic modeling, and 3D trajectory planning for the ABB-IRB 140 industrial robot using MATLAB.

# Kinematic and Dynamic Analysis of ABB-IRB 140 🤖

![MATLAB](https://img.shields.io/badge/MATLAB-Simulation-blue?logo=mathworks)
![Robotics Toolbox](https://img.shields.io/badge/Toolbox-Robotics_by_Peter_Corke-orange)
![Status](https://img.shields.io/badge/Status-Completed-success)

## 📌 Overview
This repository contains the complete kinematic and dynamic modeling of the **ABB-IRB 140**, a 6-axis industrial robot. The project was developed in MATLAB using the Robotics Toolbox and includes comprehensive simulations of a pick-and-place payload operation, analyzing the robot's physical constraints, singularities, and joint torques.

## 🚀 Key Features

* **Kinematics:** 
  * Forward and Inverse Kinematics derived via Denavit-Hartenberg (DH) parameters.
  * Validation of transformation matrices and mathematical models against the Robotics Toolbox `fkine` and `ikunc` functions.
* **Differential Kinematics & Singularities:** 
  * Jacobian matrix evaluation using `jacob0`.
  * Identification of singularity points (manipulability analysis via $\vert{}det(J)\vert{}$).
* **Dynamic Modeling:**
  * Implementation of the robot's dynamic model ($\tau = M(q)\ddot{q} + C(q,\dot{q})\dot{q} + G(q)$).
  * Analysis of gravitational forces, Coriolis components, and inertial forces across all joints using `rne`, `gravload`, `itorque`, and `coriolis`.
  * Evaluation of an external 2.0 kg payload effect on the joint actuators.
* **Trajectory Planning:** 
  * Simulation of a full *pick & place* task.
  * Integration of joint trajectories (5th-order polynomials via `jtraj`) and Cartesian trajectories (trapezoidal velocity profiles via `ctraj`).

## 🛠️ Technologies & Tools
* **Language:** MATLAB
* **Library:** Robotics Toolbox for MATLAB (by Peter Corke)
* **Concepts:** Denavit-Hartenberg representation, Trajectory Generation, Jacobian Analysis, Rigid Body Dynamics.

## 📁 Repository Structure
* `scripts/`: Contains the MATLAB `.m` files for kinematics, inverse kinematics, Jacobian analysis, and dynamics.
* `report/`: Comprehensive academic paper detailing theoretical foundations and plotted results (Position, Velocity, Torques).

## 👨‍💻 Author
* Sergio Gil Martínez

> *This project was developed as part of the Robotics Engineering degree at the University of Alicante.*
