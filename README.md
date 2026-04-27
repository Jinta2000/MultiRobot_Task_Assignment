# Distributed Multirobot Task Assignment via Consensus ADMM

## 📌 Overview

This repository contains the simulation and implementation of the paper:

**"Distributed Multirobot Task Assignment via Consensus ADMM"**

The project focuses on solving multirobot task allocation problems using distributed optimization techniques. Specifically, it implements the Consensus ADMM (Alternating Direction Method of Multipliers) framework for coordinating multiple robots in persistent surveillance scenarios.

---

## 🎯 Objectives

* Implement distributed task allocation for multiple robots
* Ensure continuous surveillance of designated nodes
* Integrate battery constraints and charging station management
* Compare performance of different algorithms:

  * MUR-TAP
  * MURD-TAP
  * MURID-TAP

---

## 🧠 Problem Description

Given:

* **N robots**
* **Surveillance stations (m_s)**
* **Charging stations (m_c)**

The goal is to:

* Assign robots to surveillance tasks efficiently
* Ensure no surveillance node is left unattended
* Manage battery levels with charging constraints
* Maintain distributed decision-making without central control

---

## ⚙️ Algorithms Implemented

### 1. MUR-TAP

* Multi-robot task assignment problem using ADMM
* Ensures consensus among robots

### 2. MURD-TAP

* Distributed extension of MUR-TAP
* Each robot performs local updates and communicates with neighbors

### 3. MURID-TAP

* Improved distributed algorithm *

---

## 🛠️ Implementation Details

* Language: **MATLAB**
* Optimization Method: **Consensus ADMM**

### Key Features

* Dynamic robot-task assignment
* Battery level tracking
* Charging station allocation
* Constraint handling (one robot per surveillance node)

---

## ⚠️ Constraints Considered

* Only one robot per surveillance node
* Charging station occupancy constraints
* Battery threshold conditions
* Continuous surveillance requirement

---

## 🔍 Results Summary

* ADMM-based approaches converge efficiently
* Distributed methods reduce centralized computation
* MURID-TAP improves robustness under battery constraints

---

## 🚀 Future Work

* Incorporate obstacle avoidance
* Real-world robot deployment
* Extend to heterogeneous robot teams
* Improve energy optimization strategies

---

## 📚 Reference

Paper: *Distributed Multirobot Task Assignment via Consensus ADMM, Ola Shorinwa , Ravi N. Haksar , Patrick Washington, IEEE TRANSACTIONS ON ROBOTICS, VOL. 39, NO. 3, JUNE 2023*

## 👤 Author: 

Jinta Gladson

## ⭐ Notes

* This project is intended for academic and research purposes
* Code may require tuning for different simulation setups
