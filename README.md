# Fault-Tolerant Control on Subsea Injection System

This repository contains simulations for faulted cascade pump on the subsea seawater injeciton system using proposed hierarchical resilient (fault-tolerant) control (HRC) framework (supervisory RTO, stochastic offset-free NMPC, and nonlinear observer), supporting our submission to the SICE Festival 2026 with Annual Conference.

## Purpose

Sharing MATLAB code and data for reproducing results and facilitating review.

## Repository Structure

* **`PumpHRC.m`**: MATLAB simulation code for HRC framework of the cascade pump system.
* **`FeedDist.mat`**: Sequential data to accompany `PumpHRC.m`.

## Requirement

* [Matlab Optimization Toolbox](https://se.mathworks.com/products/optimization.html)
* [YALMIP: A Toolbox for Modeling and Optimization in MATLAB](https://yalmip.github.io)
* [Gurobi Optimization](https://www.gurobi.com)

## SICE Festival 2026 with Annual Conference

This repository directly supports data and code for our paper submission.
