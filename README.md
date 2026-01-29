<!-- PROJECT SHIELDS -->
[![MIT License][license-shield]][license-url]
[![ReseachGate][researchgate-shield]][researchgate-url]
[![LinkedIn][linkedin-shield]][linkedin-url]
[![GIT][git-shield]][git-url]
<!-- [![Scholar][scholar-shield]][scholar-url] -->
<!-- [![Webpage][webpage-shield]][webpage-url] -->

# Spacecraft Attitude Dynamics and Control
This repository contains the code and models of the project of Attitude Dynamics and Control Systems. 

> M. Tavakoli
> Professors: Hamid Moeenfard, Vahid Bohlouri


1. BFN_00_main.m should be run in first place. It initializes all variables needed. 
2. BFN_00_zfull_system.slx contains full model 
    a. BFN_03_Disturbs.m plots disturbance torques of (2) 
    b. BFN_04_Detumbl.m plots detumbling results using previously acquire data in detumbl.mat 
    c. BFN_05_total_impulse.m calculates total impulse given by thrusters of (2) 
    d. BFN_06_slew.m plots slew maneuver results of (2) 
    e. BFN_07_stabilization.m plots stabilization results of (2) 
3. BFN_01_Gyro_filter.slx contains optic fiber gyro model 
    a. BFN_01_Gyro_filter_results.m plots the results of (3) 
4. BFN_02_Thrusters.slx contains the model of thrusters and calculates the firing sequence. 
    a. BFN_02_Thrusters1.m plots results of (4) 

## Contact
🧑‍💻 Meisam Tavakoli \
📧 [meisam.tavakoli@studio.unibo.it](mailto:meisam.tavakoli@studio.unibo.it)

[git-shield]: https://img.shields.io/badge/Github-fjakob-white?logo=github
[license-shield]: https://img.shields.io/badge/License-MIT-T?style=flat&color=blue

