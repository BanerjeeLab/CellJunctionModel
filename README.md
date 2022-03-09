# Cell Junction Model
Code for theoretical modeling of epithelial junction dynamics as used in "Force-dependent intercellular adhesion strengthening underlies asymmetric adherens junction contraction", Current Biology (2022)


## Code structure

The simulations are run using the "solve" function to numerically integrate the equations of the junction. The script contains examples showing

* Assymetric tension in the model with no feedback between tension and friction and stiffness
* Asymmetric tension in a model with feedback
* A half junction activation

### "solve" function Arguments

| Argument  | Description | Data Type|
| ------------- |:-------------:|:-------------:|
| t          | Timepoints for simulation | array |
| nx         | Number of junction subdivisions | int |
| tension    | Tension along the junction | float or array |
| mu         | Friction along the junction | float or array |
| E          | Junction Young's modulus | float |
| k_top      | Top (x=1) shoulder junction stiffness | float |
| k_bot      | Bottom (x=0) shoulder junction stiffness | float |
