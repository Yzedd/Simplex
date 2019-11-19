# Simplex

A python implementation coupling revised simplex method, revised dual simplex method and two phase method to solve linear programming problem is introduced. This python script provides an interactive manipulation for the convenience of users. Rather than solving the optimal solution, it is also able to identify case for infeasible solution and unbounded solution. Moreover, the cycling of the simplex method caused by degeneracy can be avoided. 

This python implementation is built under Python 3.6.8.
To use this python script to solve linear optimization problems, required steps are as follow with an example linear optimization problem:

Example:

min 𝑧(𝑥) = 3𝑥1 + 2𝑥2

𝑠𝑢𝑏𝑗𝑒𝑐𝑡 𝑡𝑜:

3𝑥1 + 𝑥2 ≥ 3

4𝑥1 + 3𝑥2 ≥ 6

𝑥1 + 2𝑥2 ≤ 3

𝑥1, 𝑥2 ≥ 0



Users' Manual

Step 1. Open the python script Simplex.py with Anaconda Spyder or other cross-platform integrated development environment (IDE) in the Python language. Click run file or press F5 to run the python script.

Step 2. After running the python script, a statement “Make sure it’s a minimization problem =>” is shown in the IPython console. Press Enter.

Step 3. “Objective coefficients =>” is shown in the IPython console. Enter the objective coefficients, only single spaces should be used. e.g. “3”, “space”, “2”.

Step 4. “Constraints =>” is shown in the IPython console. Enter the constraints, e.g. “-3”, “space”, “-6”, “space”, “3”.

Step 5. “Signs (only ‘<’ or ‘=’) =>” is shown in the IPython console. Enter the signs of constraints, e.g. “<”, “space”, “<”, “space”, “<”.

Step 6. “Constraint coefficients =>” is shown in the IPython console. Enter the constraint coefficients, e.g. “-3”, “space”, “-1”, “space”, “-4”, “space”, “-3”, “space”, “1”, “space”, “2”.

Step 7. The final optimization result is shown in the IPython console, e.g. “Iteration 2”, “Optimal solution: x = [6.0000000e-01 1.2000000e+00 8.8817842e-16 0.0000000e+00 0.0000000e+00], z = 4.200000000000001”. The example result shows that after two iterations, the targeted linear programming problem reaches optimal solution, with 𝑥1 = 0.6, 𝑥2 = 1.2, 𝑚𝑖𝑛 𝑧(𝑥) = 4.2.
