Checker
====================
This algorithm is built to check user’s solution compatibility and validity.
This algorithm parses user’s solution to build a tree form and check prefixed constraints.
Checked constraints:
  - Batch, Defects, Solution, Optimization parameters files existence.
  - Solution nodes (apart from branches, wastes and residual) identical to batch items.
  - Solution nodes sequencing compared to batch fixed sequence.
  - Solution nodes (apart from branches, wastes and residual) duplication and/or miss.
  - Solution nodes dimensions, compared to its parents and other imposed parameters.
  - Solution nodes (apart from branches, wastes and residual) superposing with defects.

The Checker visualization is a tool to visualize a user generated solution on virtual plates.


Environment:
--------------
This algorithm is coded in C++ language.
The checker is a CodeBlocks project, to recompile it you can download CodeBlocks IDE, 
or with a g++ compile command: g++ -O2 -Iinclude src/*.cpp -o checker

To use it and check instance A0, you can call: ./checker A0
To use it and check instance A0 and just getting a yes/or answer, you can call: ./checker A0 0

The Checker visualization requires a browser that supports HTML5.


Input / Output
---------------
Inputs:
  - Used instance index (ex: for A0_batch.csv the index is A0).

Note:
  - All csv files (batches, defects, solutions, optimization parameters) should be in the same directory named "instances_checker".

Outputs:
  - X_log.txt file of parsing, constraints verifying trace. X equals to the index of used batch file.
  - X_statistics.csv that contains used batch index under instanceId label, 0 or 1 if solution is invalid or valid under validSolution label, used plates under nPlates label, total loss under totalGeoLoss label and residual width under widthResidual label.

Checker visualization tool takes provided batch & defects files and user generated solution file as inputs and has no output.
Always refresh the html page before reuse the checker visualization.


Checker algorithm Utilization
------------------------------
Once solution file is ready, execute the checker generated binary.
Then:
  - Enter Used Batch csv file index X, the checker will search and parse X_batch.csv, X_solution.csv, X_defects.csv and global_param.csv files in "instances_checker" directory.

The Checker algorithm generates a X_log.txt file automatically.

The Checker algorithm gives the user the following menu.
    Press to verify constraint:
    1 - Items production
    2 - Defects superposing
    3 - Identity and Sequence
    4 - Items dimensions
    5 - All constraints (1-4)
    6 - Display plates wasted & used surface (%)
    0 - Exit
    Pressed key:

For each menu option, the Checker algorithm will display the constraint verification result.
Errors have a "ERROR –" prefix. For example:
      "ERROR -- Node 1: x:500 is not equal to its parent node 1 x:450"

Note:
  A Solution generating errors in its output (terminal and log.txt) is an invalid solution.


Checker visualization Utilization
----------------------------------
Once solution file is ready, open the checker visualization tool in a browser, open solution file and used defects file.
Checker visualization tool will visualize plate 0 and gives the ability to navigate between all solutions used plates.
