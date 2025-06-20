Coupled Power Functions
=======================

A coupled power function sets the RF power at every mesh point for a particular species.
In general it does this basing itself on the electron density, examples are included with
the tomator code and may also be created by the user.

Adding a New Coupled Power Function
-----------------------------------

You may want to add your own coupled power function, to do this: 

#. Step 1: Create the Coupling Function File 
    Create a C++ source file named `coupledB<name>.cpp` and include the `coupledpower.h` header,
    then define your coupling functions within this file. 

#. Step 2: Update `coupledpower.cpp` 
    Integrate your function by adding an if-statement to call it based on the `<name>` boolean identifier. 

#. Step 3: Define the New Boolean Parameter 
    Declare a boolean parameter in `simparam.h` and set its default value to `false` in `simparam.cpp`. 

#. Step 4: Modify the Python Interface 
    Add a new entry for your function's boolean in `ChooseParameters.py` -
    `Type` group, and update `self.type_to_parameters_mapping` with the
    parameters that aren't important for your function. 

#. Step 5: Configure the JSON File Using the Simulation Interface 
    Use `SimulationInterface.py` to generate a JSON configuration file that includes your new parameter. 

#. Step 6: Update the Tomator Code 
    In `extractor.cpp`, add an if-condition to the `extract_type` function for your parameter. 
