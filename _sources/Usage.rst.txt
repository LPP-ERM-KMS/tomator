Usage
=====
Tomator 1D is a standalone binary which may be 
interfaced using python through software located in the 'gui' folder. Their functionalities are discussed below

Defining and running a simulation: SimulationInterface.py
---------------------------------------------------------

We will describe the software in this section but for those interested there is
also a video introduction that can be found here: https://youtu.be/qyPUr26huhY

To run the software you need a python venv containing pandas. 
With this environment you may run the software as::

    python SimulationInterface

This will show an interface with two options:

* Define Parameters
* Select Pre-defined Parameters

The parameters are defined in a json file which may be created using the first option or loaded using
the second option, example files are located in examples/InputFiles. When clicking on "define parameters" you will be presented with a dialog where you may specify the various parameters such as:

* Main discharge
    * Magnetic field
        Define the tangential, vertical and horizontal magnetic field size

    * Toroidal machine geometry
        Define the major radius R, and the two minor radii a and b. You may also define
        the limiter extension and the HFS (High Field Side), at the LFS and how many there are
        (nlimitors). Finally specify the plasma volume (vpl)

    * Neutral pressure
        The pressures of the various gasses (Helium and dihydrogen)

* Coupled (RF) power
    Coupled power functions specify where the power goes, for more info and how
    to add your own see :doc:`RF power`

    * RF power
        Specify the amount of power (Prf), it's frequency (freq) and the time it takes
        to ramp up the power (dtpramp)

    * Types of RF power
        * proptone
            deposits power proportional to the electron density, it is also the simplest example.

        * nefix
            Tunes the coupled power for a fixed density at IC

        * fixpowerfrac
            Couples a fixed fraction of the launched power

        * blhr
            Couples a fixed fraction of the launched power (lower hybrid resonance)

        * bTOMAS
            deposition on ECR, built to fit the TOMAS machine, using the budden method.

        * bmanuel
            Deposition on ECR, done by Manuel

    * ECH
        * general EC
            Define parameters which are considered "generally applicable" such as
            the position of power deposition Rdep, the power of EC that was
            absorbed pecabs0, the width of the ECH widthech, the background ECH
            echbackground, the harmonic number (e.g 2) and the losses to the wall
            muw.

        * necfix
            parameters relevant if nefix was chosen prior, the parameters are the vertex of interest 'ic', 
            the fixed density necfix, the initial power Pini and a tunable error variable tauP.

        * TOMAS
            Multiple power locations and the corresponding fractions.

    * ICH
        * general IC
            As with general EC, parameters which are considered "generally applicible" which here is the
            antenna position Rant, the antenna vacuum resistance avlr and wether to use it bantlr.

        * LHR
            Parameters concerning the blhr type, the fraction of power going to the electrons fracpne,
            the fraction going to the lower hybrid fraclhr, the with of the lower hybrid widthlhr, the
            background lower hybrid and HtoHD, how much power gets transferred from hydrogen to deuterium.

    * OTHER
        a fixed temperature for which a csv file containing the dependency may be given.

* Transport
    * Diffusion
        Here you may choose between bDfix,bDbohm and bDscaling as outlined in the original paper,
        defining either the diffusion value if fixed or the factor if bohm or scaling.
    * Convection
        Similar to diffusion, here you may choose inbetween fixed vertical diffusion or scaling 
        diffusion with veq being which equation you choose for this (5 or 8).
    * Tune D and V
        This is the opposite of the previous two options, instead of defining
        the diffusion and convection directly, you tune them so at vertex il
        the electron density should be a specified value and same with ir

* Conditions
    * Physics to include
        All the molecular physics mentioned in the paper, hydrogen, dihydrogen,
        helium, ions, charge exchange, elastic scattering, coulomb scattering,
        impurities, transport ions, transport neutrals, edge particles, poloidal
        drift, vertical drift and collisionality.
    * Initial conditions
        From what the simulation is started, such as the various densities of the species
        present in the plasma and the initial temperature and densities of the electrons.
    * Edge conditions
        if edge physics was selected in the physics tab this governs the edge physics equations
        from the paper.

* Simulation control settings
    * Simulation grid
        The amount of meshpoints
    * Input file
        Continue on from a previous tomator simulation, this might be useful as an initial simulation
        most often starts from non-physical initial conditions as defined by the user, if you start
        from a possible scenario the stable solution arrives faster.
    * Time step
        the initial time, at what time the main ends (if no -t flag was given directly to the binary)
        and the parameters needed for adaptive time stepping
    * Time step for RF coupling
        As RF simulations can be slow, it is useful to run them in a different time scenario.
    * Output parameters
        Specify to either save every N loops or every dtsave seconds.
    * Solver parameters
        The simulation will adaptively refine until the difference between
        densities is less then this tolerance.

You may also modify an existing json file by first choosing "Load Defaults", then modifying the required variables. Having filled in all the necessary entries you may save the inputs with "Save to JSON".
Which may then be loaded with "Select Pre-defined Parameters" and ran.

The output will be generated in a folder Data, one layer above the folder in
which you started the simulation, as you ran the gui this will be in the tomator folder.

Running a simulation: using the binary directly
-----------------------------------------------

The binary takes as a required argument the json file, a simulation can thus be
ran in the terminal as::

    Tomator1D TCV5151X_fixneDV.json

Overview of the simulation: PlotterInterface.py
-----------------------------------------------

A simulation may take quite some time, to track the progress a python script called PlotterInterface
was developed, which also has a video introduciton: https://youtu.be/1ATl7nQellM, you run it as::

    python PlotterInterface

Now you have the option "Plot Simulation" which you may direct to the generated csv file in 
Data/yourjsonfilename/, having done this a browser will open showing the current status of your
simulation, to terminate this plotting server click on "Terminate Server", whereby you are
given a list of active servers you may terminate.
