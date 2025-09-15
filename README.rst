Surface-MSD: Mean Square Displacement at Liquid Interfaces
==========================================================

This package provides tools to compute the **mean square displacement (MSD)** 
of molecules restricted to the **surface molecular layer** of a liquid. 
It is based on intrinsic surface analysis methods (e.g., ITIM, GITIM) and 
was originally developed to support the analysis in 
*Fábián et al., J. Phys. Chem. B 2017, 121, 5582−5594* 
DOI: 10.1021/acs.jpcb.7b02220

The code allows you to analyze the mobility of interfacial molecules directly 
from molecular dynamics trajectories, either through a command-line script 
or programmatically via Python classes.

Installation
------------

To install the package locally from source:

.. code-block:: bash

    pip install --user .

Usage
-----

There are **two ways** to use the package: through the command-line script or 
via the Python API.

Command-line script
~~~~~~~~~~~~~~~~~~~

The package installs the script ``layer-msd``, which can be run directly 
on your trajectories.

**Example:**

.. code-block:: bash

    layer-msd config.gro traj.xtc SOL --output surface_msd.dat --algorithm ITIM --direction xy --pytim-select='name OW' --pytim-params 'alpha=1.5, cluster_cut=3.5'

This computes the in-plane (``xy``) MSD of the residue type ``SOL`` (e.g., water molecules) 
at the liquid surface computed from water oxygen atoms (``'name OW'``), specifying the 
``alpha`` and ``cluster_cut`` pytim parameters and and writes the result to ``surface_msd.dat``.

Python API
~~~~~~~~~~

You can also use the main classes in Python for more flexible analysis.

**Example:**

.. code-block:: python

    import MDAnalysis as mda
    from surface_msd import LayerMSD
    import pytim

    u = mda.Universe("config.gro", "traj.xtc")
    group = u.select_atoms("resname SOL")

    msd = LayerMSD(group, direction="xy", npoints=200)

    # detect interfacial molecules with ITIM
    inter = pytim.ITIM(u, group=u.select_atoms("name OW"))

    for ts in u.trajectory:
        msd.sample(inter.atoms)

    print(msd.data)

Command-line options
--------------------

The script accepts the following flags:

**Positional arguments:**

- ``config``: Input configuration file (e.g., ``.gro``).
- ``traj``: Input trajectory file (``.xtc`` or ``.trr``).
- ``residue``: Residue name to compute the MSD for.

**Optional arguments:**

- ``--output <file>``: Output file for MSD data (default: ``msd.dat``).
- ``--algorithm <ITIM|GITIM>``: Surface detection algorithm (default: ``ITIM``).
- ``--pytim-select <selection>``: Atom selection string for surface analysis (default: ``"name OW"``).
- ``--pytim-params <params>``: Comma-separated list of parameters for the chosen algorithm.
- ``--direction D``: A combination of any of the ``xyz`` components to include in MSD 
  (e.g., ``xy``, ``xz`` , ``z``, ``xyz``).
- ``--length <N>``: Window size (in timesteps) over which MSD is sampled (default: 100).
- ``--start <frame>``: First frame to analyze (default: 0).
- ``--stop <frame>``: Last frame to analyze (default: until trajectory end).
- ``--stride <N>``: Frame stride for reading trajectory (default: 1).
- ``--verbose``: Print progress information.

Classes and API
---------------

The main classes are defined in ``surface_msd.msd_core`` and exposed via ``surface_msd``:

- **``LayerMSD(refgroup, direction="xy", npoints=100)``**  
  Computes the mean square displacement of the center of mass of residues in the surface layer.

  - ``refgroup``: An ``MDAnalysis.AtomGroup`` of residues to follow.
  - ``direction``: Which components to compute (``xy``, ``z``, ``xyz``).
  - ``npoints``: Length of the time window for MSD statistics.

  **Methods:**
  
  - ``sample(inp)``: Sample MSD using the atoms ``inp`` currently identified as interfacial.

  **Attributes:**
  
  - ``data``: Array with MSD values.
  - ``norm``: Normalization factors.

- **``ResidueCenterOfMassPosition``**  
  Helper observable for computing center-of-mass positions of residues.

These classes allow you to integrate MSD analysis into larger workflows using 
``pytim`` and ``MDAnalysis``.

