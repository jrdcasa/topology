Usage Guide
===========

This guide shows common usage patterns for ``topology``.

Basic Example
-------------

.. code-block:: python

    from topology import Molecule

    mol = Molecule.from_xyz('water.xyz')
    print(mol.atoms)
    print(mol.bonds)

Loading molecular data
----------------------

Explain loading methods here…

Performing topology analysis
----------------------------

Explain typical workflows…

More advanced examples appear in :doc:`examples`.

