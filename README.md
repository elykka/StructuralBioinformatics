# StructuralBioinformatics

Project no. I - Identification of structure clusters in
Molecular Dynamics trajectories from Residue Interaction
Networks

Residue Interaction Networks are derived from proteins structures estimating contacts from distance
measures. RING is a command line tool implemented in C++ which takes in input a PDB file and returns
the list of contacts in a protein complex. RING is able to classify different types of contacts based on
geometrical and physico-chemical properties of the amino acids. Types of contacts are:
● Hydrogen bonds (HBOND)
● Van der Waals interactions (VDW)
● Disulfide bridges (SBOND)
● Salt bridges (IONIC)
● π-π stacking (PIPISTACK)
● π-cation (PICATION)
● Inter-Atomic Contact (IAC), generic contact simply based on distance
RING generates two files. An “edge” file containing the contacts and a “node” file containing the list of
residues which are in contacts and their attributes. A node is identified by a string like this
“A:159:_:PRO”, in which the chain, residue index, insertion code and residue name are column (“:”)
separated. An insertion code equal to “_” indicates there is no insertion code for that residue.
Project goals


The *output of a Molecular Dynamics* (MD) simulation is a trajectory file which describes the change of
atomic coordinates from an initial state to a final state (after a certain amount of time) when a forcefield is
applied. The full trajectory can be described by a subset of intermediate conformations / atomic
coordinates / snapshots. Each team of students is requested to develop a software that identifies
clusters of related (similar) structural conformations (snapshots) in a MD simulation. Additionally, the
software should be able to identify those residues and contacts (pairs of residues) which are relevant to
describe the transition between clusters.
