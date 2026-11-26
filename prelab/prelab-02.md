# Prelab 02: Lecture 3 Summary (Protein Visualization and Structural Analysis)

**Course:** BT338IU – Practice in Bioinformatics | **Instructor:** Dr. Tran Minh Nhựt

## 1. Objectives and overview of protein visualization

Lecture 3 covers protein visualization and structural analysis in `PyMOL`. The lab has three objectives: understanding the four levels of protein structure and composition, learning how physical and chemical properties shape a protein, and using `PyMOL` to display those features in three dimensions.

A protein's function depends on its 3D shape. Inspecting a 3D model lets researchers examine active sites, ligand binding pockets, and conformational changes that are hard to infer from an amino acid sequence alone. Programs such as `PyMOL`, `UCSF Chimera`, and `VMD` turn structural measurements into interactive 3D models used in drug design, protein engineering, and teaching.

These models come from the Protein Data Bank (`RCSB PDB`, `www.rcsb.org`). A `.pdb` file stores atomic coordinates, which are the 3D positions of every atom in a solved structure. Each entry has a unique four-character PDB ID, usually one digit followed by three letters. The lab writes PDB IDs in uppercase so characters like `1` and `L`, or `0` and `O`, do not get mixed up. The first tutorial structure is `1HEW`, hen egg white lysozyme bound to a trisaccharide inhibitor called tri-(N-acetylglucosamine), or tri-NAG. You can download a `.pdb` file from `RCSB PDB` and open it with `File → Open`, or fetch it directly inside `PyMOL`.

## 2. PyMOL setup and graphical interface

The course uses the free Educational version of `PyMOL` (v3.1.5.1 in the manual), available for Windows, macOS, and Linux. Setup requires registering for an educational account, installing the software, and loading the emailed license file through the `Help` menu. `PyMOL` includes both a graphical user interface (GUI) and a text command line for typing commands or running scripts. Lecture 3 starts with the GUI, then combines menu actions with short commands for analysis.

Four main areas surround the 3D workspace. The Action Toolbar holds quick camera and style tools: `Zoom` fits selected atoms to the screen, `Orient` centers the molecule along its principal axes, `Rock` gently oscillates the model for depth perception, and `Preset` applies standard views such as cartoon, sticks, surface, or publication styles. The Panel Toolbar opens four toolsets: `Builder` (to construct or edit chemical groups, amino acids, proteins, nucleic acids, bonds, and charges), `Scenes` (to save camera angles, colors, and representations for movies), `Draw/Ray` (to adjust rendering appearance), and `More`.

To the left of the workspace, the Content Panel lists every loaded object and selection in the order added. Each entry has five menus: `Action (A)`, `Show (S)`, `Hide (H)`, `Label (L)`, and `Color (C)`. At the bottom, the Toggle Toolbar switches mouse controls between `3-Button Viewing` and `3-Button Editing`, opens `Wizards` for guided edits, and toggles three panels: the `Sequence (SEQ)` bar at the top, the `Timeline` panel for animations, and the `Command` line below the workspace.

## 3. Examining the four levels of protein structure

A protein's structure is organized into four levels, and `PyMOL` can isolate each one with a few commands. Primary structure is the linear sequence of amino acids linked by peptide bonds. After loading lysozyme (`fetch 1LYZ`), switching to `cartoon`, and enabling `Display → Sequence`, hovering over the sequence bar shows residue names and positions. Specific regions can be selected and colored, such as `select my-segment, resi 1-10` followed by `color red, my-segment`.

Secondary structure is local folding driven by backbone hydrogen bonding, which produces alpha-helices (`α-helices`), beta-sheets (`β-sheets`), and loops or coils. In `PyMOL`, the `ss` selector colors each type directly: `color red, ss H` for helices, `color yellow, ss S` for sheets, and `color cyan, ss L+` for loops and coils.

Tertiary structure is the overall 3D fold of a single polypeptide chain, held together by disulfide bonds, ionic interactions, and hydrophobic packing. Because disulfide bridges are a major stabilizing force, the tutorial displays cysteine residues (`show spheres, resn CYS`), selects the connecting sticks, and colors them orange.

Quaternary structure is how multiple polypeptide subunits assemble into one complex. Hemoglobin (`1A3N`) has four subunits; coloring `chain A`, `chain B`, `chain C`, and `chain D` in red, blue, green, and yellow separates them visually. In Section 4.5, students identify all four structural levels in five PDB entries: `1A3N`, `1GZX`, `3NIR`, `2MNR`, and `1QTI`.

## 4. Analyzing protein properties, interactions, and 3D alignment

Beyond backbone geometry, Lecture 3 covers five ways to display chemical properties, molecular interactions, and structural comparisons in 3D. To map hydrophobicity onto `1HEW`, you run an external script (`color_h.py`), execute `color_h` (red for hydrophobic areas) or `color_h1` (green), and switch to `Surface` view to see nonpolar patches against polar regions. For electrostatic potential, you submit a PDB file or ID to the `PDB2PQR` web server (`server.poissonboltzmann.org/pdb2pqr`), which returns a `.pqr` file (coordinates with partial charges and atomic radii) and a `.dx` potential grid. Loading both into `PyMOL` and running `ramp_new elec_potential, potential, [-10, 0, 10]` with `set surface_color, elec_potential, protein` projects positive and negative surface charges onto the model.

To study contacts between two protein chains, the tutorial loads `1h4g` (`fetch 1h4g, async=0`), runs `InterfaceResidues.py`, and calls `interfaceResidues 1h4g`. Switching to surface view and removing water (`Action → remove water`) leaves the binding interface unobstructed. Hydrogen bonds are inspected on `1h4g` by deleting duplicate chains (`B` and `D`), selecting chain `C`, and running `Action → find → polar contacts → to other atoms in object`. Showing the involved residues as sticks, hiding non-interacting waters, and displaying bridging waters as non-bonded spheres (`nb-sphere`) isolates the polar network.

Finally, 3D structural alignment superimposes two structures directly. The `1h4g` dimer is split into two monomer `.pdb` files exported with bond-connect records and segment identifiers, then reloaded and superimposed with `align 1h4g-1, 1h4g-2`. In the final exercises (Section 6), students display hydrophobicity, electrostatic potential, and hydrogen bonds for `1A2C`, `1HCL`, `1MOE`, and `4HHB`, plus 3D subunit alignment for `1MOE` and `4HHB`.
