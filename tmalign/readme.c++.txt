TM-align: sequence-independent structure alignment of monomer proteins by
TM-score superposition. Please report issues to yangzhanglab@umich.edu
This document is for C++ version of TM-align only.

References to cite:
Y Zhang, J Skolnick. Nucl Acids Res 33, 2302-9 (2005)

DISCLAIMER
Permission to use, copy, modify, and distribute the Software for any
purpose, with or without fee, is hereby granted, provided that the
notices on the head, the reference information, and this copyright
notice appear in all copies or substantial portions of the Software.
It is provided "as is" without express or implied warranty.

==========================
How to install the program
==========================
The following command compiles the program in your Linux computer:

    g++ -static -O3 -ffast-math -lm -o TMalign TMalign.cpp

The '-static' flag should be removed on Mac OS, which does not support
building static executables.

======================
How to use the program
======================
You can run the program without argument to obtain the document.
Briefly, you can compare two structures by:

    TMalign PDB1.pdb PDB2.pdb [Options]

Options:
    -u    TM-score normalized by user assigned length (the same as -L)
          warning: it should be >= minimum length of the two structures
          otherwise, TM-score may be >1

    -a    TM-score normalized by the average length of two structures
          T or F, (default F)

    -i    Start with an alignment specified in fasta file 'align.txt'

    -I    Stick to the alignment specified in 'align.txt'

    -m    Output TM-align rotation matrix

    -d    TM-score scaled by an assigned d0, e.g. 5 Angstroms

    -o    Output the superposition to 'TM_sup*'
            $ TMalign PDB1.pdb PDB2.pdb -o TM_sup
          View superposed C-alpha traces of aligned regions by RasMol or PyMOL:
            $ rasmol -script TM_sup
            $ pymol -d @TM_sup.pml
          View superposed C-alpha traces of all regions:
            $ rasmol -script TM_sup_all
            $ pymol -d @TM_sup_all.pml
          View superposed full-atom structures of aligned regions:
            $ rasmol -script TM_sup_atm
            $ pymol -d @TM_sup_atm.pml
          View superposed full-atom structures of all regions:
            $ rasmol -script TM_sup_all_atm
            $ pymol -d @TM_sup_all_atm.pml
          View superposed full-atom structures and ligands of all regions
            $ rasmol -script TM_sup_all_atm_lig
            $ pymol -d @TM_sup_all_atm_lig.pml

 -fast    Fast but slightly inaccurate alignment by fTM-align algorithm

   -cp    Alignment with circular permutation

    -v    Print the version of TM-align

    -h    Print the full help message, including additional options

    (Options -u, -a, -d, -o will not change the final structure alignment)

Example usages:
    TMalign PDB1.pdb PDB2.pdb
    TMalign PDB1.pdb PDB2.pdb -u 100 -d 5.0
    TMalign PDB1.pdb PDB2.pdb -a T -o PDB1.sup
    TMalign PDB1.pdb PDB2.pdb -i align.txt
    TMalign PDB1.pdb PDB2.pdb -m matrix.txt
    TMalign PDB1.pdb PDB2.pdb -fast
    TMalign PDB1.pdb PDB2.pdb -cp

Additional options only available in C++ version of TMalign:
    -dir     Perform all-against-all alignment among the list of PDB
             chains listed by 'chain_list' under 'chain_folder'. Note
             that the slash is necessary.
             $ TMalign -dir chain_folder/ chain_list

    -dir1    Use chain2 to search a list of PDB chains listed by 'chain1_list'
             under 'chain1_folder'. Note that the slash is necessary.
             $ TMalign -dir1 chain1_folder/ chain1_list chain2

    -dir2    Use chain1 to search a list of PDB chains listed by 'chain2_list'
             under 'chain2_folder'
             $ TMalign chain1 -dir2 chain2_folder/ chain2_list

    -suffix  (Only when -dir1 and/or -dir2 are set, default is empty)
             add file name suffix to files listed by chain1_list or chain2_list

    -atom    4-character atom name used to represent a residue.
             Default is " CA " for proteins
             (note the spaces before and after CA).

    -ter     Strings to mark the end of a chain
             3: (default) TER, ENDMDL, END or different chain ID
             2: ENDMDL, END, or different chain ID
             1: ENDMDL or END
             0: (default in the first C++ TMalign) end of file

    -split   Whether to split PDB file into multiple chains
             0: (default) treat the whole structure as one single chain
             1: treat each MODEL as a separate chain (-ter should be 0)
             2: treat each chain as a seperate chain (-ter should be <=1)

    -outfmt  Output format
             0: (default) full output
             1: fasta format compact output
             2: tabular format very compact output
            -1: full output, but without version or citation information

    -byresi  Whether to assume residue index correspondence between the
             two structures.
             0: (default) sequence independent alignment
             1: (same as TMscore program) sequence-dependent superposition,
                i.e. align by residue index
             2: (same as TMscore -c, should be used with -ter <=1)
                align by residue index and chain ID
             3: (similar to TMscore -c, should be used with -ter <=1)
                align by residue index and order of chain

    -TMcut   -1: (default) do not consider TMcut
             Values in [0.5,1): Do not proceed with TM-align for this
                 structure pair if TM-score is unlikely to reach TMcut.
                 TMcut is normalized is set by -a option:
                 -2: normalized by longer structure length
                 -1: normalized by shorter structure length
                  0: (default, same as F) normalized by second structure
                  1: same as T, normalized by average structure length

    -mirror  Whether to align the mirror image of input structure
             0: (default) do not align mirrored structure
             1: align mirror of chain1 to origin chain2

    -het     Whether to align residues marked as 'HETATM' in addition to 'ATOM  '
             0: (default) only align 'ATOM  ' residues
             1: align both 'ATOM  ' and 'HETATM' residues

    -infmt1  Input format for chain1
    -infmt2  Input format for chain2
            -1: (default) automatically detect PDB or PDBx/mmCIF format
             0: PDB format
             1: SPICKER format
             2: xyz format
             3: PDBx/mmCIF format

Example usage: (search query.pdb against I-TASSER PDB library by fTM-align)
    wget https://zhanglab.dcmb.med.umich.edu/library/PDB.tar.bz2
    tar -xjvf PDB.tar.bz2
    TMalign -dir1 PDB/ PDB/list -suffix .pdb query.pdb -outfmt 2 -fast
Here, "-dir1 PDB/ PDB/list -suffix .pdb" means reading a set of structure
files from the folder "-dir1"; the list of entries are listed by the file
"PDB/list"; a file name extension ".pdb" should be appended to the entry
name to get the full file name. "-fast" means using fast version of TM-align
(fTM-align) to perform the structure alignment. "-outfmt 2" means using
compact tabular output format.
