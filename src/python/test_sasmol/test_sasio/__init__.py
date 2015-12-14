'''
General Contract:


contract for pdb files io tests:

1. 1ATM.pdb: a pdb with one atom only
2. 1ATM-1to2.pdb: a pdb of two frames based on the above 1ATM.pdb
3. 2AAD.pdb: a pdb with 2 amino acids
4. 2AAD-1to3.pdb: a pdb of three frames based on the above 2AAD.pdb
5. rna.pdb: a pdb of rna molecule
6. rna-1to10.pdb: a pdb of 10 frames based on the above rna.pdb
7. rna-1to1250_1.0g.pdb: a pdb of 1250 frames with a size of 1.0gb based on the above rna.pdb

contract for dcd files io tests:

1. 1ATM.dcd: two frames based on 1ATM.pdb (one atom)
2. 2AAD.dcd: three frames based on 2AAD.pdb (two amino acids)
3. rna-1to10.dcd: 10 frames based on rna.pdb
4. rna-1.0g.dcd: 1.0gb dcd based on rna.pdb

For rna tests, the resulting coordinates are partially selected to compare against the expected coordinates.

All original pdb/dcd files are located at ./data/
All generated pdb/dcd files by this test suite are located at ./data/test-results/

Note: Testing code should limited to files < 1.0 GB and thus files this size and smaller should be the qualifier for a passing score. Files > 1 GB are not supported and thus no tests should be noted to fail. (see ticket #86 for details)
'''

