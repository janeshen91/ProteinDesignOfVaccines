##checking the RMSD at the epitopes
cd /OneDrive\ -\  Inc/202409_RFDiffusion
/Applications/PyMOL.app/Contents/bin/pymol -cq ./scripts/pymolAlignToantibodyepitope.py \
--ref_obj ./antibody/HADDOCK3/antibody_epitopeOnly.pdb \
--ref_chain "A" \
--ref_motif "235-236+307-309+478-483+510-511" \
--input ./antibody/prompt1PDBFilesAndMotifs.csv \
--output_csv ./antibody/prompt1PDBFilesAndMotifsResults.csv \
--outdir ./antibody/prompt1_aligned_structures

/Applications/PyMOL.app/Contents/bin/pymol -cq ./scripts/pymolTryAlign.py

/Applications/PyMOL.app/Contents/bin/pymol -cq ./scripts/pymolStuff.py -- \
--ref_obj ./antibody/HADDOCK3/antibody_epitopeOnly.pdb \
--ref_sel "chain A and resi 235-236+307-309+478-483+510-511" \
--input ./antibody/prompt1Rest.txt \
--output_csv ./antibody/prompt1_rmsd_output.csv

/Applications/PyMOL.app/Contents/bin/pymol -cq ./scripts/pymolStuff.py -- \
--ref_obj ./antibody/HADDOCK3/antibody_epitopeOnly.pdb \
--ref_sel "chain A and resi 235-236+307-309+478-483+510-511" \
--input ./antibody/bestPDBs/prompt2iiPDBsFilesAndMotifs.txt \
--output_csv ./antibody/bestPDBs/prompt2ii_rmsd_output.csv

/Applications/PyMOL.app/Contents/bin/pymol -cq ./scripts/pymolStuff.py -- \
--ref_obj ./antibody/HADDOCK3/antibody_epitopeOnly.pdb \
--ref_sel "chain A and resi 235-236+307-309+478-483+510-511" \
--input ./antibody/bestPDBs/prompt2PDBFilesAndMotifs.txt \
--output_csv ./antibody/bestPDBs/prompt2_rmsd_output.csv


##last one

/Applications/PyMOL.app/Contents/bin/pymol -cq ./scripts/pymolStuff.py -- \
--ref_obj ./antibody/HADDOCK3/antibody_epitopeOnly.pdb \
--ref_sel "chain A and resi 235-236+307-309+478-483+510-511" \
--input ./antibody/bestPDBs/prompt1iiPDBsFilesAndMotifs.txt \
--output_csv ./antibody/bestPDBs/prompt1ii_rmsd_output.csv


/Applications/PyMOL.app/Contents/bin/pymol -cq ./scripts/pymolStuff.py -- \
--ref_obj ./antibody/HADDOCK3/antibody_epitopeOnly.pdb \
--ref_sel "chain A and resi 235-236+307-309+478-483+510-511" \
--input ./antibody/bestPDBs/prompt1ii_20FilesAndMotifs.txt \
--output_csv ./antibody/bestPDBs/prompt1ii_20_rmsd_output.csv

/Applications/PyMOL.app/Contents/bin/pymol -cq ./scripts/pymolStuff.py -- \
--ref_obj ./antibody/HADDOCK3/antibody_epitopeOnly.pdb \
--ref_sel "chain A and resi 235-236+307-309+478-483+510-511" \
--input ./antibody/prompt1PDBFilesAndMotifs.txt \
--output_csv ./antibody/prompt1_LAST_rmsd_output.csv

/Applications/PyMOL.app/Contents/bin/pymol -cq ./scripts/pymolStuff.py -- \
--ref_obj ./antibody/HADDOCK3/antibody_epitopeOnly.pdb \
--ref_sel "chain A and resi 235-236+307-309+478-483+510-511" \
--input ./antibody/prompt1ii_last_PDBFilesAndMotifs.txt \
--output_csv ./antibody/prompt1ii_LAST_rmsd_output.csv


/Applications/PyMOL.app/Contents/bin/pymol -cq ./scripts/pymolStuff.py -- \
--ref_obj ./antibody/HADDOCK3/antibody_epitopeOnly.pdb \
--ref_sel "chain A and resi 235-236+307-309+478-483+510-511" \
--input ./antibody/prompt2ii_last_PDBFilesAndMotifs.txt \
--output_csv ./antibody/prompt2ii_LAST_rmsd_output.csv

/Applications/PyMOL.app/Contents/bin/pymol -cq ./scripts/pymolStuff.py -- \
--ref_obj ./antibody/HADDOCK3/antibody_epitopeOnly.pdb \
--ref_sel "chain A and resi 235-236+307-309+478-483+510-511" \
--input ./antibody/prompt2_last_PDBFilesAndMotifs.txt \
--output_csv ./antibody/prompt2_LAST_rmsd_output.csv