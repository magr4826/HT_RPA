# High Throughput RPA (HT_RPA)
Code to run High Throughput RPA calculations on the Alexandria database [1].             
This code is an expansion on the workflows we used for High Throughput IPA calculations [2].            

**REQUIREMENTS:**
The code was tested on the makalu4 HPC cluster at the TU Ilmenau, running CentOS 7.                 

All version numbers are the versions that were tested - other version might work, but success is not guaranteed.               
- Python > v3.9
    Installed with miniconda and the following packages (plus their dependencies):
    - numpy = v1.26.4
    - pandas = v2.2.1
    - scipy = v1.13.0
    - matplotlib = v3.8.4
    - pymatgen = v2024.5.1
    - ase = v3.22.1
- QuantumEspresso = v7.1 ("/bin" directory loaded into the path, compiled for MPI)              
- Yambo = v5.2 ("/bin" directory loaded into the path, compiled for MPI)               
- Wannier90 = v3.1 ("/bin" and "/utility" directory loaded into the path, compiled for MPI, used only for the kmesh.pl utility)                

Under the condition that you already have QuantumEspresso, Yambo and Wannier90 installed and properly set up, the installation of the python packages should finish in a few minutes.                    
No further installation is necessary (though changes might have to be made to some scripts, see below).            

**USAGE:**

1. Convert the Alexandria database into .pckl-files using the conversion.ipynb notebook and place them in the input_data folder (already done for this version)
2. Specify the name of the .pckl-file and the workflows to run in the HT_watcher.py script
3. If necessary, adapt the HT_watcher.py script and the start_calc function in basic_utils.py to your cluster
4. Run the HT_watcher.py script
5. If necessary, restart failed calculations using HT_restart.py or HT_full_restart.py

Various QuantumEspresso and Yambo files should appear in the calc folder (with separate folders per compound), and .json files per compound should appear in the database folder (see below).            
Running all calculations likely takes multiple weeks/months (depending on your computational resources), but individiual calculations should finish within hours.            

The runtimechecker.ipynb notebook implements various analyses of the calc folder, e.g., of the total computational time.                
It has been checked on the calc-folder supplied by us, and parts of it might not work on other setups (e.g., if your cluster uses a different format for log-files).              
It should finish in a few minutes.              

**DATABASE:**
The main results of the calculations are saved in .json-files in the database folder.              
The database follows the structure of the original database (Alexandria and IPA).             
The IPA calculations add the following data and parameters:

data:
- ipa_epsI_0:         Imaginary part of the xx-component of the dielectric tensor
- ipa_epsR_0:         Real part of the xx-component of the dielectric tensor
- ipa_indirect_gap:   Indirect gap determined on the converged k-point grid for the IPA calculation
- ipa_direct_gap:     Direct gap determined on the converged k-point grid for the IPA calculation
- ipa_eps_similarity: Similarity coefficient between ipa_epsI_0 using the final and the penultimate k-point density
- ipa_epsI_1:         Imaginary part of the yy-component of the dielectric tensor
- ipa_epsR_1:         Real part of the yy-component of the dielectric tensor
- ipa_epsI_2:         Imaginary part of the zz-component of the dielectric tensor
- ipa_epsR_2:         Real part of the zz-component of the dielectric tensor

The latter four properties are only present if they are different from eps_I_0 and eps_R_0 based on symmetry.                
The gaps differ from the original gaps in the Alexandria database mostly due to the use of different kgrids, used codes, and pseudopotentials.

params:
- pw_conv_k:          Converged K-point density in inverse Angstrom for the ground state calculation
- pw_conv_cutoff:     Converged cutoff in Rydberg for the ground state calculation
- ipa_eps_kppa:       Converged K-point density in inverse Angstrom for the IPA calculation
- ipa_nbands:         Number of bands used for the IPA calculation
- ipa_broad:          Value of the broadening in eV used for the IPA calculation

The RPA calculations add the following data and parameters:

data:
- ipa_epsI:           Imaginary part of the dielectric tensor under the IPA (RPA with G=0)
- ipa_epsR:           Real part of the dielectric tensor under the IPA (RPA with G=0)
- rpa_epsI:           Imaginary part of the dielectric tensor under the RPA
- rpa_epsR:           Real part of the dielectric tensor under the RPA 
- rpa_similarity:     Componentwise SC between rpa_epsI using the final and the penultimate number of G-vectors         
params:
- rpa_kppa:           K-point density used for the RPA calculation
- rpa_G:              G-vector cutoff used for the RPA calculation, in mRy
- rpa_broad:          Broadening parameter used for the RPA calculation, in eV

The dielectric functions and rpa_G are dictionaries containing the diagonal elements of the dielectric tensor/the converged number of G_vectors for that direction.                   
The IPA dielectric functions using the IPA workflow and the RPA workflow might differ slightly as the used k-point grid can vary.

[1] https://alexandria.icams.rub.de/
[2] https://github.com/magr4826/HT_IPA 