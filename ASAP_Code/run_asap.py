import time
start = time.time()
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
import AutomaticSurfaceAdsorbateStructureProvider as ASAP
import numpy as np
from pymatgen.core.structure import Structure


############# Fuctions for determineing unique CEs ###########################

# # Function to pull uniq atoms for total equiv_atoms list)
def pull_usites(equiv_atoms):
    uniq = []
    [uniq.append(x) for x in equiv_atoms if x not in uniq]
    return uniq

# # Function to create list of CEs from uniq atoms
def envi(usites,structure,only_cations):
    ce_list = []
    ce = finder.compute_coordination_environments(structure,only_cations=only_cations)
    rm_list = []
    cn_list = []
    for x in usites:
        try:
            ce_sym = ce[x][0]['ce_symbol']
            ce_list.append(ce_sym) 
            Cn = finder.allcg[ce_sym]
            cn_list.append(Cn.coordination)
        except TypeError:
            # This will happen if there is no CE because 1) only_cations=False or 2) an error occured with finding the CE
            rm_list.append(x)
        except IndexError:
            rm_list.append(x)
    # usites are adjusted if no CE is found
    [usites.remove(y) for y in rm_list]
    return ce_list , usites , cn_list

############# Variables used by  ASAP and SpaceGroupAnalyzer #################
 
# Make Adjustements here based on desired structure
cutoff = 0.5        # the cutoff threshold for adsorbate atoms to lattice atoms
bonding_dist = 2.0  # the distance of the bonding adsorbent
bonding_atom = 1    # the atom in the adsorbate that is bonding
dope_list = []      # List of subsitution to be made at each adsorbate site (pass an empty array if no subtitutions are desired)
Struc = Structure.from_file('POSCAR')  # Adsorbent pymatgen structure file
Ads = Structure.from_file('POSCAR_O2') # Adsorbate pymatgen structure file

# If the adsorbate is diatomic and you desire a bindentate bond formation:
# set the diatomic atoms bond distance here, set to "None" if not desired
diatomic = None

############# Find CEs using SpacegroupAnalyzer and LocalGeometryFinder ######

SGA = SpacegroupAnalyzer(Struc)
sym_data = SGA.get_symmetry_dataset()
equiv_atoms= sym_data['equivalent_atoms']
finder = LocalGeometryFinder()

# Run functions to detemrine unique sites and CEs comment this if you are setting these mannually
usites = pull_usites(equiv_atoms)
# NOTE: setting only_cations=False will allow ASAP to treat both the cations and anions as adsorbent sites
ce_list , usites , cns = envi(usites,Struc,only_cations=True)

############# Set CEs Yourself ###############################################
# # Note you can also just set the atoms and thier CE mannually like this:
# # The CE names are the IUPAC noemnclature. 
# # They are listed in the SI of https://doi.org/10.1021/acs.chemmater.7b02766
# # Or can be printer with `finder.allcg.get_geometry_from_IUPAC_symbol`

# usites = [0,12]  # Python is a 0 based code meaning Atom 1 is actually Atom 0
# ce_list = ['T:4', 'O:6'] 
# cns = [4,6]

##############################################################################

# Save CE list and Usites for refrence
np.savetxt('Unique_Atoms.txt',usites,fmt='%1.0f')
np.savetxt('CE_of_Unique_Atoms.txt',ce_list,fmt="%s")

############# Build and run ASAP #############################################

# Build ASAP class
asap = ASAP.ASAP(ce_list,cns,usites,Ads,Struc,bonding_atom,cutoff,bonding_dist,dope_list,diatomic)

# Run and Create all unique adsorabte sites based on usites and ce_list lists
asap.create_all_adsorbates()

############# For Timing #####################################################
end = time.time()
total_time = end - start
print("This took {} sec".format(total_time))
