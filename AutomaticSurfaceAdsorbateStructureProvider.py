import operator
from pymatgen.core.sites import PeriodicSite
from pymatgen.io.vasp.inputs import Poscar
import numpy as np
import math
from sympy.geometry import Plane, Point3D
import sympy as sym    
from sympy import re
from itertools import combinations


class ASAP():
    
    # ASAP can currently handle all coordnation enviroments >= 2
    # Still need to test the TL:3 Trigonal Plane and S:4 Square Plane environmetns
    
    def __init__(
            self,
            CEs: list,
            CNs,
            Usites: list,
            Ads,
            Struc,
            bonding_atom: int,
            cutoff: float,
            bonding_dist: float,
            dope_list: list,
            diatomic: float = None ):
        
        self.Usites = Usites             # list of unique sites in the adsorbing structer
        self.CNs = CNs                   # Cordnation numbers ascoited with the CEs
        self.CEs = CEs                   # list of coordnatiosn environemets for each Usite
        self.Adsorbate = Ads             # geometry file of adsorbate
        self.Adsorbent = Struc           # geometry file of adsorbign structure
        self.bonding_atom = bonding_atom # the atom in the adsorbate that is bonding
        self.bonding_dist = bonding_dist # the distance of the bonding adsorbent
        self.cutoff = cutoff             # the cutoff threshold for adsorbate atoms to lattice atoms
        self.dope_list = dope_list       # List of subsitution to be made at each adsorbate site (pass an empty array if no subtitutions are desired)
        self.diatomic = diatomic         # Bond distance for a diatomic bindentate bond, if inlcuded
    
    def get_neighbors(self):
        nn = []
        [nn.append(self.Adsorbent.get_neighbors(self.Adsorbent.sites[x],5.0)) for x in self.Usites]
        return nn
    
    def create_all_adsorbates(self):
        spot = 1            
        NNs = self.get_neighbors()                                                                                    
        for counter,x in enumerate(self.CEs):         
                if x == 'L:2' or x == 'A:2':
                    inserts, non_per = self.L2_A2(self.Usites[counter],NNs[counter],self.CNs[counter])
                    if self.diatomic:
                        spot = self.get_bidentate_adsorbate(inserts,non_per,self.Usites[counter],spot)
                    else:
                        spot = self.get_adsorbate(inserts,non_per,self.Usites[counter],spot)
                    continue
                else:
                    inserts, non_per = self.New_Sites(self.Usites[counter],NNs[counter],self.CNs[counter])
                    if self.diatomic:
                        spot = self.get_bidentate_adsorbate(inserts,non_per,self.Usites[counter],spot)
                    else:
                        spot = self.get_adsorbate(inserts,non_per,self.Usites[counter],spot)
                
    def dope_me(self,usite,spot,name):
        og_species = self.Adsorbent[usite].species_string
        str_name = '{}{}'.format(spot,name)
        for d in self.dope_list:
            if (d != og_species):
                self.Adsorbent.replace(usite,d)
                pos = Poscar(self.Adsorbent.get_sorted_structure())
                pos.write_file('POSCAR_{}_Site_{}'.format(self.Adsorbent[usite].species_string,str_name))
            else:
                pass
        self.Adsorbent.replace(usite,og_species) 
    
    def New_Sites(self,usite,nn,cn):

        X = self.Adsorbent[usite].x                                                      
        Y = self.Adsorbent[usite].y                                                      
        Z = self.Adsorbent[usite].z                                                  
        sort_nn = sorted(nn, key=operator.itemgetter(1))                            
        neighbors = sort_nn[:cn] 
        
        vectors = []
        for f in range(cn):        
            coords = neighbors[f][0].coords          
            V = coords - [X,Y,Z]                 
            vectors.append(V) 
        vectors = np.asarray(vectors)

        combos = list(combinations(range(cn),3))
        combos = np.asarray(combos)
        sites = []
        n_cords = []
        for c in combos:
            vec_check = [] 
            n_cords = []
            Vj = [0, 0, 0]
            for x in c:
                Vj += vectors[x]
                vec_check.append(vectors[x])
                n_cords.append(neighbors[x].coords)
                
            vec_check = np.asarray(vec_check)
            
            v1 = n_cords[0] - n_cords[1]
            v2 = n_cords[0] - n_cords[2]
        
            cp1 = np.cross(v2, v1)            
            cp1 = cp1*(np.linalg.norm(Vj)/np.linalg.norm(cp1))
            
            angle1 = np.arccos(np.dot(Vj,cp1)/(np.linalg.norm(Vj)*np.linalg.norm(cp1)))
            angle1 = np.degrees(np.abs(angle1))

            newsite = 'bad'
            if np.abs(angle1) < 10 or np.abs(angle1-180) < 10 :
                newsite = 'good'
            Ang = np.zeros((len(vectors),))
            for count,V in enumerate(vectors):
                Ang[count] = np.degrees(np.arccos(np.dot(Vj,V)/(np.linalg.norm(Vj)*np.linalg.norm(V))))
            # For large cn CEs the diff in angles can be <10 and still be a wrong Vj
            if (Ang <= 8).any():
                newsite =' bad'
            if (np.linalg.norm(Vj) < vec_check.all()):
                newsite =' bad'
            if newsite == 'good':
                # adjust length of Vj to desired bond distance
                factor_adjust = self.bonding_dist/np.linalg.norm(Vj)
                Vj = factor_adjust*Vj
                # add back adsorbent site XYZ for cartesien points of adsorbate
                sites.append(Vj+[X,Y,Z])
   
        combos = list(combinations(range(cn),4))
        combos = np.asarray(combos)
        for c in combos:
            vec_check = []
            n_cords = []
            Vj = [0, 0, 0]
            for x in c:
                Vj += vectors[x]
                vec_check.append(vectors[x])
                n_cords.append(neighbors[x].coords)
                
            vec_check = np.asarray(vec_check)
            
            v1 = n_cords[0] - n_cords[1]
            v2 = n_cords[0] - n_cords[2]
            v3 = n_cords[0] - n_cords[3]
       
            cp1 = np.cross(v1, v2)
            cp1 = cp1*(np.linalg.norm(Vj)/np.linalg.norm(cp1))
            cp2 = np.cross(v1, v3)
            cp2 = cp2*(np.linalg.norm(Vj)/np.linalg.norm(cp2))
            angle1 = np.arccos(np.dot(Vj,cp1)/(np.linalg.norm(Vj)*np.linalg.norm(cp1)))
            angle2 = np.arccos(np.dot(Vj,cp2)/(np.linalg.norm(Vj)*np.linalg.norm(cp2)))
            
            angle1 = np.degrees(angle1)
            angle2 = np.degrees(angle2)

            newsite = 'bad'
            if (np.abs(angle1) < 10 or np.abs(angle1-180) < 10) and (np.abs(angle2) < 10 or np.abs(angle2-180) < 10) :
                newsite = 'good'
            
            Ang = np.zeros((len(vectors),))
            for count,V in enumerate(vectors):
                Ang[count] = np.degrees(np.arccos(np.dot(Vj,V)/(np.linalg.norm(Vj)*np.linalg.norm(V))))
            # For large cn CEs the diff in angles can be < 10 and still be a wrong Vj
            # Thus we impose a check to ensure that Vj is not close to any vector
            # inline with a NN 
            if (Ang <= 8).any():
                newsite =' bad'        
            if (np.linalg.norm(Vj) < vec_check.all()):
                newsite =' bad'
            if newsite == 'good':
                # adjust length of Vj to desired bond distance
                factor_adjust = self.bonding_dist/np.linalg.norm(Vj)
                Vj = factor_adjust*Vj
                # add back adsorbent site XYZ for cartesien points of adsorbate
                sites.append(Vj+[X,Y,Z])
        SITES = []
        for s in sites:
            SITES.append(PeriodicSite('O',s,self.Adsorbent.lattice,to_unit_cell=True,coords_are_cartesian=True))
        inserts = SITES
        non_per = sites                
        return inserts, non_per

    def L2_A2(self,usite,nn,cn):
        # Liner CE, we find sites around thins CE usisng a rotaion around the 
        # adsorabte site in a plane that is normal to the neighboring atoms that
        # make up the CE
        X = self.Adsorbent[usite].x                                                      
        Y = self.Adsorbent[usite].y                                                      
        Z = self.Adsorbent[usite].z                                                 
        sort_nn = sorted(nn, key=operator.itemgetter(1))                            
        neighbors = sort_nn[:cn]
        vectors = []
        for f in range(cn):        
            coords = neighbors[f][0].coords          
            V = coords - [X,Y,Z]               
            vectors.append(V) 
        vectors = np.asarray(vectors)
        pie = math.pi
        val = [0, pie/4,pie/2,pie*3/4,pie,pie*5/4,pie*6/4,pie*7/4]#,-pie]
        if (np.abs(np.cross(vectors[0],vectors[1])) < .1).all():  #L2
            # If the corss product is close to zero for the two corridnated 
            # atoms then the CE is L2 
            t = sym.Symbol('t', real=True)
            # Create plane through unique site with normal vector = the new oxygen site
            pln = Plane(Point3D(0,0,0), normal_vector=vectors[0])
            # Create arbitrary point OBJECT as function of t from new plane
            pln_ap = pln.arbitrary_point(t)
            sites = []                             
            axis = []                                                           
            for x in val:                                                     
                axis.append(pln_ap.subs(t,x))
            axis = np.asarray(axis,dtype=float)
            for y in axis:
                factor_adjust = self.bonding_dist/np.linalg.norm(y)
                Vj = factor_adjust*y
                # add back adsorbent site XYZ for cartesien points of adsorbate
                sites.append(Vj+[X,Y,Z])
            SITES = []
            for s in sites:
                SITES.append(PeriodicSite('O',s,self.Adsorbent.lattice,to_unit_cell=True,coords_are_cartesian=True))
            inserts = SITES
            non_per = sites
            
        else:  #A2
            pln1 = Plane(Point3D(0,0,0), normal_vector=vectors[0])
            pln2 = Plane(Point3D(0,0,0), normal_vector=vectors[1])
            # Create arbitrary point OBJECT as function of t from new plane
            t = sym.Symbol('t', real=True)
            pln_ap1 = pln1.arbitrary_point(t)
            pln_ap2 = pln2.arbitrary_point(t)
            sites = []                             
            axis = []                                                           
            for x in val:                                                       
                ax1 = pln_ap1.subs(t,x)
                ax2 = pln_ap2.subs(t,x)
                xo = (ax1[0] + ax2[0])/2
                yo = (ax1[1] + ax2[1])/2
                zo = (ax1[2] + ax2[2])/2
                axis.append([xo,yo,zo])
            axis = np.asarray(ax2,dtype=float)
            for y in axis:  
                factor_adjust = self.bonding_dist/np.linalg.norm(y)
                Vj = factor_adjust*y
                # add back adsorbent site XYZ for cartesien points of adsorbate
                sites.append(Vj+[X,Y,Z])
            SITES = []
            for s in sites:
                SITES.append(PeriodicSite('O',s,self.Adsorbent.lattice,to_unit_cell=True,coords_are_cartesian=True))
            inserts = SITES
            non_per = sites
        return inserts, non_per
        
    def make_pln(self,usite,vec):
        a = vec[0]
        b = vec[1]
        c = vec[2]
        #  to make plane go therough new site and not adsorbent
        xo = self.Adsorbent[usite].x + a
        yo = self.Adsorbent[usite].y + b
        zo = self.Adsorbent[usite].z + c
        t = sym.Symbol('t', real=True)
        # Create plane through unique site with normal vector = the new oxygen site
        pln = Plane(Point3D(xo,yo,zo), normal_vector=(a,b,c))
        # Create arbitrary point OBJECT as function of t from new plane
        ap = pln.arbitrary_point(t)
        return ap
    
    # we must loop through t = 0 to pi to find best spot of arbitrary line rotation
    def arbitraty_axis(self,ap,val,vec,usite):
        a = vec[0]
        b = vec[1]
        c = vec[2]
        xo = self.Adsorbent[usite].x + a
        yo = self.Adsorbent[usite].y + b
        zo = self.Adsorbent[usite].z + c
        t = sym.Symbol('t', real=True)
        axis = [re(ap.x.subs(t, val)) - xo, re(ap.y.subs(t, val)) - yo, re(ap.z.subs(t, val)) - zo]
        return axis
    
    def rotation_matrix(self,axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        """
        #dot = sym.Function('dot')
        axis = np.asarray(axis)
        axis = axis / math.sqrt(np.dot(axis, axis))
        a = math.cos(theta / 2)
        b, c, d = -axis * math.sin(theta / 2)
        aa, bb, cc, dd = a * a, b * b, c * c, d * d
        bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
        return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                         [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                         [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])
    
    
    def dist_from_lattice(self,period_sites):                                                  
        dist = 0
        for x in period_sites:                                                   
            self.Adsorbent.append(x.species_string,x.coords,coords_are_cartesian=True)
            nn = self.Adsorbent.get_neighbors(self.Adsorbent.sites[-1],5.0)                       
            sort_nn = sorted(nn, key=operator.itemgetter(1))                            
            try:                                                                        
                dist += sort_nn[0][1]                                                  
            except IndexError:                                                          
                dist += 5.0                                                            
            self.Adsorbent.remove(self.Adsorbent.sites[-1])                                        
        return dist 
    
    
    def shortest_dist_from_lattice(self,period_sites):                                                  
        dist = 100
        for x in period_sites:                                                   
            self.Adsorbent.append(x.species_string,x.coords,coords_are_cartesian=True)
            nn = self.Adsorbent.get_neighbors(self.Adsorbent.sites[-1],5.0)                       
            sort_nn = sorted(nn, key=operator.itemgetter(1))  
            try:
                if sort_nn[0][1] < dist:
                    dist = sort_nn[0][1] 
            except IndexError:                                                          
                dist += 5.0                                                                                                          
            self.Adsorbent.remove(self.Adsorbent.sites[-1])                                                                                           
        return dist 
    
    
    def get_adsorbate(self,inserts,non_per,usite,spot):                       
        start_spot = spot  
        remove = []  
        remove_non_per = []
        for counter,site in enumerate(inserts):
            for K,y in enumerate(self.Adsorbent):                                                   
                if (K != usite) and (site.distance(y) < self.cutoff): # 2.4 is a cutoff 
                    remove.append(site)
                    remove_non_per.append(non_per[counter])
                    break
        new_inserts = [x for x in inserts if x not in remove]
        new_non_per = [x for x in non_per if x not in np.asarray(remove_non_per)]
        sorb_vecs,adsorbate_types = self.adsorbate_vectors()
        if len(sorb_vecs) > 1:                               
            for counter,site in enumerate(new_inserts): # enumerte is used so .remove can be used
                non_p_site = new_non_per[counter]                                                                                       
                pie = math.pi                                                       
                val = [0, pie/4, pie/2, pie*3/4, pie, -pie/4, -pie/2, -pie*3/4]#,-pie]               
                vec = non_p_site - self.Adsorbent[usite].coords                          
                pln_ap = self.make_pln(usite,vec)                              
                axis = []                                                           
                for x in val:                                                       
                    axis.append(self.arbitraty_axis(pln_ap, x, vec, usite))        
                dist = 0  
                rot_angles = np.asarray(np.linspace(0,360,9,dtype='int'))
                np.delete(rot_angles,-1) # dont need 360 degree rotations
                for arb_rot in rot_angles: 
                    for angle in rot_angles:
                        # count10 = 0  # used for vizualizing all rotations
                        # First rotaion matrix is a rotation around the vec
                        rotate_mx = self.rotation_matrix(vec, np.radians(angle))
                        for x in axis:
                            atom_spots = []
                            # Second rotaion is around an arbitrary axis through
                            # the new adsorbate site 
                            rotate_mx2 = self.rotation_matrix(x, np.radians(arb_rot))
                            for z in sorb_vecs:
                                # z is the arbitrary axis (do not add vec here)
                                first_rot = np.dot(rotate_mx2,z)
                                atom_spots.append(np.dot(rotate_mx,first_rot+vec))
                            period_sites = []
                            for count,new in enumerate(atom_spots):
                                atom = np.array(new, dtype=np.float64) + self.Adsorbent[usite].coords
                                ATOM = PeriodicSite(adsorbate_types[count],atom,self.Adsorbent.lattice,to_unit_cell=True,coords_are_cartesian=True)
                                period_sites.append(ATOM)
                            ##################################################    
                            # # Strickly for visualizing all rotation of adsorbate tested
                            # for g in period_sites:                                                  
                            #     self.Adsorbent.append(g.species_string,g.coords,coords_are_cartesian=True)
                            # pos = Poscar(self.Adsorbent)                                               
                            # pos.write_file('POSCAR_{}_Site_{}_{}_{}_{}'.format(self.Adsorbent[usite].species_string,spot,arb_rot,angle,count10))
                            # count10 += 1
                            # for h in period_sites: # just usign x as a counter 
                            #     self.Adsorbent.remove(self.Adsorbent.sites[-1]) 
                            ##################################################
                            tot = self.dist_from_lattice(period_sites)                    
                            if (tot > dist):                                                
                                best_sites = period_sites                                                                           
                                dist = tot
                # # check one more time to ensure adsrobate is not sitting too cloase to any oher atoms
                final_check = self.shortest_dist_from_lattice(best_sites) 
                if final_check > 0.5:
                    for x in best_sites:                                                  
                        self.Adsorbent.append(x.species_string,x.coords,coords_are_cartesian=True)
                        pos = Poscar(self.Adsorbent)                                             
                        pos.write_file('{}_{}_Site_{}'.format('POSCAR',self.Adsorbent[usite].species_string,spot))
                        if len(self.dope_list) > 0:
                            self.dope_me(usite,spot,'') # Create same site but doped      
                    for x in best_sites: # just usign x as a counter 
                        self.Adsorbent.remove(self.Adsorbent.sites[-1])                                                            
                else:
                    pass   
                spot = spot + 1                                               
            if len(self.dope_list) > 0:
        	    self.dope_me(usite,start_spot,'_{}_no_O2'.format(int(spot)-1)) #creates doped usite with no O2
        elif len(sorb_vecs) == 1:
            for x in new_inserts:                                                  
                self.Adsorbent.append(adsorbate_types[0],x.coords,coords_are_cartesian=True)
                pos = Poscar(self.Adsorbent)                                             
                pos.write_file('{}_{}_Site_{}'.format('POSCAR',self.Adsorbent[usite].species_string,spot))  
                if len(self.dope_list) > 0:
                    self.dope_me(usite,spot,'') # Create same site but doped
                spot = spot + 1
                self.Adsorbent.remove(self.Adsorbent.sites[-1])       
            if len(self.dope_list) > 0:
                self.dope_me(usite,start_spot,'_{}_no_O2'.format(int(spot)-1)) #creates doped usite with no adsorbate                                                  
        return spot 
    
    def get_bidentate_adsorbate(self,inserts,non_per,usite,spot):                       
        start_spot = spot  
        remove = []  
        remove_non_per = []
        for counter,site in enumerate(inserts):
            for K,y in enumerate(self.Adsorbent):                                                 
                if (K != usite) and (site.distance(y) < self.cutoff): # 2.4 is a cutoff 
                    remove.append(site)
                    remove_non_per.append(non_per[counter])
                    break
        new_inserts = [x for x in inserts if x not in remove]
        new_non_per = [x for x in non_per if x not in np.asarray(remove_non_per)]
        sorb_vecs,adsorbate_types = self.adsorbate_vectors()
        if len(sorb_vecs) > 1:                               
            for counter,site in enumerate(new_inserts): #enumerte is used so .remove can be used                                                     
                non_p_site = new_non_per[counter]                                                                                      
                pie = math.pi                                                       
                val = [0, pie/4, pie/2, pie*3/4, pie, -pie/4, -pie/2, -pie*3/4]#,-pie]               
                vec = non_p_site - self.Adsorbent[usite].coords                          
                pln_ap = self.make_pln(usite,vec)                              
                axis = []                                                           
                for x in val:                                                       
                    axis.append(self.arbitraty_axis(pln_ap, x, vec, usite))        
                dist = 0  
                rot_angles = np.asarray(np.linspace(0,360,9,dtype='int'))
                np.delete(rot_angles,-1) # dont need 360 degree rotations
                # count10 = 0 # used for vizualizing all rotations
                for x in axis:
                    atom_spots = []                                                      
                    theta_up = np.arcsin(0.5*(self.diatomic/self.bonding_dist))                           
                    theta_down = np.arcsin(-0.5*(self.diatomic/self.bonding_dist))                             
                    rotate_up = self.rotation_matrix(x, theta_up)                        
                    rotate_down = self.rotation_matrix(x, theta_down) 
                    atom_spots.append(np.dot(rotate_up,vec))
                    atom_spots.append(np.dot(rotate_down,vec)) 

                    period_sites = []
                    for count,new in enumerate(atom_spots):
                        atom = np.array(new, dtype=np.float64) + self.Adsorbent[usite].coords
                        ATOM = PeriodicSite(adsorbate_types[count],atom,self.Adsorbent.lattice,to_unit_cell=True,coords_are_cartesian=True)
                        period_sites.append(ATOM)
                    ##########################################################
                    # # Strickly for visualizing all rotation of adsorbate tested
                    # for g in period_sites:                                                  
                    #     self.Adsorbent.append(g.species_string,g.coords,coords_are_cartesian=True)
                    # pos = Poscar(self.Adsorbent)                                             
                    # pos.write_file('POSCAR_{}_Site_{}_{}'.format(self.Adsorbent[usite].species_string,spot,count10))
                    # count10 += 1
                    # for h in period_sites: # just usign x as a counter 
                    #     self.Adsorbent.remove(self.Adsorbent.sites[-1]) 
                    ##########################################################
                    tot = self.dist_from_lattice(period_sites)                    
                    if (tot > dist):                                                
                        best_sites = period_sites                                                                           
                        dist = tot
                ## check one more time to ensure adsrobate is not sitting too cloase to any oher atoms
                final_check = self.shortest_dist_from_lattice(best_sites) 
                if final_check > 0.5:
                    for x in best_sites:                                                  
                        self.Adsorbent.append(x.species_string,x.coords,coords_are_cartesian=True)
                        pos = Poscar(self.Adsorbent)                                             
                        pos.write_file('{}_{}_Site_{}'.format('POSCAR',self.Adsorbent[usite].species_string,spot))
                        if len(self.dope_list) > 0:
                            self.dope_me(usite,spot,'') # Create same site but doped      
                    for x in best_sites: # just usign x as a counter 
                        self.Adsorbent.remove(self.Adsorbent.sites[-1])                                                            
                else:
                    pass   
                spot = spot + 1                                               
            if len(self.dope_list) > 0:
	            self.dope_me(usite,start_spot,'_{}_no_O2'.format(int(spot)-1)) #creates doped usite with no O2

        return spot    
 
    def adsorbate_vectors(self):
        origin_site = self.Adsorbate[self.bonding_atom-1] # this is assiging the pymatgen PeriodicSite
        origin = origin_site.coords   # this is getting the cart coords for the site
        sorb_vecs = np.zeros((len(self.Adsorbate),3))  # intilize vector
        # loop through adsorbate and create a list of vectros from the origin site 
        # which will be the site that bonds to the absorbent
        adsorbate_types = []
        for count,x in enumerate(self.Adsorbate):
            sorb_vecs[count,:] = x.coords - origin
            adsorbate_types.append(x.species_string)
        return sorb_vecs , adsorbate_types

