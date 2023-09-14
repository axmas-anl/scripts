import atomman as am
from diffpy.structure.parsers import getParser


def read_cif(cif_file):
    p = getParser('cif')
    cif = p.parseFile(cif_file)

    elements, atype = np.unique(cif.element, return_inverse=True)
    atype += 1
    pos = cif.xyz

    atoms = am.Atoms(atype=atype, pos=pos)
   
    box = am.Box(a=cif.lattice.a, b=cif.lattice.b, c=cif.lattice.c,
                 alpha=cif.lattice.alpha, beta=cif.lattice.beta, gamma=cif.lattice.gamma)

    return am.System(atoms=atoms, box=box, scale=True, symbols=elements)
    

def find_vector(s, v, scale=1.0):
    abc = np.matrix((s.box.avect, s.box.bvect, s.box.cvect)).T
    xyz = np.array(abc * np.matrix(v).T).T[0] / scale
    neighbors = am.NeighborList(system=s, cutoff=np.linalg.norm(xyz)*1.1)
    def iv(vect):
        vv = scale * np.array((vect.T*np.linalg.inv(abc.T)))[0]
        return '(%g, %g, %g)' % (vv[0], vv[1], vv[2])
    ivs = []
    for (i,n) in enumerate(neighbors):
        a1 = s.symbols[s.atoms[i].atype.item()-1]
        for j in n:
            a2 = s.symbols[s.atoms[j].atype.item()-1]
            dv = s.dvect(i,j)
            if ((np.linalg.norm(dv-xyz) < 0.5 or np.linalg.norm(xyz-dv) < 0.5)):
                print(a1, iv(s.atoms[i].pos[0]), a2, iv(s.atoms[j].pos[0]), iv(dv))
                ivs.append(list(eval(iv(dv))))
    return ivs

def no_really_find_all_the_vectors(s):
    abc = np.matrix((s.box.avect, s.box.bvect, s.box.cvect)).T
    neighbors = am.NeighborList(system=s, cutoff=np.linalg.norm(abc)*1.1)
    # can't iterate s.atoms directly
    atom_labels = [s.symbols[s.atoms[i].atype.item()-1] for i in range(len(s.atoms))]
    list_of_numbers = list(range(len(s.atoms)))
    xx,yy = np.meshgrid(list_of_numbers, list_of_numbers)
    vectors = s.dvect(xx.ravel(), yy.ravel())

    labels = np.empty((len(s.atoms),len(s.atoms)), dtype="object")

    for i in range(len(s.atoms)):
        for j in range(len(s.atoms)):
            labels[i,j] = "{}-{}".format(atom_labels[i], atom_labels[j])
            

    return vectors.reshape((len(s.atoms),len(s.atoms),3)), labels

    

def find_all_the_vectors(s, plane, offset, scale=1.0, ax=None):
    # plane should be "hk0" or similar even if 0 is not zero
    # and z (i.e. the offset dir)
    vectors, labels = no_really_find_all_the_vectors(s)

    # what is the first vertical axis?
    offset_idx = plane.find("0")

    # connvert offset to Angstroms HACK
    unit_cell = [s.box.a, s.box.b, s.box.c]

    vertical_angstroms = offset * unit_cell[offset_idx]

    the_plane_of_interest = vectors[:,:,offset_idx]

    distance_from_offset = np.abs(the_plane_of_interest - vertical_angstroms)

    whether_to_plot = distance_from_offset < 0.5
    
    if plane == "xy0":
        plane_x = vectors[:,:,0]/s.box.a
        plane_y = vectors[:,:,1]/s.box.b
    elif plane == "x0z":
        plane_x = vectors[:,:,0]/s.box.a
        plane_y = vectors[:,:,2]/s.box.c
    elif plane == "0yz":
        plane_x = vectors[:,:,1]/s.box.b
        plane_y = vectors[:,:,2]/s.box.c
    
    px = plane_x[whether_to_plot].ravel()
    py = plane_y[whether_to_plot].ravel()
    ll = labels[whether_to_plot].ravel()
    
    
    if not ax:
        ax = plt.gca()

    annotations = []
    for x,y,l in zip(px,py,ll):
        if x<1e-6 and y<1e-6:
            continue
        _ox = (np.random.random()-0.5 ) /5
        _oy = (np.random.random()-0.5 ) /5
        annotations.append(ax.annotate(l, (x,y), xytext=(x+_ox,y+_oy), arrowprops={"width":0.5, "headwidth":0.5}))

    return annotations


class AtomAnnotations():
    def __init__(self, atomman_structure, axis_handle):
        self.s = atomman_structure
        vectors, labels = no_really_find_all_the_vectors(atomman_structure)
        self.vectors = vectors
        self.labels = labels
        self.annotations = []
        self.axis_handle = axis_handle
    def clear(self):
        [aa.remove() for aa in self.annotations]
    def plot(self, plane, offset, threshold=0.5):
        # what is the first vertical axis?
        offset_idx = plane.find("0")

        # connvert offset to Angstroms HACK
        unit_cell = [self.s.box.a, self.s.box.b, self.s.box.c]

        vertical_angstroms = offset * unit_cell[offset_idx]

        the_plane_of_interest = self.vectors[:,:,offset_idx]

        distance_from_offset = np.abs(the_plane_of_interest - vertical_angstroms)

        whether_to_plot = distance_from_offset < threshold
    
        if plane == "xy0":
            plane_x = self.vectors[:,:,0]/s.box.a
            plane_y = self.vectors[:,:,1]/s.box.b
        elif plane == "x0z":
            plane_x = self.vectors[:,:,0]/s.box.a
            plane_y = self.vectors[:,:,2]/s.box.c
        elif plane == "0yz":
            plane_x = self.vectors[:,:,1]/s.box.b
            plane_y = self.vectors[:,:,2]/s.box.c
    
        px = plane_x[whether_to_plot].ravel()
        py = plane_y[whether_to_plot].ravel()
        ll = self.labels[whether_to_plot].ravel()
    
        ax = self.axis_handle

        for x,y,l in zip(px,py,ll):
            if x<1e-6 and y<1e-6:
                continue
            _ox = (np.random.random()-0.5 ) /5
            _oy = (np.random.random()-0.5 ) /5
            self.annotations.append(ax.annotate(l, (x,y), xytext=(x+_ox,y+_oy), arrowprops={"width":0.5, "headwidth":0.5}))
