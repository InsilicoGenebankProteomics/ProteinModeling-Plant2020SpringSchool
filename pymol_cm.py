import numpy as np
import matplotlib.pyplot as plt
from pymol import cmd, stored

def plot_contact_map( arg1, arg2, max_dist ):

    print ("Calculating contact map.")
    print ("You passed in %s and %s" % (arg1, arg2))
    xyz1 = cmd.get_coords(arg1+' and name CA', 1)
    xyz2 = cmd.get_coords(arg2+' and name CA', 1)
    stored.ls1=[]
    cmd.iterate(arg1+' and name CA', 'stored.ls1.append(int(resi))')
    stored.ls2=[]
    cmd.iterate(arg2+' and name CA', 'stored.ls2.append(int(resi))')

    cm = []
    for i in range(len(xyz1)):
        cm.append([])
        for j in range(len(xyz2)):
            diff = xyz2[j] - xyz1[i]
            dist = np.sqrt(np.sum(diff * diff))
            cm[-1].append(dist)
         
    plt.imshow(cm, cmap="Greens_r", vmin=0, vmax=max_dist,extent=[stored.ls1[0],stored.ls1[-1],stored.ls2[0],stored.ls2[-1]])
    plt.colorbar(label="Distance (Ångströms)")
    plt.savefig('contact_map.png',dpi=300)
    plt.close()
