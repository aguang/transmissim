import numpy as np
from numpy.random import gamma
import math
import ete2

"""
Uses a birth-death process to generate gene trees... only allow 1 birth rate, death rate across branches
"""

def generate(E, S, lam, mu, genetree):
    """
    Generates birth and death of a gene tree lineage along a single branch
    in the species tree drawn from an exponential distribution.
    """
    tbegin = 0
    tend = E.dist
    k = S.copy_num
    C = genetree.search_nodes(name=E.name)

    while tbegin < tend and C and k != 0:
        W=np.random.exponential(scale=1/(k*(lam + mu)), size=1)[0]
        if W<tend-tbegin:
            tbegin = tbegin+W
            X=np.random.uniform(0.0,1.0)
            if X<(lam*k/((lam+mu)*k)):
                copy_num = np.random.randint(k);
                parent = C[copy_num].up
                original = C[copy_num].detach()
                copy = original.copy()
                new_child = parent.add_child(name='B', dist=W)
                new_child.add_child(original, dist=tend-tbegin)
                new_child.add_child(copy, dist=tend-tbegin)
                k = k+1
                C.append(copy)
            else:
                copy_num = np.random.randint(k)
                parent = C[copy_num].up
                C[copy_num].detach()
                parent.add_child(name='D', dist=tbegin)
                k = k - 1
                del C[copy_num]
        else: break
        """
        print "tbegin is: ", tbegin
        print "tend is: ", tend
        tbirth=np.random.exponential(scale=1/(k*lam), size=1)[0]
        tdeath=np.random.exponential(scale=1/(k*mu), size=1)[0]
        if tbirth < tdeath and tbirth < tend and tbegin < tbirth:
            copy_num = np.random.randint(k);
            parent = C[copy_num].up
            original = C[copy_num].detach()
            new_child = parent.add_child(name='B', dist=tbirth)
            new_child.add_child(original, dist=tend-tbirth)
            copy = original.copy()
            new_child.add_child(copy,dist=tend-tbirth)
            k = k + 1
            C.append(copy)
            tbegin = tbirth
        elif tdeath < tbirth and tdeath < tend and tdeath > tbegin:
            copy_num = np.random.randint(k)
            parent = C[copy_num].up
            C[copy_num].detach()
            parent.add_child(name='D', dist=tdeath)
            k = k - 1
            tbegin = tdeath
            del C[copy_num]
        else:
            break
        """
    E.add_feature('copy_num', k)

def filter_tree(tree):
    """
    removes deletions from the generated gene tree
    """
    tree.resolve_polytomy(recursive=True)
    D = tree.search_nodes(name="D")
    for death in D:
        P = death.up
        if P.is_root():
            death.detach()
            children = P.children
            if children:
                child = P.children[0] #assumes bifurcating tree
                child.dist = P.dist + child.dist
                P.delete()
            else:
                return 0
        else:
            G = P.up
            G.dist=G.dist+P.dist
            death.detach()
            P.delete()
    return tree    
    

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description="generates gene trees with incorporated missing data rate")

    ds = ' [%(default)s]'
    parser.add_argument('-sp', '--species_tree', help='shared species tree')
    parser.add_argument('-n', '--num_genefams')
    parser.add_argument('--out_dir', help='out directory')

    opts = parser.parse_args()
    species_tree_file = opts.species_tree
    out_dir = opts.out_dir

    sp = ete2.Tree(species_tree_file)
    num_genefams = int(opts.num_genefams)
    genetrees = []

    change_rate = gamma(.625, .03, num_genefams) # 30 mya
#    change_rate = gamma(4, .0004, num_genefams) 1mya
    for i in range(num_genefams):
        # attribute set up for simulation
        lam = change_rate[i]/3
        mu = change_rate[i] - lam
        """
        print "tree number: ", i
        print "change rate: ", change_rate[i]
        print "birth rate: ", lam
        print "death rate: ", mu
        print "----------"
        """
        genetree = sp.copy()
        descendants = genetree.get_descendants('preorder')
        internal_order=0
        for node in descendants:
            node.add_feature('copy_num', 1)
            if node.name=='NoName': node.name='I'+str(internal_order)
            internal_order=internal_order+1
        root = genetree.get_tree_root()
        root.add_feature('copy_num', 1)

        # actual simulation
        for node in descendants:
            E = node
            S = node.up
            generate(E, S, lam, mu, genetree)
            filter_tree(genetree)
        genetrees.append(genetree)

        # gene tree $i will be written into out_dir/species_tree_file${i}
        out = out_dir + species_tree_file + str(i)
        genetree.write(format=1, outfile=out)

        
