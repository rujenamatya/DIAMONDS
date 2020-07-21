# Query code for Chebi Ontology
# Rujen Amatya, Yan Zhou Chen, Xinmeng Li


# install owlready 2 in your command shell: pip install owlready2

# on your python shell:

# import the owlready2 package
from owlready2 import *

# import fnmatch for wildcard query search
import fnmatch

# retrieve the ontology file (.owl) from your source
chebi_ont = get_ontology("http://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.owl").load()

# create the namespace for the ChEBI ontology. You need this to use ChEBI id
obo = chebi_ont.get_namespace("http://purl.obolibrary.org/obo/")

# function print_list
# parameter: a list of chebi classes
# returns: nothing
# function: prints the chebi id and name of the chebi class
def print_list(lst):  
    for x in range(len(lst)): 
        if (lst[x].name == "Thing"): 
            print ("none") 
        else: 
            print (lst[x], lst[x].label[0])
            
# function children_list
# parameter: a chebi class
# returns: list of the descendants of the specified levels of ancestors
def children_list(x):
  
    # get has_part and functional_parent relation 
    has_part_list = x.BFO_0000051
    has_func_par_list = x.has_functional_parent
    child_of_anc_with_smiles = []

    for x in range(len(has_part_list)): 
        if (has_part_list[x].smiles and has_part_list[x].monoisotopicmass): 
            child_of_anc_with_smiles.append(has_part_list[x])

    for x in range(len(has_func_par_list)): 
        if (has_func_par_list[x].smiles and has_func_par_list[x].monoisotopicmass): 
            child_of_anc_with_smiles.append(has_func_par_list[x])

    # get the parent list
    parent_list = ont.get_parents_of(x)
    
    # get the grandparent list
    gparent_list = []
    for x in range(len(parent_list)):
        if (parent_list[x].ancestors()):
            gparent_list = gparent_list + ont.get_parents_of(parent_list[x])
    gparent_list = list(set(gparent_list))
    
    # get the great grandparent list
    ggparent_list = []
    for x in range(len(gparent_list)):
        if (gparent_list[x].ancestors()):
            ggparent_list = ggparent_list + ont.get_parents_of(gparent_list[x])
    ggparent_list = list(set(ggparent_list))
    
    # get the great great grandparent list
    gggparent_list = []
    for x in range(len(ggparent_list)):
        if (ggparent_list[x].ancestors()):
            gggparent_list = gggparent_list + ont.get_parents_of(ggparent_list[x])
    gggparent_list = list(set(gggparent_list))
            
    # here change to the following:
    # gparent_list: for 3 levels
    # ggparent_list: for 4 levels
    # gggparent_list: for 5 levels
    child_of_anc = []
    for x in range(len(ggparent_list)):
        if (ggparent_list[x].ancestors()):
            child_of_anc = child_of_anc + list(ggparent_list[x].descendants())
    child_of_anc = list(set(child_of_anc))
    
    
    for x in range(len(child_of_anc)): 
        if (child_of_anc[x].smiles and child_of_anc[x].monoisotopicmass): 
            child_of_anc_with_smiles.append(child_of_anc[x])
    return child_of_anc_with_smiles

# function mass_query_list
# parameter: query mass and the children list obtained from children_list function
# returns: list of the chebi classes with the specified mass in the range -H and +3H
def mass_query_list(qm, child_of_anc_with_smiles):
    mass_child = []
    for x in range(len(child_of_anc_with_smiles)):
        mass_child.append(child_of_anc_with_smiles[x].monoisotopicmass[0])
        
    
    qmass_int = int(qm)
    qmass = str(qmass_int) + '.*'
    
    qmass1_int = qmass_int + 1
    qmass1 = str(qmass1_int) + '.*'
    
    qmass2_int = qmass_int + 2
    qmass2 = str(qmass2_int) + '.*'
    
    qmass3_int = qmass_int + 3
    qmass3 = str(qmass3_int) + '.*'

    qmass4_int = qmass_int + 4
    qmass4 = str(qmass4_int) + '.*'
    
    qmassn1_int = qmass_int - 1
    qmassn1 = str(qmassn1_int) + '.*'

    qmassn2_int = qmass_int - 2
    qmassn2 = str(qmassn2_int) + '.*'

    qmassnH_float = qm - 1.007825
    qmass3H_float = qm + 3 * 1.007825
    qmassnH = str(qmassnH_float)
    qmass3H = str(qmass3H_float)
    
    # query mass
    filtered = fnmatch.filter(mass_child, qmass)
    index_list = []
    for y in filtered:
        indexes = [i for i, x in enumerate(mass_child) if x == y]
        index_list = index_list + indexes
    index_list = list(set(index_list))
    query_index_list = []
    query_index_list = [int(i) for i in index_list]
    query_list = []
    for x in query_index_list:
        query_list.append(child_of_anc_with_smiles[x])

    # query mass+1   
    filtered1 = fnmatch.filter(mass_child, qmass1)
    index_list1 = []
    for y in filtered1:
        indexes = [i for i, x in enumerate(mass_child) if x == y]
        index_list1 = index_list1 + indexes
    index_list1 = list(set(index_list1))
    query_index_list1 = []
    query_index_list1 = [int(i) for i in index_list1]
    query_list1 = []
    for x in query_index_list1:
        query_list1.append(child_of_anc_with_smiles[x])
    
    # query mass+2 
    filtered2 = fnmatch.filter(mass_child, qmass2)
    index_list2 = []
    for y in filtered2:
        indexes = [i for i, x in enumerate(mass_child) if x == y]
        index_list2 = index_list2 + indexes
    index_list2 = list(set(index_list2))
    query_index_list2 = []
    query_index_list2 = [int(i) for i in index_list2]
    query_list2 = []
    for x in query_index_list2:
        query_list2.append(child_of_anc_with_smiles[x])
    
    # query mass+3  
    filtered3 = fnmatch.filter(mass_child, qmass3)
    index_list3 = []
    for y in filtered3:
        indexes = [i for i, x in enumerate(mass_child) if x == y]
        index_list3 = index_list3 + indexes
    index_list3 = list(set(index_list3))
    query_index_list3 = []
    query_index_list3 = [int(i) for i in index_list3]
    query_list3 = []
    for x in query_index_list3:
        if (not float(child_of_anc_with_smiles[x].monoisotopicmass[0]) >= qmass3H_float):
          query_list3.append(child_of_anc_with_smiles[x])

    # query mass+4  
    filtered4 = fnmatch.filter(mass_child, qmass4)
    index_list4 = []
    for y in filtered4:
        indexes = [i for i, x in enumerate(mass_child) if x == y]
        index_list4 = index_list4 + indexes
    index_list4 = list(set(index_list4))
    query_index_list4 = []
    query_index_list4 = [int(i) for i in index_list4]
    query_list4 = []
    for x in query_index_list4:
        if (not float(child_of_anc_with_smiles[x].monoisotopicmass[0]) >= qmass3H_float):
          query_list4.append(child_of_anc_with_smiles[x])

    # query mass-1  
    filteredn1 = fnmatch.filter(mass_child, qmassn1)
    index_listn1 = []
    for y in filteredn1:
        indexes = [i for i, x in enumerate(mass_child) if x == y]
        index_listn1 = index_listn1 + indexes
    index_listn1 = list(set(index_listn1))
    query_index_listn1 = []
    query_index_listn1 = [int(i) for i in index_listn1]
    query_listn1 = []
    for x in query_index_listn1:
        if (not float(child_of_anc_with_smiles[x].monoisotopicmass[0]) <= qmassnH_float):
          query_listn1.append(child_of_anc_with_smiles[x])

    # query mass-2  
    filteredn2 = fnmatch.filter(mass_child, qmassn2)
    index_listn2 = []
    for y in filteredn2:
        indexes = [i for i, x in enumerate(mass_child) if x == y]
        index_listn2 = index_listn2 + indexes
    index_listn2 = list(set(index_listn2))
    query_index_listn2 = []
    query_index_listn2 = [int(i) for i in index_listn2]
    query_listn2 = []
    for x in query_index_listn2:
        if (not float(child_of_anc_with_smiles[x].monoisotopicmass[0]) <= qmassnH_float):
          query_listn2.append(child_of_anc_with_smiles[x])

    all_query = []
    all_query = all_query + query_list + query_list1 + query_list2 + query_list3 + query_list4 + query_listn1 +query_listn2
    all_query = list(set(all_query))
    return all_query

