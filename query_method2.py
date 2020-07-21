#install owlready 2 in your command shell: pip install owlready2

#on your python shell:

#import the owlready2 package
from owlready2 import *

import numpy as np

import fnmatch

#retrieve the ontology file (.owl) from your source
chebi_ont = get_ontology("http://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.owl").load()

#create the namespace for the ChEBI ontology. You need this to use ChEBI id
obo = chebi_ont.get_namespace("http://purl.obolibrary.org/obo/")

#list of classes (these are the molecules)
#all_classes = list(chebi_ont.classes())

#list of properties
#all_properties = list(chebi_ont.properties())

#simple query with mass
#query_list = chebi_ont.search(mass = "155.15468")
#for x in range(len(query_list)):
#   print (query_list[x])
    

#testing
# =============================================================================
# for x in range(len(anc_of_caps)): 
#     if (anc_of_caps[x].name == "Thing"): 
#         print ("none") 
#     else: 
#         print (anc_of_caps[x].label[0])
#         
# for x in range(len(anc_of_caps)): 
#     if (anc_of_caps[x].name == "Thing"): 
#         print ("owl Thing") 
#     elif (not anc_of_caps[x].smiles) :
#         print ("no smiles")
#     else: 
#         print (anc_of_caps[x].smiles[0])
# =============================================================================
def print_list(lst):  
    for x in range(len(lst)): 
        if (lst[x].name == "Thing"): 
            print ("none") 
        else: 
            print (lst[x], lst[x].label[0])
            

def children_list(x):
  
    parent_list = ont.get_parents_of(x)
    
    gparent_list = []
    for x in range(len(parent_list)):
        if (parent_list[x].ancestors()):
            gparent_list = np.concatenate((gparent_list, ont.get_parents_of(parent_list[x])))
    gparent_list = list(set(gparent_list))
    
    ggparent_list = []
    for x in range(len(gparent_list)):
        if (gparent_list[x].ancestors()):
            ggparent_list = np.concatenate((ggparent_list, ont.get_parents_of(gparent_list[x])))
    ggparent_list = list(set(ggparent_list))
    
    gggparent_list = []
    for x in range(len(ggparent_list)):
        if (ggparent_list[x].ancestors()):
            gggparent_list = np.concatenate((gggparent_list, ont.get_parents_of(ggparent_list[x])))
    gggparent_list = list(set(gggparent_list))
            
    child_of_anc = []
    for x in range(len(ggparent_list)):
        if (ggparent_list[x].ancestors()):
            child_of_anc = np.concatenate((child_of_anc, list(ggparent_list[x].descendants())))
    child_of_anc = list(set(child_of_anc))
    
    child_of_anc_with_smiles = []
    for x in range(len(child_of_anc)): 
        if (child_of_anc[x].smiles and child_of_anc[x].monoisotopicmass): 
            child_of_anc_with_smiles.append(child_of_anc[x])
    return child_of_anc_with_smiles


def mass_query_list(x, child_of_anc_with_smiles):
    mass_child = []
    for x in range(len(child_of_anc_with_smiles)):
        mass_child.append(child_of_anc_with_smiles[x].monoisotopicmass[0])
        
    
    qmass_int = x
    qmass = str(qmass_int) + '.*'
    
    qmass1_int = qmass_int + 1
    qmass1 = str(qmass1_int) + '.*'
    
    qmass2_int = qmass_int + 2
    qmass2 = str(qmass2_int) + '.*'
    
    qmass3_int = qmass_int + 3
    qmass3 = str(qmass3_int) + '.*'
    
    qmassn1_int = qmass_int - 1
    qmassn1 = str(qmassn1_int) + '.*'
    
    # query mass
    filtered = fnmatch.filter(mass_child, qmass)
    index_list = []
    for y in filtered:
        indexes = [i for i, x in enumerate(mass_child) if x == y]
        index_list = np.concatenate((index_list, indexes))
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
        index_list1 = np.concatenate((index_list1, indexes))
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
        index_list2 = np.concatenate((index_list2, indexes))
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
        index_list3 = np.concatenate((index_list3, indexes))
    index_list3 = list(set(index_list3))
    query_index_list3 = []
    query_index_list3 = [int(i) for i in index_list3]
    query_list3 = []
    for x in query_index_list3:
        query_list3.append(child_of_anc_with_smiles[x])
        
    # query mass-1  
    filteredn1 = fnmatch.filter(mass_child, qmassn1)
    index_listn1 = []
    for y in filteredn1:
        indexes = [i for i, x in enumerate(mass_child) if x == y]
        index_listn1 = np.concatenate((index_listn1, indexes))
    index_listn1 = list(set(index_listn1))
    query_index_listn1 = []
    query_index_listn1 = [int(i) for i in index_listn1]
    query_listn1 = []
    for x in query_index_listn1:
        query_listn1.append(child_of_anc_with_smiles[x])
     
    all_query = []
    all_query = np.concatenate((all_query, query_list, query_list1, query_list2, query_list3, query_listn1))
    all_query = list(set(all_query))
    return all_query


print(len(all_query))
#print(all_query.index(obo.CHEBI_64723))

