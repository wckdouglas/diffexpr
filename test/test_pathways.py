from diffexpr.py_pathway import pathwayview
from Bio.KEGG.KGML.KGML_pathway import Pathway
import os

def test_plotpathway():
    enriched_genes = 'GADD45A,PLK1,TTK,CDC6,CDC25C,CDC25A'.split(',')
    pv = pathwayview()
    pathway = pv.plot_pathway(enriched_genes = enriched_genes, 
                    pathway_id = 'hsa04110',
                    figurename = 'pathway.pdf')
    assert os.path.isfile('pathway.pdf') 
    assert ( isinstance(pathway, Pathway) )
