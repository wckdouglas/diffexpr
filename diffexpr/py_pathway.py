from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
import numpy as np

def gene_is_enriched(enriched_genes, possible_names):
    '''
    check if any gene name in the box matched the enriched genes
    return the first matched gene name
    '''
    gene_names = possible_names.replace(' ','').split(',')
    for gn in gene_names:
        if gn in enriched_genes:
            return gn


class pathwayview:
    '''
    Like pathwayview in R, plotting a pathway and labeling the changed genes

    Usage:
        args:
            enriched_genes: a list of enriched genes that are in the pathway
            pathway_id: KEGG pathway ID (e.g. hsa04110)
            figurename: Must be pdf file
        
        returns:
            pathway: a biopython KEGG KGML object

    Example:

        enriched_genes = 'GADD45A,PLK1,TTK,CDC6,CDC25C,CDC25A'.split(',')
        pv = pathwayview()
        pathway = pv.plot_pathway(enriched_genes = enriched_genes, 
                        pathway_id = 'hsa04110',
                        figurename = 'pathway.pdf')


        # changing labeling colors:
        pv.enriched_color = '#00000' #must be hex codes


    '''


    def __init__(self, fontsize = 12):
        self.enriched_box_color = '#FFB600' #oragne
        self.non_enriched_box_color = '#FFFFFF' #white
        self.enriched_text_color = '#000000' #black
        self.non_enriched_text_color = '#000000'
        self.fontsize = fontsize

    
    def set_color(self, 
            enriched_box_color = '#FFB600',
            non_enriched_box_color = '#FFFFFF', 
            enriched_text_color = '#000000',
            non_enriched_text_color = '#000000'):
        '''
        Color control, inputs must be HEX color code
        '''
        self.enriched_box_color = enriched_box_color #oragne
        self.non_enriched_box_color = non_enriched_box_color #white
        self.enriched_text_color = enriched_text_color #black
        self.non_enriched_text_color = non_enriched_text_color


    def plot_pathway(self, 
                    enriched_genes, 
                    pathway_id = 'hsa05322', 
                    figurename = None):

        # config figure name
        if not figurename:
            figurename = '%s.pdf' %pathway_id
        assert(figurename.endswith('.pdf'))
        
        # fetch pathway
        pathway = KGML_parser.read(kegg_get(pathway_id, "kgml"))
        
        # change color for pathway elements
        for entry in pathway.entries.values():
            possible_gene_names = entry.graphics[0].name
            matched_name = gene_is_enriched(enriched_genes, possible_gene_names)
            if matched_name:
                entry.graphics[0].bgcolor = self.enriched_box_color #set box color
                entry.graphics[0].fgcolor = self.enriched_text_color # set text color
                entry.graphics[0].name = matched_name
            else:
                entry.graphics[0].bgcolor = self.non_enriched_box_color
                entry.graphics[0].fgcolor = self.non_enriched_text_color
                entry.graphics[0].name = entry.graphics[0].name.split(',')[0]
        
        canvas = KGMLCanvas(pathway, import_imagemap=True, fontsize=self.fontsize)
        canvas.draw(figurename)
        print('Drawn: ', figurename)
        return pathway