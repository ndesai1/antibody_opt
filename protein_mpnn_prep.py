'''Functions to extract non-CDR residues that will be kept fixed
while CDR regions are varied using protein MPNN'''
import os
from get_ab_seq import Ab_Ag_Seq_Extractor


class Protein_MPNN_prep:
    def __init__(self, pdb_id, prot_mpnn_dir):
        self.pdb_id = pdb_id
        self.prot_mpnn_dir = prot_mpnn_dir

    def prot_prep(self):
        '''
        extract residue information, residue numbers to be kept fixed as
        well as CDR sequences that will be varied by protein MPNN
        '''
        if not os.path.exists(self.prot_mpnn_dir):
            os.makedirs(self.prot_mpnn_dir)

        extractor = Ab_Ag_Seq_Extractor(self.pdb_id)
        antibody_info = extractor.get_cdr_sequences()
        light_fixed_ranges = antibody_info['light_fixed_ranges']
        heavy_fixed_ranges = antibody_info['heavy_fixed_ranges']
        fixed_list_heavy = list(heavy_fixed_ranges[0]) + list(heavy_fixed_ranges[1]) + \
                            list(heavy_fixed_ranges[2])
        fixed_list_light = list(light_fixed_ranges[0]) + list(light_fixed_ranges[1]) + \
                            list(light_fixed_ranges[2])
        path_for_parsed_chains ="{}/parsed_pdbs.jsonl".format(self.prot_mpnn_dir)
        path_for_assigned_chains ="{}/assigned_pdbs.jsonl".format(self.prot_mpnn_dir)
        path_for_fixed_positions = self.prot_mpnn_dir
        chains_to_design = "H L"
        heavy_positions = ' '.join([str(element) for element in fixed_list_heavy])
        light_positions = ' '.join([str(element) for element in fixed_list_light])
        design_only_positions = "{}, {}".format(heavy_positions, light_positions)







