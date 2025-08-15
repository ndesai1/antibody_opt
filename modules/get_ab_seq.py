"""Module for extracting the sequence and structures for the
antibody and antigen, and getting the CDR regions
from the heavy and light chains of the antibodies"""

import os

from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.SeqUtils import seq1
from Bio import SeqIO

from pdb_split import AbPDBSplit
from modules.ab_seq_annotation import Annotate


class Ab_Ag_Seq_Extractor:
    def __init__(self,
                 pdb_id,
                 ab_chains=("H", "L"),
                 ):
        """ Class for loading in an antibody-antigen PDB file, split up antibody
        and antigen structures, and get sequence for the heavy and light chains and
        CDR regions

        Args:
            pdb_file (str): path to the PDB file a known antibody-antigen complex
            chains (list): heavy and light chain identifiers
        """
        self.pdb_id = pdb_id
        self.input_dir = os.getcwd()
        self.pdb_file = '{}/{}.pdb1'.format(self.input_dir, self.pdb_id)
        self.ab_chains = ab_chains

    def get_ag_chains(self)->list[str]:
        """
        get the list of chain_id(s) for antigen chains

        Returns:
            ag_chains (list): list of antibody chain IDs

        """
        ag_chains = []
        ab_parser = PDBParser()
        io = PDBIO()
        structure = ab_parser.get_structure(self.pdb_id, self.pdb_file)
        pdb_chains = structure.get_chains()
        for chain in pdb_chains:
            if chain.get_id() not in ['H', 'L']:
                ag_chains.append(chain.get_id())

        return ag_chains

    def split_ab_ag(self)->dict:
        """ split a PDB file containing antibody-antigen structures into
        two files: one containing the antibody chains and one containing the antigen
        chain(s)
        """
        ag_chains = self.get_ag_chains()
        # create file with antibody structure
        pdb_split = AbPDBSplit(self.pdb_id, self.pdb_file)
        ab_file_path = pdb_split.make_pdb(suffix='ab', chain_list=['H', 'L'])
        ag_file_path = pdb_split.make_pdb(suffix='ag', chain_list=ag_chains)
        file_list = {'antibody_file':ab_file_path, 'antigen_file': ag_file_path}
        return file_list

    def get_cdr_sequences(self)->dict:
        """
        Extract CDR loop regions for heavy and light antibody chains

        Returns:
            CDR_list (dict): dictionary containing CDR regions of antibody structure
        """
        file_list = self.split_ab_ag()
        ab_parser = PDBParser()
        structure = ab_parser.get_structure(self.pdb_id, self.pdb_file)
        chains = {chain.id: seq1(''.join(residue.resname for residue in chain)) for chain in structure.get_chains()}

        heavy_chain = chains['H']
        light_chain = chains['L']

        scheme = "imgt"
        heavy_seq_annote = Annotate(heavy_chain, scheme)
        heavy_annotation = heavy_seq_annote.retrieve()
        heavy_annotation = heavy_annotation[0]
        light_seq_annote = Annotate(light_chain, scheme)
        light_annotation = light_seq_annote.retrieve()
        light_annotation = light_annotation[0]

        self.CDR_seqs = {'H-CDR1':heavy_annotation['H-CDR1'],
               'H-CDR2':heavy_annotation['H-CDR2'],
               'H-CDR3':heavy_annotation['H-CDR3'],
               'L-CDR1':light_annotation['L-CDR1'],
               'L-CDR2':light_annotation['L-CDR2'],
               'L-CDR3':light_annotation['L-CDR3']}

        heavy_region_FR1 = len(heavy_annotation['H-FR1'])
        heavy_region_CDR1 = heavy_region_FR1 + len(heavy_annotation['H-CDR1'])
        heavy_region_FR2 = heavy_region_CDR1 + len(heavy_annotation['H-FR2'])
        heavy_region_CDR2 = heavy_region_FR2 + len(heavy_annotation['H-CDR2'])
        heavy_region_FR3 = heavy_region_CDR2 + len(heavy_annotation['H-FR3'])
        heavy_region_CDR3 = heavy_region_FR3 + len(heavy_annotation['H-CDR3'])


        self.heavy_fixed_ranges = [range(0, heavy_region_FR1),
                             range(heavy_region_CDR1+1, heavy_region_FR2),
                             range(heavy_region_FR3+1, heavy_region_CDR3),
                             ]

        light_region_FR1 = len(light_annotation['L-FR1'])
        light_region_CDR1 = light_region_FR1 + len(light_annotation['L-CDR1'])
        light_region_FR2 = light_region_CDR1 + len(light_annotation['L-FR2'])
        light_region_CDR2 = light_region_FR2 + len(light_annotation['L-CDR2'])
        light_region_FR3 = light_region_CDR2 + len(light_annotation['L-FR3'])
        light_region_CDR3 = light_region_FR3 + len(light_annotation['L-CDR3'])

        self.light_fixed_ranges = [range(0, light_region_FR1),
                                   range(light_region_CDR1 + 1, light_region_FR2),
                                   range(light_region_FR3 + 1, light_region_CDR3),
                                   ]

        self.anchor_slices = {
                "H1": ('H', (25, 40)),
                "H2": ('H', (54, 67)),
                "H3": ('H', (103, 119)),
                "L1": ('L', (25, 40)),
                "L2": ('L', (54, 67)),
                "L3": ('L', (103, 119))}

        antibody_info: dict = {'CDR_seqs' : self.CDR_seqs,
                               'heavy_fixed_ranges': self.heavy_fixed_ranges,
                               'light_fixed_ranges': self.light_fixed_ranges}

        return antibody_info




