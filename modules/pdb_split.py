"""Modules for splitting up PDB file containing antibody structure
and extracting the heavy and light chains for the antibody"""

import os
from Bio import PDB

from pathlib import Path

class AbPDBSplit:
    def __init__(self,
                 pdb_id:str,
                 pdb_path: str,
                 out_dir=None,
                 ):
        """Class to parse antibody-antigen structures from
        PDB files, split them based on chain identifiers, and output
        split files to output directory

        Args:
            pdb_id (str): PDB ID to be used for the output antibody and antigen files
            pdb_path (str): path to the PDB file a known antibody-antigen complex
        """
        self.pdb_id = pdb_id
        self.pdb_path = pdb_path
        self.ab_parser = PDB.PDBParser()
        self.ab_writer = PDB.PDBIO()
        if out_dir is None:
            out_dir = os.path.join(os.getcwd(), "ab_ag_files")
        self.out_dir = out_dir

    def make_pdb(self,
                 suffix: str,
                 chain_list: list,
                 overwrite=False,
                 ):
        """ Takes in pdb structure and list of chains, outputs new
        PDB file only containing specified chains

        Args:
            suffix (str): suffix to be added to pdb output file
            chain_list (list): list of chains to be extracted to output file
            overwrite (bool): should output file be overwritten if it exists?
            struct: PDB structure object

        """
        #convert list of chains to upper case
        chain_list = [chain.upper() for chain in chain_list]

        # Input/output files
        out_name = "{}_{}.pdb".format(self.pdb_id, suffix)
        out_path = os.path.join(self.out_dir, out_name)

        # if overwrite is false, skip PDB chain file generation
        if (not overwrite) and (os.path.isfile(out_path)):
            print("output file {} already exists".format(out_path))
            return out_path
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)


        #parse input PDB structure, output file only with chains in chain_list
        struct = self.ab_parser.get_structure(self.pdb_id, self.pdb_path)
        self.ab_writer.set_structure(struct)
        self.ab_writer.save(out_path, select=SelectChains(chain_list))

        return out_path


class SelectChains(PDB.Select):
    def __init__(self,
                 chain_list
                 ):
        """Get PDB.Select list of PDB chains"""
        self.chain_list = chain_list

    def accept_chain(self, chain):
        selected_chains = (chain.get_id() in self.chain_list)
        return selected_chains