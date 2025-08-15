Repo for the extraction of antibody sequence and structure information for optimization of CDR regions with Protein MPNN

To install the repo:

1. Clone the git repo:
> git clone https://github.com/ndesai1/antibody_opt.git

2. Create conda environment using environment.yml file:
   
> cd antibody_opt

> conda env create --name antibody_opt --file environment.yml

> conda activate antibody_opt

3. Anarci and ProteinMPNN have to be installed to use antibody_opt:
   
> git clone https://github.com/oxpig/ANARCI.git

> cd ANARCI

> python setup.py install


> git clone https://github.com/dauparas/ProteinMPNN.git

> conda install pytorch torchvision torchaudio cudatoolkit=11.3 -c pytorch

Once these steps are completed, antibody_opt is ready to be used
