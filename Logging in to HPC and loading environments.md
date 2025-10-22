## Log in to HPC
Log in to HPC with [OpenOnDemand GUI](https://docs.ycrc.yale.edu/clusters-at-yale/access/ood/) or via SSH on your computer's terminal. 

To log in with SSH you need to create a specific key do this with the following:
```
ssh-keygen
cat ~/YaleSSHkey.pub #I named my key "YaleSSHkey" but you can name it whatever
```
Upload the key to Yale: https://sshkeys.ycrc.yale.edu/ 

Check that you can SSH into the Yale HPC and have the proper settings: 
```
ssh -i ~/YaleSSHkey {your ID}@mccleary.ycrc.yale.edu
nano ~/.ssh/config: Host mccleary
    HostName mccleary.ycrc.yale.edu
    User {your ID}
    IdentityFile ~/YaleSSHkey
chmod 600 ~/.ssh/config
chmod 600 ~/YaleSSHkey
chmod 700 ~/.ssh
```
Now you can use: `ssh -i ~/YaleSSHkey {your ID}@mccleary.ycrc.yale.edu`

## Load conda environment
To work in the HPC you need to activate an environment with packages will we use for RNA seq qnalysis.
Conda environment: rnaseq_tools (Packages: sra-tools, entrez-direct, fastqc, multiqc, samtools, etc...)

`module load miniconda`
`conda activate rnaseq_tools`

This was installed via: 

`conda install -c bioconda {package}` `conda create -n{name} -c bioconda {packages separated by spaces}`

You can check what packages are you the environment with `conda list`.

With the environment loaded you can install things needed in it later with the `conda install -c bioconda {package}` code

## Set up/go to working directory
`mkdir {path to your directory here, make an RNA seq folder and then sub folders for each project}`
`cd {path to your directory}`
