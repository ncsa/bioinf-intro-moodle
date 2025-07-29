# Define install directory

mkdir -p ~/tools
cd ~/tools

# Fastp
module load anaconda3
conda create -n my.anaconda python
source activate my.anaconda

conda install -c bioconda fastp

# bwa
conda install -c bioconda bwa

# samtools
wget https://sourceforge.net/projects/samtools/files/samtools/1.22.1/samtools-1.22.1.tar.bz2 -O samtools.tar.bz2
tar -xjvf samtools.tar.bz2
cd samtools-1.22.1
make
make prefix = $HOME/tools/samtools/bin install
echo PATH=$PATH:$HOME/tools/samtools/bin/bin >> ~/bash_profile
source ~/.bash_profile
cd ..

# bcftools
wget https://sourceforge.net/projects/samtools/files/samtools/1.22/bcftools-1.22.tar.bz2 -O bcftools.tar.bz2
tar -xjvf bcftools.tar.bz2
cd bcftools-1.22
make
make prefix =$HOME/tools/bcftools/bin install
echo PATH=$PATH:$HOME/tools/bcftools/bin/bin >> ~/.bash_profile
source ~/.bash_profile
cd ..

# GATK
wget https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip
unzip gatk-4.5.0.0.zip
echo PATH=$PATH:$HOME/tools/gatk-4.5.0.0/ >> ~/.bash_profile
source ~/.bash_profile
cd ..
