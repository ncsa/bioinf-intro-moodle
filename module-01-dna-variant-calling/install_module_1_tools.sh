# Define install directory
mkdir -p ~/tools/bin
cd ~/tools

# Add to PATH
echo 'export PATH=$HOME/tools/bin:$PATH' >> ~/.bashrc
source ~/.bashrc

# 1. Install BWA
wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2
tar -xvjf bwa-0.7.17.tar.bz2
cd bwa-0.7.17 && make
cp bwa ../bin/
cd ..

# 2. Install Samtools
wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2
tar -xvjf samtools-1.18.tar.bz2
cd samtools-1.18
./configure --prefix=$HOME/tools && make && make install
cd ..

# 3. Install bcftools
wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2
tar -xvjf bcftools-1.20.tar.bz2
cd bcftools-1.20
make && make prefix=$HOME/tools install
cd ..

# 4. Download GATK (prebuilt Java jar)
mkdir gatk && cd gatk
wget https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip
unzip gatk-4.5.0.0.zip
cp gatk-4.5.0.0/gatk $HOME/tools/bin/
cd ..

# 5. Install fastp
wget http://opengene.org/fastp/fastp
chmod +x fastp
mv fastp ~/tools/bin/

