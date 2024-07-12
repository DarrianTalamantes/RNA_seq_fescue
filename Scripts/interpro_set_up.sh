version_major=5.68
version_minor=100.0
CONDA_PREFIX=/home/drt83172/.conda/envs/interproscan

wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${5.68}-${version_minor}/interproscan-${version_major}-${version_minor}-64-bit.tar.gz.md5
wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${version_major}-${version_minor}/interproscan-${version_major}-${version_minor}-64-bit.tar.gz
md5sum -c interproscan-${version_major}-${version_minor}-64-bit.tar.gz.md5
tar xvzf interproscan-${version_major}-${version_minor}-64-bit.tar.gz
rm -rf $CONDA_PREFIX/share/InterProScan/data/
cp -r interproscan-${version_major}-${version_minor}/data $CONDA_PREFIX/share/InterProScan/
