DATE=$(date +%Y%m%d)
cd /home2/db/blast
mkdir nr_$DATE
cd nr_$DATE
ncftpget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr*
cat *md5 > all.md5
md5sum -c all.md5
for t in *.tar.gz; do echo $t; tar xf $t; done
cd /home2/db/blast
rm -f nr
ln -s nr_$DATE nr

mkdir nt_$DATE
cd nt_$DATE
ncftpget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt*
cat *md5 > all.md5
md5sum -c all.md5
for t in *.tar.gz; do echo $t; tar xf $t; done
cd /home2/db/blast
rm -f nt
ln -s nt_$DATE nt


