DATE=`date +'%Y%m%d'`
for i in anthill.sdsu.edu edwards-data.sdsu.edu rambox phantome.org edwards-dna; do
	echo $i;
	ssh $i 'lastlog' > $i.$DATE.lastlog
done


python2.7 ~/EdwardsLab/bin/merge_last_logs.py -l anthill.sdsu.edu.$DATE.lastlog -l edwards-data.sdsu.edu.$DATE.lastlog -l rambox.$DATE.lastlog -l phantome.org.$DATE.lastlog -l edwards-dna.$DATE.lastlog > lastlog.$DATE

