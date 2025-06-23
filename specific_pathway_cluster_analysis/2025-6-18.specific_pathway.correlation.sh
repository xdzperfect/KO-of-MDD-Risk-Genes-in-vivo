#2025-6-3.specific_pathway analysis for cluster analysis
#
cd /home/xuzz01/mv/zls_project/2024-4-1-dual-sg-SPD/single2024-9-9/script_neurons/specific_pathway/prepare/2025-6-3_neuron_DEGs
cp /home/kongxr/kongxr/zls_project/analysis_from_zhb/analysis20250219/DEGs*.edgeR.txt ./
rename.pl 's/zls_//g' *

cd ../
#cutoff: pvalue< 0.05 and |logfc|>0
(base) -bash-4.2$ cat pipeline2.sh
pipeline=$1
sed -n '2,$p' 2025-6-3_neuron_DEGs/DEGs."$pipeline"*edgeR.txt|awk '{print $6,$0}'|getLoci_inf.pl -i1 alluniq -i2 - -o -|awk '{if($5<0.05 && $2>0){print $1,"1"};if($5<0.05 && $2<0){print $1,"-1"}}' > "$pipeline".raw
cat "$pipeline".raw|awk '{print $1}'|sort > "$pipeline".raw.onlyID
cp "$pipeline".raw "$pipeline".raw2
cat "$pipeline".raw|awk '{print $1}'|sort|comm -13 - alluniq |awk '{print $1,"0"}' >> "$pipeline".raw2
sort -k1,1  "$pipeline".raw2 |awk 'BEGIN{print"geneID '$pipeline'"}{print}'> "$pipeline".txt_tmp
rm "$pipeline".raw
rm "$pipeline".raw2
rm "$pipeline".raw.onlyID

##
#final LDA matched MDD risk genes
(base) -bash-4.2$ cat >gene44_part.list.txt
Dennd1a
Chmp3
Pacrg
Erbb4
Tank
Kirrel3
Nrxn1
Nkain2
Tox
Pclo
Snrk
Ascc3
Suds3
Kcnq5
Etv1
Dagla
Ap3b1
Zfp804a
Exoc4
Gigyf2
Cnnm2
Fhit
Fam120a
Ccser1
Usp3
Slc30a9
Sppl3
Cntln
Plcl2
B4galnt4
Bltp1
Galnt2
Luzp2
Foxp2
Ltbp3
Cstpp1
Baz2b
Tenm2
Rere
Rbks
Cdh9
Nlgn1
Lhpp
Astn2

##running pipeline2
cat gene44_part.list.txt|while read i; do sh pipeline2.sh $i; done &

mv *txt_tmp gene44_part
#
cd gene44_part/
mkdir mediant
ls *txt_tmp | while read i; do
    cat $i | awk '{print $2}' > mediant/"$i".reduce;
done
cat Dennd1a.txt_tmp |awk '{print $1}' > mediant/Dennd1a.txt
cd mediant/
ls *reduce|RowColumConvertion.pl -i - -o - -f " "|awk '{print "paste Dennd1a.txt "$0" |sed '\''s/\\s/,/g'\'' >final.csv"}'|sh


#R script see "D:\data_analysis\zls_project\0SPD_dual-sg_pertubseq_project\2025-Nature_genetics_respone_prepare\2025-6-3-specific_pathway\specific-pathway.r"