pipeline=$1
sed -n '2,$p' 2025-6-3_neuron_DEGs/DEGs."$pipeline"*edgeR.txt|awk '{print $6,$0}'|getLoci_inf.pl -i1 alluniq -i2 - -o -|awk '{if($5<0.05 && $2>0){print $1,"1"};if($5<0.05 && $2<0){print $1,"-1"}}' > "$pipeline".raw
cat "$pipeline".raw|awk '{print $1}'|sort > "$pipeline".raw.onlyID
cp "$pipeline".raw "$pipeline".raw2
cat "$pipeline".raw|awk '{print $1}'|sort|comm -13 - alluniq |awk '{print $1,"0"}' >> "$pipeline".raw2
sort -k1,1  "$pipeline".raw2 |awk 'BEGIN{print"geneID '$pipeline'"}{print}'> "$pipeline".txt_tmp
rm "$pipeline".raw
rm "$pipeline".raw2
rm "$pipeline".raw.onlyID
