#!/usr/bin/env bash
RATES=${1?Error: UniquelyMappingRates path not given}
READS=${2?Error: UniquelyMappingReads path not given}
MATRIX=${3?Error: STARmatrix path not given}

# To get uniq mapping rates to txt file from STAR logs:
#if [ ! -f deliverables/UniquelyMappingRates.txt ]; then
      echo -e 'sample\tUniquelyMappingRates' > $RATES
      echo -e 'sample\tUniquelyMappingReads' > $READS
      grep 'Uniquely mapped reads %' star/*Log.final.out |cut -d ' ' -f1,29|awk -F 'star/|.Log.final.out: \\|' '{print $2$3}' >> $RATES
      grep 'Uniquely mapped reads number' star/*Log.final.out |cut -d ' ' -f1,24|awk -F 'star/|.Log.final.out: \\|' '{print $2$3}' >> $READS
#fi

cat > program.awk << EOF
BEGIN {FS = "\t"} # field separator
#{for (i = 2; i <= NF; i += 4) printf ("%s%c", \$i, i + 3 <= NF ? "\t" : "\n");}					# for unstranded
{for (i = 4; i <= NF; i += 4) printf ("%s%c", \$i, i + 3 <= NF ? "\t" : "\n");}					# for stranded
EOF

# get name of each sample and add to list called header.txt
find star/*ReadsPerGene.out.tab | awk -F 'star/|.ReadsPerGene' '{print $2}' |tr '\n' '\t' > header.txt && sed -i -e '$a\' header.txt && sed -i -e 's/^/\t/' header.txt
# paste field '1' (stranded read counts) of each 'ReadsPerGene.out.tab' starting with the 5th line and add to rownames (get all gene ids)
paste star/*ReadsPerGene.out.tab | tail -n +5 | cut -f 1 > rownames
# add gene counts for each sample to a new count.matrix
paste star/*ReadsPerGene.out.tab | tail -n +5 | awk -f program.awk > deliverables/joinedSTAR2pass.count.matrix
# add rownames (gene IDs) to the count.matrix
paste rownames deliverables/joinedSTAR2pass.count.matrix > deliverables/joinedSTAR2pass.withrownames.count.matrix
# add header with sample IDs to the top of the count.matrix
cat header.txt deliverables/joinedSTAR2pass.withrownames.count.matrix > $MATRIX

#remove temp files
rm header.txt \
rownames \
program.awk \
deliverables/joinedSTAR2pass.count.matrix \
deliverables/joinedSTAR2pass.withrownames.count.matrix \
