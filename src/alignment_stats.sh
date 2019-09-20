#!/usr/bin/env bash
RATES=${1?Error: UniquelyMappingRates path not given}
READS=${2?Error: UniquelyMappingReads path not given}


# To get uniq mapping rates to txt file from STAR logs:
#if [ ! -f deliverables/UniquelyMappingRates.txt ]; then
      echo -e 'sample\tUniquelyMappingRates' > $RATES
      echo -e 'sample\tUniquelyMappingReads' > $READS
      grep 'Uniquely mapped reads %' analysis/star/*Log.final.out |cut -d ' ' -f1,29|awk -F 'star/|.Log.final.out: \\|' '{print $2$3}' >> $RATES
      grep 'Uniquely mapped reads number' analysis/star/*Log.final.out |cut -d ' ' -f1,24|awk -F 'star/|.Log.final.out: \\|' '{print $2$3}' >> $READS
#fi
