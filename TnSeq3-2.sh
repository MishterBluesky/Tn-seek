#!/usr/bin/env bash

usage () {
  echo "usage: $0 [-i <IR seq>] [-g <path to genome files>] <pfx> "
  echo "Required parameters:"
  echo "-i     This is your IR sequence"
  echo "-g     The location of the genome you're using"
  echo "must load these modules prior to using script: 
        python/2.7, cutadapt/1.8.1, bowtie2/2.3.2"
  echo "To load modules: module load python/2.7 , etc."
  echo ""
  echo "Example:"
  echo "$0 -i TATAAGAGTCAG -g \$HOME/ref_genome/PA14/PA14 condition1"
}

# Read options
while getopts "i:g:" option; do
  case "$option" in
    i) IR="$OPTARG" ;;
    g) GENOME="$OPTARG" ;;
    h) usage; exit 0 ;;
    :) echo "Error: -$option requires an argument"; usage; exit 1 ;;
    ?) echo "Error: unknown option -$option"; usage; exit 1 ;;
  esac
done
shift $(( OPTIND - 1 ))

if [ -z "$IR" ] || [ -z "$GENOME" ] || [ $# -lt 1 ]; then
  usage
  exit 1
fi

PREFIX=$1
BOWTIEREF=$GENOME
MAPFILE="${PREFIX}.fragment_map.tsv"

echo "Performing TnSeq analysis on $PREFIX..."
echo "TnSeq processing stats for $PREFIX" > $PREFIX-TnSeq.txt
echo "Total sequences: " >> $PREFIX-TnSeq.txt
egrep -c '^@' $PREFIX.trim.fastq >> $PREFIX-TnSeq.txt

# Mapping
echo "$PREFIX: Mapping with Bowtie2..."
echo "Bowtie2 report:" >> $PREFIX-TnSeq.txt
bowtie2 --end-to-end --very-sensitive -R 6 -p 16 -a -x $GENOME -U $PREFIX.trim.fastq -S $PREFIX.sam 2>> $PREFIX-TnSeq.txt

# Extract mapped fragments
cat $PREFIX.sam | grep -v '^@' | awk '$2 !~ /4/' | sort -u -k1,1 > $PREFIX-mapped.sam

echo "Number of reads mapping at high enough score:" >> $PREFIX-TnSeq.txt
cat $PREFIX-mapped.sam | wc -l >> $PREFIX-TnSeq.txt

# Collapse fragment hits to unique parent reads
if [ -f "$MAPFILE" ]; then
  echo "$PREFIX: Collapsing mapped fragments to parent reads..."

  cut -f1 $PREFIX-mapped.sam | sort | uniq > ${PREFIX}-mapped_fragments.txt

  join -1 1 -2 1 \
    <(sort ${PREFIX}-mapped_fragments.txt) \
    <(sort -k1,1 $MAPFILE) \
    | cut -f2 | sort | uniq > ${PREFIX}-mapped_parents.txt

  echo "Number of parent reads with mapped fragments:" >> $PREFIX-TnSeq.txt
  wc -l ${PREFIX}-mapped_parents.txt >> $PREFIX-TnSeq.txt
else
  echo "WARNING: Fragment map file $MAPFILE not found. Skipping parent deduplication." >> $PREFIX-TnSeq.txt
fi

# Directional site analysis
echo "$PREFIX: Tallying mapping results with directionality..."
grep -v '^@'  $PREFIX-mapped.sam | while IFS=$'\t' read -r f1 f2 f3 f4 f5 f6 f7 f8 f9 f10; do
  if ! [[ "$f2" =~ "100" ]]; then
    if ! [[ "$f2" =~ "10" ]]; then
      echo "$f4 $f2"
    else
      echo "$((f4 + ${#f10} - 2)) $f2"
    fi
  fi
done | grep '[0-9]' | sort | uniq -c | sort -n -r > $PREFIX-directional-sites.txt

# Position-only site analysis
echo "$PREFIX: Tallying mapping results..."
grep -v '^@' $PREFIX-mapped.sam | while IFS=$'\t' read -r f1 f2 f3 f4 f5 f6 f7 f8 f9 f10; do
  if ! [[ "$f2" =~ "100" ]]; then
    if ! [[ "$f2" =~ "10" ]]; then
      echo "$f4"
    else
      echo "$((f4 + ${#f10} - 2))"
    fi
  fi
done | grep '[0-9]' | sort | uniq -c | sort -n -r > $PREFIX-sites.txt

echo "Number of insertion sites identified:" >> $PREFIX-TnSeq.txt
wc -l $PREFIX-sites.txt >> $PREFIX-TnSeq.txt
echo "Most frequent sites:" >> $PREFIX-TnSeq.txt
head -10 $PREFIX-sites.txt >> $PREFIX-TnSeq.txt

# Cleanup
echo "$PREFIX: Cleaning up..."
mkdir $PREFIX 2> /dev/null
mv $PREFIX.trim.fastq $PREFIX/
mv $PREFIX-TnSeq.txt $PREFIX/
mv $PREFIX.sam $PREFIX/
mv $PREFIX-mapped.sam $PREFIX/
mv $PREFIX-directional-sites.txt $PREFIX/
mv $PREFIX-sites.txt $PREFIX/
[ -f "${PREFIX}-mapped_parents.txt" ] && mv ${PREFIX}-mapped_parents.txt $PREFIX/
[ -f "${PREFIX}-mapped_fragments.txt" ] && mv ${PREFIX}-mapped_fragments.txt $PREFIX/
