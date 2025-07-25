#!/usr/bin/env bash

usage () {
  echo "usage: $0 [-i <IR seq>] [-g <path to genome files>] <pfx> "
  echo "Required parameters:"
  echo "-i     This is your IR sequence"
  echo "-g     The location of the genome you're using"
  echo ""
  echo "Example:"
  echo "$0 -i TATAAGAGTCAG -g \$HOME/ref_genome/PA14/PA14 condition1"
}

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
TRIMM="${PREFIX}.trimm.fastq"
TRIMMED="${PREFIX}.trim.fastq"
COLLAPSED="${PREFIX}.collapsed.trim.fastq"
COLLAPSED_SAM="${PREFIX}.collapsed.sam"

echo "Performing TnSeq analysis on $PREFIX..."
echo "TnSeq processing stats for $PREFIX" > $PREFIX-TnSeq.txt

# Total parent reads
PARENT_TOTAL=$(grep -c '^@' "$TRIMM")
echo "Total parent reads: $PARENT_TOTAL" >> $PREFIX-TnSeq.txt

# Total fragments (15mers)
echo -n "Total 15mer fragments: " >> $PREFIX-TnSeq.txt
grep -c '^@' "$TRIMMED" >> $PREFIX-TnSeq.txt

# Primary mapping
echo "$PREFIX: Mapping all fragments with Bowtie2..."
bowtie2 --end-to-end --very-sensitive -R 6 -p 16 -a -x "$GENOME" -U "$TRIMMED" -S "$PREFIX.raw.sam" 2>> $PREFIX-TnSeq.txt

# Extract best-scoring fragment per parent
echo "$PREFIX: Collapsing to best-scoring fragment per parent..."
grep -v '^@' "$PREFIX.raw.sam" | awk '
{
  split($1, parts, "_");
  parent = parts[1];
  score = $5;
  line = $0;
  if (score > best[parent]) {
    best[parent] = score;
    mapline[parent] = line;
  }
}
END {
  for (p in mapline) print mapline[p];
}
' > "$PREFIX.best.sam"

# Get fragment headers
cut -f1 "$PREFIX.best.sam" > "$PREFIX.best_fragments.txt"

# Reconstruct collapsed FASTQ from .trim.fastq
echo "$PREFIX: Reconstructing collapsed.fastq from best fragments..."
grep -A 3 -Ff "$PREFIX.best_fragments.txt" "$TRIMMED" | grep -v '^--$' > "$COLLAPSED"

# Remap only best fragments
echo "$PREFIX: Mapping collapsed reads with Bowtie2..."
bowtie2 --end-to-end --very-sensitive -R 6 -p 16 -a -x "$GENOME" -U "$COLLAPSED" -S "$COLLAPSED_SAM" 2>> $PREFIX-TnSeq.txt

# Filter mapped reads
grep -v '^@' "$COLLAPSED_SAM" | awk '$2 != 4' | sort -u -k1,1 > "$PREFIX.mapped.sam"

echo "Number of reads mapping at high enough score:" >> $PREFIX-TnSeq.txt
wc -l < "$PREFIX.mapped.sam" >> $PREFIX-TnSeq.txt

# Tally sites (directional)
echo "$PREFIX: Generating directional site list..."
grep -v '^@' "$PREFIX.mapped.sam" | while IFS=$'\t' read -r f1 f2 f3 f4 f5 f6 f7 f8 f9 f10; do
  if ! [[ "$f2" =~ "100" ]]; then
    if ! [[ "$f2" =~ "10" ]]; then
      echo "$f4 $f2"
    else
      echo "$((f4 + ${#f10} - 2)) $f2"
    fi
  fi
done | grep '[0-9]' | sort | uniq -c | sort -n -r > "$PREFIX-directional-sites.txt"

# Tally sites (position only)
echo "$PREFIX: Generating positional site list..."
grep -v '^@' "$PREFIX.mapped.sam" | while IFS=$'\t' read -r f1 f2 f3 f4 f5 f6 f7 f8 f9 f10; do
  if ! [[ "$f2" =~ "100" ]]; then
    if ! [[ "$f2" =~ "10" ]]; then
      echo "$f4"
    else
      echo "$((f4 + ${#f10} - 2))"
    fi
  fi
done | grep '[0-9]' | sort | uniq -c | sort -n -r > "$PREFIX-sites.txt"

# Site stats
echo "Number of insertion sites identified:" >> $PREFIX-TnSeq.txt
wc -l < "$PREFIX-sites.txt" >> $PREFIX-TnSeq.txt
echo "Most frequent sites:" >> $PREFIX-TnSeq.txt
head -10 "$PREFIX-sites.txt" >> $PREFIX-TnSeq.txt

# Parent read mapping stats
MAPPED_PARENTS=$(cut -f1 "$PREFIX.mapped.sam" | sed 's/_.*//' | sort | uniq | wc -l)
UNMAPPED_PARENTS=$((PARENT_TOTAL - MAPPED_PARENTS))
PCT_MAPPED=$(awk -v a=$MAPPED_PARENTS -v b=$PARENT_TOTAL 'BEGIN { printf "%.2f", (a/b)*100 }')
PCT_UNMAPPED=$(awk -v a=$UNMAPPED_PARENTS -v b=$PARENT_TOTAL 'BEGIN { printf "%.2f", (a/b)*100 }')

echo "Parent reads with one or more mapped fragments: $MAPPED_PARENTS ($PCT_MAPPED%)" >> $PREFIX-TnSeq.txt
echo "Parent reads with no mapped fragments: $UNMAPPED_PARENTS ($PCT_UNMAPPED%)" >> $PREFIX-TnSeq.txt

# Organize outputs
mkdir -p "$PREFIX"
mv $PREFIX.*.fastq "$PREFIX/" 2>/dev/null
mv $PREFIX.*.sam "$PREFIX/" 2>/dev/null
mv $PREFIX.*.txt "$PREFIX/" 2>/dev/null
mv $PREFIX-sites.txt "$PREFIX/"
mv $PREFIX-directional-sites.txt "$PREFIX/"
