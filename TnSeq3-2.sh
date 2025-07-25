#!/usr/bin/env bash

usage () {
  echo "usage: $0 [-i <IR seq>] [-g <path to genome files>] <prefix>"
  echo "Required parameters:"
  echo "-i     IR sequence (not used in this version)"
  echo "-g     Path to genome index for Bowtie2"
  echo ""
  echo "Example:"
  echo "$0 -i TATAAGAGTCAG -g \$HOME/ref_genome/PA14/PA14 sample1"
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
TRIMM="${PREFIX}.trimm.fastq"
TRIMMED="${PREFIX}.trim.fastq"
COLLAPSED="${PREFIX}.collapsed.trim.fastq"
MAPFILE="${PREFIX}.fragment_map.tsv"
ALLSAM="${PREFIX}.all.sam"
BESTSAM="${PREFIX}.best.sam"

echo "Performing TnSeq analysis on $PREFIX..."
echo "TnSeq processing stats for $PREFIX" > "$PREFIX-TnSeq.txt"

# Count parent reads
echo -n "Total parent reads: " >> "$PREFIX-TnSeq.txt"
egrep -c '^@' "$TRIMM" >> "$PREFIX-TnSeq.txt"

# Count total 15-mers
echo -n "Total 15-mer fragments: " >> "$PREFIX-TnSeq.txt"
egrep -c '^@' "$TRIMMED" >> "$PREFIX-TnSeq.txt"

# Mapping all fragments
echo "$PREFIX: Mapping with Bowtie2..."
echo "Bowtie2 report:" >> "$PREFIX-TnSeq.txt"
bowtie2 --end-to-end --very-sensitive -R 6 -p 16 -a -x "$GENOME" -U "$TRIMMED" -S "$ALLSAM" 2>> "$PREFIX-TnSeq.txt"

# Collapse to best fragment per parent
echo "$PREFIX: Collapsing best fragment per parent..."
awk '!/^@/ {
  split($1, a, "_");
  parent = a[1];
  if (!seen[parent] || $5 > score[parent]) {
    seen[parent] = $0;
    score[parent] = $5;
  }
} END {
  for (p in seen) print seen[p];
}' "$ALLSAM" > "$BESTSAM"

# Save best fragment names
awk '{print $1}' "$BESTSAM" > "$PREFIX.best_fragments.txt"

# Extract fragment reads from original trimmed FASTQ
grep -A 3 -F -w -f "$PREFIX.best_fragments.txt" "$TRIMMED" | grep -v '^--$' > "$COLLAPSED"

# Remap collapsed FASTQ
echo "$PREFIX: Remapping collapsed reads..."
bowtie2 --end-to-end --very-sensitive -R 6 -p 16 -a -x "$GENOME" -U "$COLLAPSED" -S "$PREFIX.sam" 2>> "$PREFIX-TnSeq.txt"

# Extract mapped reads
grep -v '^@' "$PREFIX.sam" | awk '$2 != 4' | sort -u -k1,1 > "$PREFIX-mapped.sam"

# Report aligned reads
echo -n "Number of reads mapping at high enough score: " >> "$PREFIX-TnSeq.txt"
wc -l < "$PREFIX-mapped.sam" >> "$PREFIX-TnSeq.txt"

# Tally insertion sites with direction
grep -v '^@' "$PREFIX-mapped.sam" | while IFS=$'\t' read -r f1 f2 f3 f4 f5 f6 f7 f8 f9 f10; do
  if ! [[ "$f2" =~ "100" ]]; then
    if ! [[ "$f2" =~ "10" ]]; then
      echo "$f4 $f2"
    else
      echo "$((f4 + ${#f10} - 2)) $f2"
    fi
  fi
done | grep '[0-9]' | sort | uniq -c | sort -n -r > "$PREFIX-directional-sites.txt"

# Tally position-only insertion sites
grep -v '^@' "$PREFIX-mapped.sam" | while IFS=$'\t' read -r f1 f2 f3 f4 f5 f6 f7 f8 f9 f10; do
  if ! [[ "$f2" =~ "100" ]]; then
    if ! [[ "$f2" =~ "10" ]]; then
      echo "$f4"
    else
      echo "$((f4 + ${#f10} - 2))"
    fi
  fi
done | grep '[0-9]' | sort | uniq -c | sort -n -r > "$PREFIX-sites.txt"

# Report site counts
echo -n "Number of insertion sites identified: " >> "$PREFIX-TnSeq.txt"
wc -l < "$PREFIX-sites.txt" >> "$PREFIX-TnSeq.txt"
echo "Most frequent sites:" >> "$PREFIX-TnSeq.txt"
head -10 "$PREFIX-sites.txt" >> "$PREFIX-TnSeq.txt"

# Cleanup
echo "$PREFIX: Organizing output..."
mkdir -p "$PREFIX"
mv "$TRIMMED" "$PREFIX/"
mv "$TRIMM" "$PREFIX/"
mv "$ALLSAM" "$PREFIX/"
mv "$BESTSAM" "$PREFIX/"
mv "$COLLAPSED" "$PREFIX/"
mv "$PREFIX.sam" "$PREFIX/"
mv "$PREFIX-mapped.sam" "$PREFIX/"
mv "$PREFIX-sites.txt" "$PREFIX/"
mv "$PREFIX-directional-sites.txt" "$PREFIX/"
mv "$PREFIX-TnSeq.txt" "$PREFIX/"
mv "$PREFIX.best_fragments.txt" "$PREFIX/"
