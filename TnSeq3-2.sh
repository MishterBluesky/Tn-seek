#!/usr/bin/env bash

usage () {
  echo "usage: $0 [-i <IR seq>] [-g <path to genome files>] <pfx> "
  echo "Required parameters:"
  echo "-i     This is your IR sequence"
  echo "-g     The location of the genome you're using"
  echo "Must load modules prior to using script: python/2.7, cutadapt/1.8.1, bowtie2/2.3.2"
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
shift $((OPTIND - 1))

if [ -z "$IR" ] || [ -z "$GENOME" ] || [ $# -lt 1 ]; then
  usage
  exit 1
fi

PREFIX=$1
BOWTIEREF=$GENOME
MAPFILE="${PREFIX}.fragment_map.tsv"

echo "Performing TnSeq analysis on fragments $PREFIX..."
echo "TnSeq processing stats for $PREFIX" > $PREFIX-TnSeq.txt
echo "Total sequences in original trimmed FASTQ: " >> $PREFIX-TnSeq.txt
egrep -c '^@' $PREFIX.trim.fastq >> $PREFIX-TnSeq.txt

# Initial Mapping (optional, can keep for stats but downstream uses collapsed)
echo "$PREFIX: Mapping with Bowtie2..."
echo "Bowtie2 report:" >> $PREFIX-TnSeq.txt
bowtie2 --end-to-end --very-sensitive -R 6 -p 16 -a -x "$GENOME" -U "$PREFIX.trim.fastq" -S "$PREFIX.sam" 2>> $PREFIX-TnSeq.txt

# Extract mapped fragments from initial mapping (optional for stats)
grep -v '^@' $PREFIX.sam | sort -u -k1,1 > $PREFIX-mapped.sam
echo "Number of reads mapping (original):" >> $PREFIX-TnSeq.txt
wc -l < "$PREFIX-mapped.sam" >> $PREFIX-TnSeq.txt

# Collapse fragment hits to best parent if fragment map available
if [ -f "$MAPFILE" ]; then
  echo "$PREFIX: Collapsing mapped fragments to parent reads..."

  awk '{
    fragment = $1;
    score = $5 + 0;

    n = split(fragment, parts, ":");
    if (index(parts[3], "_") > 0) {
      split(parts[3], subparts, "_");
      parts[3] = subparts[1];
    }
    parent = parts[1];
    for (i = 2; i <= n; i++) {
      parent = parent ":" parts[i];
    }
    if (!(parent in best) || score > best[parent]) {
      best[parent] = score;
      best_read[parent] = fragment;
    }
  }
  END {
    for (p in best_read) print best_read[p];
  }' "$PREFIX-mapped.sam" > "$PREFIX.best_fragments.txt"

  echo "Number of unique best fragments for each read (should be same as cutadapt step): $(wc -l < "$PREFIX.best_fragments.txt")"

  # Create collapsed FASTQ from best fragments
  echo "$PREFIX: Generating collapsed trimmed FASTQ..."

  TRIMMED="${PREFIX}.trim.fastq"
  COLLAPSED="${PREFIX}.collapsed.trim.fastq"

  awk '{print "@" $1}' "$PREFIX.best_fragments.txt" > "$PREFIX.best_headers.txt"
  grep -A 3 -F -f "$PREFIX.best_headers.txt" "$TRIMMED" | grep -v '^--$' > "$COLLAPSED"

  echo "Created collapsed trimmed fastq to rerun bowtie on best fragments: $COLLAPSED"

  # Map collapsed reads
  echo "$PREFIX: Mapping collapsed reads with Bowtie2..."
  bowtie2 --end-to-end --very-sensitive -R 6 -p 16 -a -x "$GENOME" -U "$COLLAPSED" -S "${PREFIX}.collapsed.sam" 2>> $PREFIX-TnSeq.txt

  # Update stats for collapsed reads
  echo "Total sequences in collapsed FASTQ:" >> $PREFIX-TnSeq.txt
  egrep -c '^@' "$COLLAPSED" >> $PREFIX-TnSeq.txt

  echo "Number of reads mapping (collapsed):" >> $PREFIX-TnSeq.txt
  grep -v '^@' "${PREFIX}.collapsed.sam" | wc -l >> $PREFIX-TnSeq.txt

else
  echo "WARNING: Fragment map file $MAPFILE not found. Skipping parent deduplication." >> $PREFIX-TnSeq.txt
fi

# Use collapsed.sam for all site analysis (if exists), else fallback
if [ -f "${PREFIX}.collapsed.sam" ]; then
  FINAL_SAM="${PREFIX}.collapsed.sam"
else
  FINAL_SAM="${PREFIX}-mapped.sam"
fi

# Directional site analysis on collapsed reads
echo "$PREFIX: Tallying mapping results with directionality..."
grep -v '^@' "$FINAL_SAM" | while IFS=$'\t' read -r f1 f2 f3 f4 f5 f6 f7 f8 f9 f10; do
  if ! [[ "$f2" =~ "100" ]]; then
    if ! [[ "$f2" =~ "10" ]]; then
      echo "$f4 $f2"
    else
      echo "$((f4 + ${#f10} - 2)) $f2"
    fi
  fi
done | grep '[0-9]' | sort | uniq -c | sort -nr > "$PREFIX-directional-sites.txt"

# Position-only site analysis on collapsed reads
echo "$PREFIX: Tallying mapping results..."
grep -v '^@' "$FINAL_SAM" | while IFS=$'\t' read -r f1 f2 f3 f4 f5 f6 f7 f8 f9 f10; do
  if ! [[ "$f2" =~ "100" ]]; then
    if ! [[ "$f2" =~ "10" ]]; then
      echo "$f4"
    else
      echo "$((f4 + ${#f10} - 2))"
    fi
  fi
done | grep '[0-9]' | sort | uniq -c | sort -nr > "$PREFIX-sites.txt"

echo "Number of insertion sites identified:" >> $PREFIX-TnSeq.txt
wc -l < "$PREFIX-sites.txt" >> $PREFIX-TnSeq.txt
echo "Most frequent sites:" >> $PREFIX-TnSeq.txt
head -10 "$PREFIX-sites.txt" >> $PREFIX-TnSeq.txt

# Cleanup
echo "$PREFIX: Cleaning up..."
mkdir -p "$PREFIX"
mv "$PREFIX.trim.fastq" "$PREFIX/"
mv "$PREFIX-TnSeq.txt" "$PREFIX/"
mv "$PREFIX.sam" "$PREFIX/"
mv "$PREFIX-mapped.sam" "$PREFIX/"
mv "$PREFIX-directional-sites.txt" "$PREFIX/"
mv "$PREFIX-sites.txt" "$PREFIX/"
[ -f "$MAPFILE" ] && mv "$MAPFILE" "$PREFIX/"
[ -f "$PREFIX.best_fragments.txt" ] && mv "$PREFIX.best_fragments.txt" "$PREFIX/"
[ -f "$PREFIX.best_headers.txt" ] && mv "$PREFIX.best_headers.txt" "$PREFIX/"
[ -f "$COLLAPSED" ] && mv "$COLLAPSED" "$PREFIX/"
[ -f "${PREFIX}.collapsed.sam" ] && mv "${PREFIX}.collapsed.sam" "$PREFIX/"

echo "$PREFIX: Done!"
