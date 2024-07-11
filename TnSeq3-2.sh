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
  echo "The required parameters must precede the file prefix for your sequence file:"
  echo "  (e.g. if your sequence file is named condition1.fastq,"
  echo "   the prefix is \"condition1\")"
  echo ""
  echo "Example:"
  echo "$0 -i TATAAGAGTCAG -g $HOME/ref_genome/PA14/PA14 condition1"
}

# Read in the important options
while getopts "i:g:" option; do
  case "$option" in
          i)  IR="$OPTARG" ;;
          g)  GENOME="$OPTARG" ;;
    h)  # it's always useful to provide some help 
        usage
        exit 0 
        ;;
    :)  echo "Error: -$option requires an argument" 
        usage
        exit 1
        ;;
    ?)  echo "Error: unknown option -$option" 
        usage
        exit 1
        ;;
  esac
done    
shift $(( OPTIND - 1 ))

# Do some error checking to make sure parameters are defined
if [ -z "$IR" ]; then
  echo "Error: you must specify the Tn end sequence using -i"
  usage
  exit 1
fi

if [ -z "$GENOME" ]; then
  echo "Error: you must specify an assembly using -g"
  usage
  exit 1
fi

# Give the usage if there aren't enough parameters
if [ $# -lt 1 ] ; then
  echo "you must provide a file prefix for analysis"
  usage
  exit 1
fi

PREFIX=$1
BOWTIEREF=$GENOME

echo "Performing TnSeq analysis on $PREFIX..."
echo "TnSeq processing stats for $PREFIX" > $PREFIX-TnSeq.txt
echo "Total sequences: " >> $PREFIX-TnSeq.txt
egrep -c '^@' $PREFIX.fastq >> $PREFIX-TnSeq.txt


# IRs
echo "$PREFIX: Searching for reads with an IR..."
echo "$PREFIX: Removing primer and IR sequences..."
# Modify the -M flag depending on your sequence length. This is currently based on 75 bp reads.

# Map and convert - feel free to change bowtie2 parameters yourself
echo "$PREFIX: Mapping with Bowtie2..."
echo "Bowtie2 report:" >> $PREFIX-TnSeq.txt
bowtie2 --end-to-end --very-sensitive -R 6 -p 16 -a -x $GENOME -U $PREFIX.trim.fastq -S $PREFIX.sam 2>> $PREFIX-TnSeq.txt
grep '^@' $PREFIX.sam > $PREFIX-mapped.sam
cat $PREFIX.sam | grep -v '^@' | awk '$2 !~ /4/' | sort -u -k1,1 >> $PREFIX-mapped.sam

echo "Number of reads mapping at high enough score:" >> $PREFIX-TnSeq.txt
cat $PREFIX-mapped.sam | wc -l >> $PREFIX-TnSeq.txt

# Tallying mapping results (corrected the calls for reads mapped to the reverse strand)
echo "$PREFIX: Tallying mapping results with directionality..."
# If it's forward strand, print start site (before TA); if it's reverse, print start plus length - 2 (TA are designated before)
grep -v '^@'  $PREFIX-mapped.sam | while IFS=$'\t' read -r f1 f2 f3 f4 f5 f6 f7 f8 f9 f10; do
  if ! [[ "$f2" =~ "100" ]]; then
    if ! [[ "$f2" =~ "10" ]]; then
      echo "$f4 $f2"
    else
      echo "$((f4 + ${#f10} - 2)) $f2"
    fi
  fi
done | grep '[0-9]' | sort | uniq -c | sort -n -r > $PREFIX-directional-sites.txt

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

# Sort output, cleanup
echo "$PREFIX: Cleaning up..."
mkdir $PREFIX 2> /dev/null
mv $PREFIX.trim.fastq $PREFIX/
mv $PREFIX-TnSeq.txt $PREFIX/
mv $PREFIX.sam $PREFIX/
mv $PREFIX-mapped.sam $PREFIX/
mv $PREFIX-directional-sites.txt $PREFIX/
mv $PREFIX-sites.txt $PREFIX/
