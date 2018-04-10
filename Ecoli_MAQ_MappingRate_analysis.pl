#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin $Script);
use Getopt::Long;

my ($indir,$outdir,$lane);
my ($barcode,@barcode, $length, $threads, $ref, $bwa, $samtools, $pid, $qid, $help);

##set software path
my $path;
my $qsub_sge = "$path/qsub-sge.pl";
my $java = "$path/java";
my $picard_dir = "$path/picard-tools-1.119";
my $fastqc = "$path/fastqc";
my $fqcheck = "$path/fqcheck33";

GetOptions(
  'i|input_dir:s' => \$indir,
  'o|output_dir:s' => \$outdir,
  'lane:s' => \$lane,
  'bc|barcode:s' => \$barcode,
  'l|read_length:i' => \$length,
  't|threads:i' => \$threads,
  'r|ref_fasta:s' => \$ref,
  'algn|algntools:s' => \$bwa,
  'samtools:s' => \$samtools,
  'h|help|?' => \$help);

## set default value
$barcode ||= "501-508,2382";
$length ||= 50;
$ref ||= "$path/Ecoli.fa";
$bwa ||= "$path/bwa";
$samtools ||= "$path/samtools";
$threads ||= 4;


#usage
my $usage = "
perl $0 -i <inputdir> -o <outputdir> -lane <lane>
-bc   STR   barcode,defualt: 501-508,2382
-l    INT   read length
-t    INT   threads
-r    STR   reference fasta
-h    help
";

##
if(!defined($indir) || !defined($outdir) || !defined($lane)){ die "$usage\n"}

##scripts
`mkdir -p $outdir/bin/shell_tmp`;
open(FQQC, ">", "$outdir/bin/shell_tmp/s1.fqqc.sh");
open(ALIGN, ">", "$outdir/bin/shell_tmp/s2.align.sh");
open(STAT, ">", "$outdir/bin/shell_tmp/s3.stat.sh");
open(RUN, ">", "$outdir/bin/shell_tmp/run.1.pl");
open(MAIN, ">", "$outdir/bin/shell_tmp/main.sh");


if (!defined($barcode)){
  @barcode = (501..508,2382);
}
else{
  my @temp = split(/,/, $barcode);
  foreach my $i (@temp){
      if($i =~ /-/){
          my @num = split("-",$i);
          push @barcode, $num[0]..$num[1];
      }
      else{
          push @barcode,$i;
      }
  }
}

foreach my $i (@barcode){
    my $od = "$outdir/align/$i";
    `mkdir -p $od`;
    `mkdir -p $outdir/fqqc/$i`;
    my $fq1 = `ls $indir/${lane}_*_${i}.fq.gz`;
    chomp $fq1;

    print FQQC "${fastqc} -f fastq $fq1;$fqcheck -r $fq1 -c ${outdir}/fqqc/$i/1.fqcheck;perl $path/fqcheck_distribute.pl ${outdir}/fqqc/$i/1.fqcheck -o ${outdir}/fqqc/$i/${lane}_${i}.\n";
    print ALIGN "${bwa} aln -e $length -t $threads $ref $fq1 > ${od}/${lane}_${i}.sai; ${bwa} samse ${ref} ${od}/${lane}_${i}.sai $fq1 | ${samtools} view -@ $threads -buhS -t ${ref}.fai - | ${samtools} sort -@ $threads -m 1000000000 -T ${od}/hg.sort -o ${od}/${lane}_${i}.1.sorted.bam -; ${java} -jar ${picard_dir}/MarkDuplicates.jar I=${od}/${lane}_${i}.1.sorted.bam O=${od}/${lane}_${i}.rmdup.bam M=${od}/${lane}_${i}.2.rmdup.mat REMOVE_DUPLICATES=false;${samtools} index ${od}/${lane}_${i}.rmdup.bam\n";
    print STAT "$samtools flagstat ${od}/${lane}_${i}.rmdup.bam >${od}/${i}.bam.flagstat;perl $Bin/align.stat.pl ${od}/${lane}_${i}.rmdup.bam ${od}/${i}_mapstat.xls\n";
}
close FQQC;
close ALIGN;
close STAT;

print RUN "#!/usr/bin/perl\nuse strict;\nuse threads;
  my \$step1 = async{
    system(\"perl $qsub_sge -resource=\\\"num_proc=1,vf=1G -P $pid -q $qid\\\" --jobprefix fqqc --convert no -reqsub $outdir/bin/shell_tmp/s1.fqqc.sh\");
  };
  my \$step2 = async{
    system(\"perl $qsub_sge -resource=\\\"num_proc=1,vf=8G -P $pid -q $qid\\\" --jobprefix align --convert no -reqsub $outdir/bin/shell_tmp/s2.align.sh\");
    system(\"perl $qsub_sge -resource=\\\"num_proc=1,vf=8G -P $pid -q $qid\\\" --jobprefix stat --convert no -reqsub $outdir/bin/shell_tmp/s3.stat.sh\");
  };
  if(\$step1->join() == 0 && \$step2->join() == 0){
    exit 0;
  }";

close RUN;
print MAIN "perl $outdir/bin/shell_tmp/run.1.pl\n";
print MAIN "perl $Bin/print_result.pl $outdir/align $outdir $barcode\n";
close MAIN;



