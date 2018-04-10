#!/usr/bin/perl -w
use Data::Dumper;

die "perl $0 <inputdir> <outputdir> <barcode>\n" unless(@ARGV==3);
my $indir = $ARGV[0];
my $outdir = $ARGV[1];
my $barcode = $ARGV[2];
my @sample;
my @barcodes = split(/,/, $barcode);
foreach my $i (@barcodes){
    if($i =~ /-/){
        my @num = split("-",$i);
        push @sample, $num[0]..$num[1];
    }
    else{
        push @sample,$i;
    }
}

open(OUT, ">$outdir/result.xls");
print OUT "sample\tmapped reads\tmapping rate\tMAQ20\tMAQ30\terr rate\n";

foreach my $i(@sample){
    my @info;
    push @info,$i;
    open IN1, "$indir/$i/$i.bam.flagstat";
    open IN2, "$indir/$i/${i}_mapstat.xls";
    my $mapreads;
    my $maprate;
    my ($err, $maq20, $maq30);
    while(<IN1>){
    #    if(/(\d+) \+ \d+ mapped \((\d+)%.*\)/){
        if(/mapped/){
            chomp;
            my @line = split;
#            print Dumper(\@line),"\n";

            $mapreads = sprintf("%.2f", $line[0] / 1000000);
            $maprate = $line[4];
            $maprate =~ s/[(%]//g; 
            push @info, $mapreads, $maprate;
            last;
        }
    }
    close IN1;
    while(<IN2>){
        if(/Mismatch rate\t(\d+.\d+)%/){$err = $1;}
        if(/MAQ20\t(\d+.\d+)%/){$maq20 = $1;}
        if(/MAQ30\t(\d+.\d+)%/){$maq30 = $1;}
    }
    push @info, $maq20, $maq30, $err;
    close IN2;
    print OUT join("\t", @info),"\n";
}
close OUT;
