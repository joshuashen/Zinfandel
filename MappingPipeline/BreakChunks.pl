#!/usr/bin/perl
use File::Basename;

#sub BreakFile($ShortName, $OutputDir)  {
    mkdir "$ARGV[1]";
    my $length=0;
    my $nChunks=1;
    ($name, $path) = fileparse("$ARGV[0]");

    open (IN1, "$ARGV[0].read1.fastq") or die "Can't open file $ARGV[0].read1.fastq";
    open (IN2, "$ARGV[0].read2.fastq") or die "Can't open file $ARGV[0].read2.fastq";
    open (OUT1, ">$ARGV[1]$name.read1.$nChunks.fastq") or die "Can't open output file $ARGV[1]$name.read1.$nChunks.fastq";
    open (OUT2, ">$ARGV[1]$name.read2.$nChunks.fastq") or die "Can't open output file $ARGV[1]$name.read2.$nChunks.fastq";

    while ($l1=<IN1> and $l2=<IN2>)  {
            if ($length >=8000000) { 
	        $length = 0;
	        $nChunks += 1;
	        close (OUT1);
		close (OUT2);
	        open (OUT1, ">$ARGV[1]$name.read1.$nChunks.fastq") or die "Can't open output file $ARGV[1]$name.read1.$nChunks.fastq";
		open (OUT2, ">$ARGV[1]$name.read2.$nChunks.fastq") or die "Can't open output file $ARGV[1]$name.read2.$nChunks.fastq";
            }
	    print OUT1 "$l1";
	    print OUT2 "$l2";
	    $length += 1;
    }

    close (IN1);
    close (IN2);
    close (OUT1);
    close (OUT2);
    print "$nChunks";
#}
