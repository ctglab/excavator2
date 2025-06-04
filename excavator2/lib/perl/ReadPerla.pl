#!/usr/bin/env perl

use strict;
use warnings;

use Pod::Usage;
use Getopt::Long;
use File::Path;
use File::Basename;
use YAML::Tiny;

######################################################################
#
#	Variables initialization with default values
#
######################################################################

my ($verbose, $help, $man) = (0, 0, 0);
my ($targetFolder, $MAPQ, $excavatorPath, $samplesFile, $outputFolder);
my ($sampleOutputFolder);

$MAPQ = 20; # Why not 0?

######################################################################
#
#	Reading user's options
#
######################################################################

GetOptions(
  'verbose|v' => \$verbose,
  'help|h' => \$help,
  'man|m' => \$man,
  'target|t=s' => \$targetFolder,
  'samples|s=s' => \$samplesFile,
  'output|o=s' => \$outputFolder,
  'mapq|q=i' => \$MAPQ
) or pod2usage ("Error in command line");

######################################################################
#
#	Defining system variables
#
######################################################################

$excavatorPath = dirname(dirname(dirname(Cwd::abs_path($0)))); # orribile
#print "Program folder is: $excavatorPath\n";

my $yaml = YAML::Tiny->read($samplesFile);
my %samples = %{$yaml->[0]};

for my $sample (keys %samples) {

  my $bam = $samples{$sample};

  # Check again for existence?
  if ($bam =~ /\.bam/) {
    die("Input BAM file $bam for sample $sample not found.\n") unless (-f $bam);
    die("Input BAM.BAI file $bam.bai for sample $sample not found.\n") unless (-f "$bam.bai");
  } elsif ($bam =~ /\.cram/) {
    die("Input CRAM file $bam for sample $sample not found.\n") unless (-f $bam);
    die("Input CRAM.CRAI file $bam.crai for sample $sample not found.\n") unless (-f "$bam.crai");
  } else {
    die("Input file $bam should be bam or cram.\n");
  }

  $sampleOutputFolder = "$outputFolder/$sample";
  mkpath ("$sampleOutputFolder/RC");
  mkpath ("$sampleOutputFolder/RCNorm");
  mkpath ("$sampleOutputFolder/Images");
  
  print STDERR "Working on sample $sample.\n";
  
  ######################################################################
  #
  #  Read count
  #
  ######################################################################
  
  print STDERR "Creating Read Count Data on sample $sample.\n";
  system(qq("$excavatorPath/lib/bash/FiltBam.sh" "$bam" "$MAPQ" "$sampleOutputFolder" "$excavatorPath" "$sample" "$targetFolder")) && die("Couldn't filter BAM/CRAM file for sample $sample");
  rmtree("$sampleOutputFolder/.tmp");
  print STDERR "Read Count on sample $sample done!\n";

  print STDERR "Normalizing Read Count Data on sample $sample.\n";
  system(qq(R --slave --args "$excavatorPath" "$sampleOutputFolder" "$sample" "$targetFolder" < "$excavatorPath/lib/R/EXCAVATORNormalizationExome.R")) && die("Couldn't normalize read count data for sample $sample");
  print STDERR "Normalization on sample $sample done!\n";

}

print STDERR "Performed RC calculations and Normalization on samples: \n";
print STDERR "$_\n" for (keys %samples);
