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
my ($targetFolder, $samplesFile, $experiment, $inputFolder, $outputFolder, $parametersFile);

my $excavatorPath = dirname(dirname(dirname(Cwd::abs_path($0)))); # orribile
#print "Program folder is: $excavatorPath\n";

$inputFolder = '.';
$outputFolder = '.';
$samplesFile = "$excavatorPath/samples.yaml";
$parametersFile = "$excavatorPath/parameters.yaml";

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
  'input|i=s' => \$inputFolder,
  'samples|s=s' => \$samplesFile,
  'output|o=s' => \$outputFolder,
  'experiment|e=s' => \$experiment,
  'parameters|p=s' => \$parametersFile
) or pod2usage ("Error in command line");

#@ARGV == 0 or pod2usage ("Syntax error: the number of arguments found at command line is incorrect.");

######################################################################
#
#	Defining system variables
#
######################################################################

die("Mode $experiment is not a valid option.\n") unless ($experiment =~ m/^(pooling|paired)$/i);
die("Target folder $targetFolder not found.\n") unless (-d $targetFolder);
die("Input folder $inputFolder not found.\n") unless (-d $inputFolder);
die("Samples file $samplesFile not found.\n") unless (-f $samplesFile);
die("Parameters file $parametersFile not found.\n") unless (-f $parametersFile);

my $settingsFile = "$inputFolder/settings.yaml";
die("Settings file $settingsFile not found.\n") unless (-f $settingsFile);

my $yaml = YAML::Tiny->read($settingsFile);
my $reference = $yaml->[0]->{'Reference'};
die("Centromeres coordinates file $reference->{'Centromeres'} not found.\n") unless -f $reference->{'Centromeres'};

$yaml = YAML::Tiny->read($samplesFile);
my %samples = %{$yaml->[0]};

for my $label (keys %samples) {
  my $sample = $samples{$label};
  if ($label =~ m/^T\d+$/i) { 
     mkpath("$outputFolder/Results/$sample");
     mkpath("$outputFolder/Plots/$sample");
  } else {
     die("Label '$label' is not valid.\n") unless ($label =~ m/^C\d+$/i);
  }
}

######################################################################
#
#	Data calculations
#
######################################################################

print STDERR "Starting Segmentation and Calling...\n";
system(qq(R --slave --args "$outputFolder" "$targetFolder" "$samplesFile" "$experiment" "$excavatorPath" "$reference->{'Assembly'}" "$inputFolder" "$parametersFile" "$reference->{'Centromeres'}" < "$excavatorPath/lib/R/EXCAVATORInferenceExome.R")) && die("Couldn't perform segmentation and calling.\n");
print STDERR "Segmentation and Calling Complete\n";

print STDERR "Starting Plots Generation...\n";
system(qq(R --slave --args "$outputFolder" "$samplesFile" "$reference->{'Assembly'}" < "$excavatorPath/lib/R/EXCAVATORPlotsExome.R")) && die("Couldn't generate plots.\n");
print STDERR "Plots Generation Complete\n";

print STDERR "Performed CNVs analysis on samples: \n";
print STDERR keys %samples, "\n";
