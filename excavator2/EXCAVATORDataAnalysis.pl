#!/usr/bin/env perl

use strict;
use warnings;

use Pod::Usage;
use Getopt::Long qw(:config gnu_getopt);
use File::Path;
use File::Basename;
use YAML::Tiny;
use POSIX;

######################################################################
#
#	Variables initialization with default values
#
######################################################################

my ($version, $verbose, $help, $man, $forced) = (0, 0, 0, 0, 0);
my ($targetFolder, $inputFolder, $outputFolder, $samplesFile, $experiment, $parametersFile, $threads);

my $excavatorPath = dirname(Cwd::abs_path($0));
### print "Program folder is: $excavatorPath\n";

$targetFolder = $excavatorPath;
$inputFolder = '.';
$outputFolder = '.';
$samplesFile = "$excavatorPath/samples.yaml";
$parametersFile = "$excavatorPath/parameters.yaml";
$threads = 1;

######################################################################
#
#	Reading user's options
#
######################################################################

GetOptions(
  'version' => \$version,
  'verbose|v' => \$verbose,
  'help|h' => \$help,
  'man|m' => \$man,
  'force|f' => \$forced,
  'input|i=s' => \$inputFolder,
  'samples|s=s' => \$samplesFile,
  'target|t=s' => \$targetFolder,
  'output|o=s' => \$outputFolder,
  'experiment|e=s' => \$experiment,
  'parameters|p=s' => \$parametersFile,
  'threads|@=s' => \$threads
) or pod2usage ("Error in command line");

$help and pod2usage (-verbose => 1, -exitval => 1, -output => \*STDOUT);
$man and pod2usage (-verbose => 2, -exitval => 1, -output => \*STDOUT);
#@ARGV == 0 or pod2usage ("Syntax error: the number of arguments found at command line is incorrect.");

open my $versionFile, '<', "$excavatorPath/version.txt" or die("Couldn't determine EXCAVATOR2 version.\n"); 
my $versionNumber = <$versionFile>; 
close $versionFile;

if ($version) {
  print "This is EXCAVATOR2 version $versionNumber\n";
  exit;
}

######################################################################
#
#	Checking options
#
######################################################################

die("Experiment mode $experiment is not a valid option.\n") unless $experiment =~ m/^(pooling|paired)$/i;
die("Target folder $targetFolder not found.\n") unless -d $targetFolder;
die("Input folder $inputFolder not found.\n") unless -d $inputFolder;
die("Samples file $samplesFile not found.\n") unless -f $samplesFile;
die("Parameters file $parametersFile not found.\n") unless -f $parametersFile;

my $settingsFile = "$inputFolder/settings.yaml";
die("Settings file $settingsFile not found.\n") unless (-f $settingsFile);

my $yaml = YAML::Tiny->read($settingsFile);
my $reference = $yaml->[0]->{'Reference'};
die("Centromeres coordinates file $reference->{'Centromeres'} not found.\n") unless -f $reference->{'Centromeres'};

if (-d $outputFolder) {
  if ($forced) {
    print STDERR "Overwriting output folder.\n- If you're trying to run two or more EXCAVATOR2 instances within the same output folder, results would be unpredictable.\n\n";
    rmtree($outputFolder);
  } else {
    die("Output folder '$outputFolder' already exists.\n");
  }
}

mkpath("$outputFolder/.tmp");
mkpath("$outputFolder/Results");
mkpath("$outputFolder/Plots");

######################################################################
#
#  When mode is pooling generates the Read Count for Controls
#
######################################################################

if ($experiment eq "pooling") { 
  print STDERR "Creating Pooling Control!\n";
  system(qq(R --slave --args "$excavatorPath" "$outputFolder" "$targetFolder" "$samplesFile" "$experiment" "$inputFolder" < "$excavatorPath/lib/R/PoolingCreateControl.R")) && die("Couldn't create Pooling Control");
  print STDERR "Pooling Control Created!\n";
}

######################################################################
#
#  Creating Multi Processor Analysis
#
######################################################################

my $Input_File_Parallel="$outputFolder/.tmp/ParallelAnalyzePerla.sh";

print STDERR "Preparing Multiprocessor Analysis...\n";
system(qq(R --slave --args "$excavatorPath" "$samplesFile" "$outputFolder" "$threads" "$experiment" "$inputFolder" "$targetFolder" "$parametersFile" < "$excavatorPath/lib/R/DataAnalysisParallel.R")) && die("Couldn't prepare jobs to execute");
print STDERR "Starting Multiprocessor Analysis!\n";

open(CHECKBOOK, "$Input_File_Parallel") || die "Couldn't open the input file $Input_File_Parallel.";
my @pids;
while(my $record = <CHECKBOOK>) {
  my $childpid = fork() or exec($record);
  push(@pids, $childpid);
}
close(CHECKBOOK);

#print STDERR "My Children: ", join(' ',@pids), "\n";
waitpid($_,0) for @pids;

rmtree("$outputFolder/.tmp");

$yaml->[2] = {
  'Module' => {
    'Name' => 'EXCAVATORDataAnalysis',
    'Version' => $versionNumber,
    'Run' => POSIX::strftime("%Y-%m-%d %H:%M:%S", localtime(time)),
    'Settings' => {
      'Mode' => $experiment,
      'Input' => $inputFolder,
      'Target' => $targetFolder,
      'Samples' => $samplesFile,
      'Output' => $outputFolder,
      'Threads' => $threads
    }
  }
};
$yaml->write("$outputFolder/settings.yaml");

print STDERR "Multiprocessor Analysis Complete!\n";




######################################################################
#
#  Documentation
#
######################################################################

=head1 SYNOPSIS 

 perl EXCAVATORDataAnalysis.pl [arguments] [options]
 
  Options:

       -h, --help                   Print help message.
       -m, --man                    Print complete documentation.
       -v, --verbose                Use verbose output.

 Function:
 
 EXCAVATORDataAnalysis.pl performs segmentation of the WMRC and classify each segmented region as one of 5 possible discrete states (2-copy deletion, 1-copy deletion, normal, 1-copy duplication and N-copy amplification).

 Example: 
 
 EXCAVATOR2> perl EXCAVATORDataAnalysis.pl ExperimentalFileAnalysis.w50K.txt --processors 6 --target MyTarget_w50K --assembly hg19 --output /.../OutEXCAVATOR2/Results_MyProject_w50K --mode pooling/paired

=head1 OPTIONS

=over 8

=item B<--help>

Print a brief usage message and detailed explanation of options.

=item B<--man>

Print the complete manual of the script.

=item B<--verbose>

Use verbose output.

=item B<--processors>

The number of thread to use for the analysis.

=item B<--output>

The output folder for resulting files.

=item B<--assembly>

The assembly exploited for read mapping and target initialization.

=item B<--target>

The "target name" used for target initialization with TargetPerla.pl.

=item B<--mode>

The experimental design mode to use. The possible options are "pooling" or "paired".

=back

=head1 DESCRIPTION

EXCAVATORDataAnalysis.pl perform the segmentation of the WMRC by means of the Shifting Level Model algorithm and exploits FastCall algorithm to classify each segmented region as one of the five possible discrete states (2-copy deletion, 1-copy deletion, normal, 1-copy duplication and N-copy amplification). The FastCall calling procedure takes into account sample heterogeneity and exploits the Expectation Maximization algorithm to estimate the parameters of a five gaussian mixture model and to provide the probability that each segment belongs to a specific copy number state.


EXCAVATOR2 is freely available to the community for non-commercial use. For questions or comments, please contact "romina.daurizio@gmail.com".
=cut
