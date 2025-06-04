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
my ($targetFolder, $samplesFile, $outputFolder, $MAPQ, $threads);

my $excavatorPath = dirname(Cwd::abs_path($0));
### print "Program folder is: $excavatorPath\n";

$outputFolder = '.';
$samplesFile = "samples-list.yaml";
$threads = 1;
$MAPQ = 0;

######################################################################
#
#  Reading user's options
#
######################################################################

GetOptions(
  'version' => \$version,
  'verbose|v' => \$verbose,
  'help|h' => \$help,
  'man|m' => \$man,
  'force|f' => \$forced,
  'target|t=s' => \$targetFolder,
  'samples|s=s' => \$samplesFile,
  'output|o=s' => \$outputFolder,
  'mapq|q=i' => \$MAPQ,
  'threads|@=i' => \$threads
) or pod2usage("Error in command line");

$help and pod2usage (-verbose => 1, -exitval => 1, -output => \*STDOUT);
$man and pod2usage (-verbose => 2, -exitval => 1, -output => \*STDOUT);

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

die("Threads must be greater than 0.\n") if ($threads < 1);

die("Target path $targetFolder not found.\n") unless (-d $targetFolder);

die("MAPQ must be greater than or equal to 0.\n") if ($MAPQ < 0);

die("Samples list file $samplesFile not found.\n") unless (-f $samplesFile);

my $settingsFile = "$targetFolder/settings.yaml";
die("Settings file $settingsFile not found.\n") unless (-f $settingsFile);

my $yaml = YAML::Tiny->read($samplesFile);
my %samples = %{$yaml->[0]};
for my $sample (keys %samples) {
  my $bam = $samples{$sample};
  if ($bam =~ /\.bam/) {
    die("Input BAM file $bam for sample $sample not found.\n") unless (-f $bam);
    die("Input BAM.BAI file $bam.bai for sample $sample not found.\n") unless (-f "$bam.bai");
  } elsif ($bam =~ /\.cram/) {
    die("Input CRAM file $bam for sample $sample not found.\n") unless (-f $bam);
    die("Input CRAM.CRAI file $bam.crai for sample $sample not found.\n") unless (-f "$bam.crai");
  } else {
    die("Input file $bam should be bam or cram.\n");
  }
}
print STDERR (scalar keys %samples), " samples to process.\n";

if (-d $outputFolder) {
  if ($forced) {
    print STDERR "Overwriting output folder.\n- If you're trying to run two or more EXCAVATOR2 instances within the same output folder, results would be unpredictable.\n\n";
    rmtree($outputFolder);
  } else {
    die("Output folder '$outputFolder' already exists.\n");
  }
}

######################################################################
#
#  Creating Multi Processor Analysis
#
######################################################################

### Creating Folder for MultiProcessor analysis ###
mkpath("$outputFolder/.tmp");
my $Input_File_Parallel = "$outputFolder/.tmp/ParallelReadPerla.sh";

print STDERR "Preparing Multiprocessor Analysis...\n";
system(qq(R --slave --args "$excavatorPath" "$samplesFile" "$outputFolder" "$targetFolder" "$threads" < "$excavatorPath/lib/R/DataPrepareParallel.R")) && die("Couldn't prepare jobs to execute");
print STDERR "Starting Multiprocessor Analysis!\n";

open(CHECKBOOK, "$Input_File_Parallel") || die "Couldn't open the input file $Input_File_Parallel.";
my @pids;
while(my $record = <CHECKBOOK>) {
  my $childpid = fork() or exec($record);
  push(@pids, $childpid);
}
close(CHECKBOOK);
#print STDERR "My Children: ", join(' ',@pids), "\n";
waitpid($_, 0) for @pids;

# Removing temp folder
rmtree("$outputFolder/.tmp");

$yaml = YAML::Tiny->read($settingsFile);
$yaml->[1] = {
  'Module' => {
    'Name' => 'EXCAVATORDataPrepare',
    'Version' => $versionNumber,
    'Run' => POSIX::strftime("%Y-%m-%d %H:%M:%S", localtime(time)),
    'Settings' => {
      'Target' => $targetFolder,
      'Samples' => $samplesFile,
      'Output' => $outputFolder,
      'MAPQ' => $MAPQ,
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

 perl EXCAVATORDataPrepare.pl [arguments] [options]

 Options:

       -h, --help                   Print help message.
       -m, --man                    Print complete documentation.
       -v, --verbose                Use verbose output.
           --mapq <integer>         Select mapping quality for .bam file filtering; if omitted default value is 0.
           
 Function:
 
 perl EXCAVATORDataPrepare.pl performs RC calculations, data normalization and data analysis on multiple .bam files.

 Example: 
 
 EXCAVATOR2> perl EXCAVATORDataPrepare.pl ExperimentalFilePrepare.w50000.txt --threads 6 --target MyTarget_w50000 --assembly hg19 


=head1 OPTIONS

=over 8

=item B<--help>

Print a brief usage message and detailed explanation of options.

=item B<--man>

Print the complete manual of the script.

=item B<--verbose>

Use verbose output.

=item B<--mapq>

Sets the numeric value of the mapping quality for .bam file filtering; must be an integer number. If omitted default value is 0.

=item B<--assembly>

The assembly exploited for read mapping and target initialization (hg19 or hg38).

=item B<--target>

The "target name" used for target initialization with TargetPerla.pl.


=back

=head1 DESCRIPTION

ReadPerla.pl is a Perl script which is part of the EXCAVATOR2 package. It performs RC calculations for In-target and Off-target regions, data normalization and data analysis on multiple .bam files.

The mapping quality value which is used by SAMtools can be set by means of the option --mapq when running ReadPerla.pl. If omitted default value is 0.

EXCAVATOR is freely available to the community for non-commercial use. For questions or comments, please contact "romina.daurizio@gmail.com".
=cut


