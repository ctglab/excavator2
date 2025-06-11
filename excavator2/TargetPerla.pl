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
my ($outputFolder, $settingsFile);

my $excavatorPath = dirname(Cwd::abs_path($0));
### print "Program folder is: $excavatorPath\n";

$outputFolder = $excavatorPath;
$settingsFile = 'config.yaml';

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
  'output|o=s' => \$outputFolder,
  'settings|s=s' => \$settingsFile
) or pod2usage("Error in command line");
### print("verbose: ${verbose}\n");
### print("help: ${help}\n");
### print("man: ${man}\n");
### print("window: ${windowSize}\n");
### print("target: ${targetName}\n");
### print("output: ${outputFolder}\n");
### print("parameters: ${parametersFile}\n");

$help and pod2usage(-verbose => 1, -exitval => 1, -output => \*STDOUT);
$man and pod2usage(-verbose => 2, -exitval => 1, -output => \*STDOUT);

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

die("Settings file $settingsFile not found.\n") unless (-f $settingsFile);

my $yaml = YAML::Tiny->read($settingsFile);
my $reference = $yaml->[0]->{'Reference'};
my $target = $yaml->[0]->{'Target'};

### print($yaml->write_string);
### print("\n");

die("Window size is too small.\n") if ($target->{'Window'} < 10);

$target->{'Name'} =~ s/^\s+|\s+$//g;
die("Must specify target name.\n") if ($target->{'Name'} eq '');

die("BigWig Genome file $reference->{'BigWig'} not found.\n") unless (-f $reference->{'BigWig'});

die("FASTA Genome file $reference->{'FASTA'} not found.\n") unless (-f $reference->{'FASTA'});

die("Target BED file $target->{'BED'} not found.\n") unless (-f $target->{'BED'});

die("Chromosomes coordinates file $reference->{'Chromosomes'} not found.\n") unless (-f $reference->{'Chromosomes'});

die("Gaps coordinates file $reference->{'Gaps'} not found.\n") unless (-f $reference->{'Gaps'});

my $finalOutputFolder = "$outputFolder/$reference->{'Assembly'}/$target->{'Name'}/w_$target->{'Window'}";
if (-d $finalOutputFolder) {
  if ($forced) {
    print STDERR "Overwriting output folder\n";
    rmtree($finalOutputFolder);
  } else {
    die("Output folder '$finalOutputFolder' already exists.\n");
  }
}

######################################################################
#
#	Checking target file format
#
######################################################################

print STDERR "Checking target file BED format... ";

# Test only first exon
my $bedTest = qx(head -1 "$target->{'BED'}" | awk '\$2 < \$3 { print "OK" }');

if ($bedTest eq "OK\n") {
  print STDERR "Target file seems properly formatted.\n";
} else {
  print STDERR "Exon #1 end value is smaller than start value.\nExons start (end) must be in column 2 (3).\n";
  die("Exiting: target seems to be uncorrectly formatted. Please check target BED file format!\n");
}

######################################################################
#
#	Target initialization
#
######################################################################

print STDERR "Filtering target from $target->{'BED'}\n";
system(qq(R --slave --args "$target->{'BED'}" "$finalOutputFolder" "$target->{'Name'}" "$reference->{'Assembly'}" "$target->{'Window'}" "$reference->{'Chromosomes'}" "$reference->{'Gaps'}" < "$excavatorPath/lib/R/FilterTarget.R")) && die("Couldn't filter targets");

print STDERR "Calculating Mappability and GC content...\n";
system(qq("$excavatorPath/lib/bash/TargetCreate.sh" "$reference->{'BigWig'}" "$finalOutputFolder/Filtered.txt" "$excavatorPath" "$reference->{'Assembly'}" "$target->{'Name'}" "$reference->{'FASTA'}" "$finalOutputFolder")) && die("Couldn't create target files");
#print "1)$reference->{'BigWig'}\n2)$outputFiltered\n3)$excavatorPath\n4)$reference->{'Assembly'}\n5)$target->{'Name'}\n6)$reference->{'FASTA'}\n7)$outputFolder";

$yaml->[0]->{'Module'} = {
  'Name' => 'TargetPerla',
  'Version' => $versionNumber,
  'Run' => POSIX::strftime("%Y-%m-%d %H:%M:%S", localtime(time)),
  'Settings' => {
    'Output' => "$outputFolder",
    'Settings' => "$settingsFile",
    'WorkDirectory' => Cwd::getcwd()
  }
};
$yaml->write("$finalOutputFolder/settings.yaml");
print STDERR "done!\n";

######################################################################
#
#	Documentation (to be updated)
#
######################################################################

=head1 SYNOPSIS 

 perl TargetPerla.pl [arguments] [options]

 Options:

       -h, --help                   Print help message.
       -m, --man                    Print complete documentation.
       -v, --verbose                Use verbose output.

 Function:
 
TargetPerla.pl initialises target files for further data processing with the EXCAVATOR2 package. It requires 5 arguments (one source files - with space-delimited paths to source data for mappability and GC-content calculations), path to target file, target name, window size and assembly to run properly. A sub-folder with the specified target name will be created under "EXCAVATOR2/data/targets/hgXX". Target input file (.bed, .txt or any plain text file) must be tab-delimited. 

 Example: 
 
 EXCAVATOR2> perl TargetPerla.pl SourceTarget.txt /Users/.../MyTarget.bed TargetName 50000  hg19
 
=head1 OPTIONS

=over 8

=item B<--help>

Print a brief usage message and detailed explanation of options.

=item B<--man>

Print the complete manual of the script.

=item B<--verbose>

Use verbose output.


=back

=head1 DESCRIPTION

TargetPerla.pl is a Perl script which is part of the EXCAVATOR2 package. It includes all of the first step operations of the EXCAVATOR2 package. It filters a target file and calculates, for a specific assembly, GC content and mappability.

It requires, as arguments, the path to a source file (the default source file is "SourceTarget.txt" which is placed in the main EXCAVATOR2 folder) containing the paths to source data (for the calculations of mappability and GC-content), the path to the target input file, a "target name", the window size and the assembly. Target input file (.bed, .txt or any plain text file) must be tab-delimited. Setting the target name as "MyTarget", all data calculated will be saved in the "MyTarget" folder in (if you are using the hg19 assembly) EXCAVATOR2/data/targets/hg19/MyTarget.

The allowed assemblies are hg19 and hg38.

EXCAVATOR2 is freely available to the community for non-commercial use. For questions or comments, please contact "romina.daurizio@gmail.com".

=cut

