#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use Table;

use lib "$ENV{FIREDIR}/SCRIPTS";



use Fire;


use strict;


use Getopt::Long;

if (@ARGV == 0) {
  die "Args: --seedfile --expfile\n";
}

my $seedfile = undef;
my $expfile  = undef;

GetOptions("seedfile=s" => \$seedfile,
           "expfile=s"  => \$expfile);


my $h_ref_seeds = Sets::getIndex($seedfile);

print STDERR "Got " . scalar(keys(%$h_ref_seeds)) . " seeds.\n";

my $filename = Sets::filename($expfile);

# summary
my $fsum    = "$expfile\_FIRE/DNA/$filename.summary";
die "$fsum does not exist\n" if (! -e $fsum);
my $a_ref_s = Fire::loadFireMotifSummaryArray($fsum);
my %motifs = ();
foreach my $r (@$a_ref_s) {
  #print Sets::jointab($r);

  if (defined($h_ref_seeds->{ $r->{SEED} })) {
    $motifs{ $r->{MOTIF} } = 1;
  }
}



my $ta = Table->new;
$ta->loadFile("$expfile\_FIRE/DNA/$filename.matrix");
my $a_ref = $ta->getArray();
my $r     = shift @$a_ref;
print     Sets::jointab($r);

foreach my $r (@$a_ref) {
  if (defined($motifs{$r->[0]})) {
    print Sets::jointab($r);
  }
}

