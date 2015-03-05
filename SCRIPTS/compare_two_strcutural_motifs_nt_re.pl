use lib "$ENV{FIREDIR}/SCRIPTS";

use Table;
use Sets;
use Getopt::Long;
use strict;

my ($m1, $m2) = @ARGV ;

my $re1 = &seq_to_re($m1) ;
my $more = $re1;
$more =~ s/^\.+//;
$more =~ s/\.+$//;

# output WMed motif to tmp file
my $mo = Sets::myre2scanace($more);
my $mn1 = Sets::getTempFile("/tmp/m1");
open OUT, ">$mn1" or die "cannot open mot file for output.\n";
print OUT "Motif 1\n$mo";
close OUT;

my $re2 = &seq_to_re($m2) ;
$more = $re2;
$more =~ s/^\.+//;
$more =~ s/\.+$//;

# output WMed motif to tmp file
my $mo = Sets::myre2scanace($more);
my $mn2 = Sets::getTempFile("/tmp/m2");
open OUT, ">$mn2" or die "cannot open mot file for output.\n";
print OUT "Motif 1\n$mo";
close OUT;

my $s_todo = "$compareace $mn1 $mn2 -ss -simple ";
my $score = `$s_todo`;
chomp $score;
print $score ;


sub seq_to_re {
  my $S = shift @_ ;
  $S =~ s/N/./g ;
  $S =~ s/Y/[TC]/g ;
  $S =~ s/R/[AG]/g ;
  $S =~ s/K/[TG]/g ;
  $S =~ s/M/[AC]/g ;
  $S =~ s/S/[GC]/g ;
  $S =~ s/W/[AT]/g ;
  $S =~ s/B/[GTC]/g ;
  $S =~ s/D/[GAT]/g ;
  $S =~ s/H/[ACT]/g ;
  $S =~ s/V/[GCA]/g ;

  return $S ;
}
