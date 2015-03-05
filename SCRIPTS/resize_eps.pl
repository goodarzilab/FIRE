use lib "$ENV{FIREDIR}/SCRIPTS";
use lib "$ENV{FIREDIR}/SCRIPTS/PostScript-Simple-0.07/lib";

use PostScript::Simple;

my $filename = shift;

if(length($filename) > 0) {

my $e = new PostScript::Simple::EPS(file => "$filename");
my $e_width = abs($e->width());
my $e_height = abs($e->height());

my $scale = 1;

my $height = $e_height;

while ($height > 5000) {
	$scale -= 0.05;
	$height = $scale * $e_height;
}

my $height = $scale * $e_height;
my $width = $scale * $e_width;

$p = new PostScript::Simple(xsize => $width, ysize => $height, colour => 1, units=>"pt");
    
$e->scale($scale);
$p->importepsfile("$filename", 0, 0, $width, $height);
$p->output("$filename");
}