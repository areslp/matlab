$argc = $#ARGV + 1;

$infile = $ARGV[0];

$decrate = 4;
if ($argc > 1) {
    $decrate = $ARGV[1];
}
    

print "$infile\n";

open FILE,"$infile";

open DECFILE,"> decimated_$infile";
open RESTFILE,"> remaining_$infile";

$deccount = 0;
$restcount = 0;

$linecount = 0;
while (<FILE>) {
    $line = $_;
    $linecount++;
    if ($linecount % $decrate == 0) {
	print DECFILE $_;
	$deccount++;
    } else {
	print RESTFILE $_;
	$restcount++;
    }
}
print "$deccount in decimated file, $restcount in remaining file\n";    
