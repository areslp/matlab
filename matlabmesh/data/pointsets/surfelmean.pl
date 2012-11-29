

open FILE,$ARGV[0];

$sumx = 0;
$sumy = 0;
$sumz = 0;
$count = 0;
while (<FILE>) {
    $line = $_;
    @data = split /\s+/,$line;
    $sumx += $data[0];
    $sumy += $data[1];
    $sumz += $data[2];
    $count++;
}
$sumx /= $count;
$sumy /= $count;
$sumz /= $count;
print "$count surfels, mean $sumx $sumy $sumz\n";
