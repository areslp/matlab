

open FILE,$ARGV[0];

while (<FILE>) {
    $line = $_;
    @data = split /\s+/,$line;
    print "vn $data[3] $data[4] $data[5]\n";
    print "vc $data[6] $data[7] $data[8]\n";
    print "v $data[0] $data[1] $data[2]\n";
}    
