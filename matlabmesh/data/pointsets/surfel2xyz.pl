

open FILE,$ARGV[0];

while (<FILE>) {
    $line = $_;
    @data = split /\s+/,$line;
    print "$data[0] $data[1] $data[2] $data[3] $data[4] $data[5]\n";
}    
