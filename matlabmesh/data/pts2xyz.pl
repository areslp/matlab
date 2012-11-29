

open FILE,$ARGV[0];

while (<FILE>) {
    $line = $_;
    @data = split /\s+/,$line;
    print "$data[1] $data[2] $data[3]\n";
}    
