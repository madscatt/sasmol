#!/usr/bin/perl

sub uniq2 {
    my %seen = ();
    my @r = ();
    foreach my $a (@_) {
        unless ($seen{$a}) {
            push @r, $a;
            $seen{$a} = 1;
        }
        else
        {
           #print "Duplicate ".$a."\n";
        }
    }
    return @r;
}

$file = $ARGV[$0];
open(FILE, $file);
@lines = <FILE>;


foreach $line (@lines)
{
   @words = split(/ /, $line);
   if ($words[0] =~ ATOM)
   {
      #print $line;
      push(@atoms, $words[1]);
   }
}

@uniqueatoms = &uniq2(@atoms);

foreach $atom (@uniqueatoms)
{
   print $atom."\n";
}
