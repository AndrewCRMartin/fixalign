#!/usr/bin/perl -s
# A quick and dirty program to extract SS information from DSSP files. 
# Doesn't work very well if there is a chain break!
#
# Usage: dssp2ss.pl [-keep] file.dssp > file.ss
#        -keep keeps detailed SS assignments (G, B, S, etc). Default is
#              to convert to 3-state assignments
#
# (c) 2010, Andrew Martin, UCL
#
# V1.0  25.02.10  Original
#
use strict;

my (@chain, @aa, @ss, $data);
my $indata = 0;

# Read the data from the DSSP file
while(<>)
{
    # See if we've reached the comment preceeding the data section
    if(/^  #/)
    {
        $indata = 1;
    }
    elsif($indata)              # If in the data section 
    {
        if(substr($_,13,2) ne "!*")
        {
            # Store the chain
            push @chain, substr($_,11,1);

            # Store the amino acid
            $data = substr($_,13,1);
            $data = 'C' if(islower($data)); # Deal with disulphides
            push @aa,   $data;

            # Store the SS assignment
            $data = substr($_,16,1);
            if(defined($::keep))
            {
                $data = "C" if($data eq " ");
            }
            else
            {
                $data = "C" if(($data eq " ")||($data eq "T")||
                               ($data eq "S")||($data eq "B"));
                $data = "H" if(($data eq "G"));
            }
            push @ss,   $data;
        }
    }
}

# Write the data for each chain
my $thechain = "XXX";           # Used to identify a chain change
my $start    = 0;               # Used to record the start of each chain
for(my $i=0; $i<@chain; $i++)
{
    # If chain has changed...
    if($chain[$i] ne $thechain)
    {
        $thechain = $chain[$i];
        # If it's not just the start of the first chain...
        if($i != 0)
        {
            # Finish printing sequence data
            print "\n";

            # Print SS data
            for(my $j=$start; $j<$i; $j++)
            {
                print $ss[$j];
            }
            print "\n";
        }
        $start = $i;
        # Print start of chain
        print ">Chain_$thechain\n";
    }
    # Print this amino acid
    print $aa[$i];
}

# Print the last chain's ss data
print "\n";
for(my $j=$start; $j<@chain; $j++)
{
    print $ss[$j];
}
print "\n";

##########################################################################
# Test if a string starts with a lower case letter
sub islower
{
    my($c) = @_;
    return(1) if($c =~ /^[a-z]/);
    return(0);
}
