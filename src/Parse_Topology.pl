#!/usr/bin/env perl
#Parsing multigenbank with bioperl
#LUDO 
#20180523



use Bio::SeqIO;

use Data::Dumper;


# get the file name, somehow 
$file = shift; 
$seqio_object = Bio::SeqIO->new(-file => $file,'-format' => 'genbank'); 
#my $seq_object = $seqio_object->next_seq;


while (my $seq = $seqio_object->next_seq) {
   #print $seq->display_id,"\n";
   #print $seq->primary_id,"\n";
#   print Dumper(\%t);
   #print values(%t),"\n";
   #print $seq->accession_number,"\n";
   #print ,"\n##########################################\n";
   # print $seq->is_circular,"\n";
   print ($seq->display_id."\t".$seq->accession_number."\t");
   print (($seq->is_circular==1)?"circular":"linear");
    print("\n");
}


