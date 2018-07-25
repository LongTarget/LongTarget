#!/usr/bin/perl

use strict;
use Getopt::Std;
use DBI;

getopts('h:u:p:d:s:', \ my %opts);

my ($dbhost, $dbuser, $dbpasswd, $dbname, $db_file, $datasource_id, $species_id);

if ( defined($opts{h}) )
{
  $dbhost=$opts{h};
}
else
{
  $dbhost="localhost";
}

if ( defined($opts{u}) )
{
  $dbuser=$opts{u};
}
else
{
  $dbuser="master";
}

if ( defined($opts{p}))
{
  $dbpasswd=$opts{p};
}
else
{
  $dbpasswd='';
}

if ( defined($opts{d}))
{
  $dbname=$opts{d};
}
else
{
  $dbname='lncrnadb';
}

if ( defined($opts{s}))
{
  $datasource_id=$opts{s};
}
else
{
  $datasource_id=1;
}

my $dbh = DBI->connect( "dbi:mysql:dbname=$dbname;host=$dbhost", $dbuser, $dbpasswd);
my ($sth1,$sth2);
if ( !defined $dbh )
{
  die "Cannot connect to database  $dbhost $dbname $dbuser.\n";
}

my ($transposon_count,$position_status)=0;
my ($transposon_id, $seq_start, $seq_end, $gene_id);

$sth1=$dbh->prepare("select transposon_id,seq_start,seq_end,gene_id from tb_transposon where datasource_id=$datasource_id")
  or die "Can't prepare SQL statement: $DBI::errstr\n";
$sth1->execute();
$transposon_count=0;
while (($transposon_id, $seq_start, $seq_end, $gene_id) = $sth1->fetchrow_array())
{
  my ($gene_start,$gene_end);
  my ($exon_start,$exon_end);
  $transposon_count++;
  $position_status=2;

  $sth2=$dbh->prepare("select start,end from tb_gene where gene_id=$gene_id")
    or die "Can't prepare SQL statement: $DBI::errstr\n";
  $sth2->execute();
  if (!(($gene_start,$gene_end) = $sth2->fetchrow_array()))
  {
  	$sth2->finish;
  	print "Can't not find the gene where gene_id=$gene_id\n";
  	next;
  }
  $sth2->finish;
  
  $sth2=$dbh->prepare("select start,end from tb_myexon where gene_id=$gene_id")
    or die "Can't prepare SQL statement: $DBI::errstr\n";
  $sth2->execute();
  while (($exon_start,$exon_end) = $sth2->fetchrow_array())
  {
  	$exon_start=$exon_start-$gene_start;
  	$exon_end=$exon_end-$gene_start;
  	if( ($seq_start>=$exon_start) && ($exon_end>=$seq_end))
  	{
  	  $position_status=1;
  	  last;
  	}
  	elsif( ($seq_start>=$exon_start) && ($exon_end<$seq_end) )
  	{
  	  $position_status=3;
  	  last;
  	}
  }
  $sth2->finish;
  
  $sth2=$dbh->prepare("update tb_transposon set position_status=$position_status where transposon_id=$transposon_id")
    or die "Can't prepare SQL statement: $DBI::errstr\n";
  $sth2->execute();
  $sth2->finish;
  print "update transposon status count:\t$transposon_count\r";
}
$sth1->finish;

$dbh->disconnect;
print "\n";

sub usage
{
  print "-h [Database Hostname], default is 'localhost'.\n";
  print "-d [Database name], default is 'lncrna'.\n";
  print "-u [Database username], default is 'master'.\n";
  print "-p [Database password], default is null.\n";
}

