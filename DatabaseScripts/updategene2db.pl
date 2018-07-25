#!/usr/bin/perl

use strict;
use Getopt::Std;
use DBI;
use Bio::EnsEMBL::Registry;

getopts('h:u:p:d:s:', \ my %opts);

my ($dbhost, $dbuser, $dbpasswd, $dbname, $datasource_id);

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
  $dbname='lncrna';
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
my ($sth1, $sth2);
if ( !defined $dbh )
{
  die "Cannot connect to database  $dbhost $dbname $dbuser.\n";
}

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);
my $slice_adaptor = $registry->get_adaptor( 'Homo_sapiens', 'Core', 'Slice' );
my $slice;

$sth1 = $dbh->prepare("select gene_id,seqname,start,end,strand,flank from tb_gene")
  or die "Can't prepare SQL statement: $DBI::errstr\n";
$sth1->execute()  or die "Can't execute SQL statement: $DBI::errstr\n";

my($gene_id, $seqname,$start, $end, $strand, $flank, $sequence, $chrom_num, $strand_num);
my $db_count=0;
while(($gene_id, $seqname, $start, $end, $strand, $flank) = $sth1->fetchrow_array())
{
  $chrom_num=substr($seqname,index($seqname,'chr')+3);
  if( $strand eq '+' )
  {
    $strand_num=1;
  }
  else
  {
    $strand_num=-1;
  }
  $slice = $slice_adaptor->fetch_by_region('chromosome', $chrom_num, $start-$flank, $end+$flank, $strand_num);
  if (defined($slice))
  {
    $sequence = $slice->seq();
  }

  if ($sth2 = $dbh->prepare( "update tb_gene set sequence=\'$sequence\' where gene_id=$gene_id"))
  {
    if( $sth2->execute())
    {
      $db_count++;
      print "Update sequence in DB: $db_count\r";
    }
    else
    {
      print "\nError in:$_\n";
      print "Can't execute SQL statement: $DBI::errstr\n";
      next;
    }
  }
  else
  {
    print "\nCan't prepare SQL statement: $DBI::errstr\n";
    next;
  }
}
print "\n";
$sth2->finish;
$sth1->finish;

$dbh->disconnect;

sub usage
{
  print "-h [Database Hostname], default is 'localhost'.\n";
  print "-d [Database name], default is 'lncrna'.\n";
  print "-u [Database username], default is 'master'.\n";
  print "-p [Database password], default is null.\n";
}

