#!/usr/bin/perl

use strict;
use Getopt::Std;
use DBI;

getopts('h:u:p:d:f:s:', \ my %opts);

my ($dbhost, $dbuser, $dbpasswd, $dbname, $repeatmasker_file, $datasource_id);

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

if ( defined($opts{f}) )
{
  $repeatmasker_file=$opts{f};
}
else
{
  &usage();
  exit 1;
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
my $sth;
if ( !defined $dbh )
{
  die "Cannot connect to database  $dbhost $dbname $dbuser.\n";
}

if ($sth=$dbh->prepare("delete from tb_transposon where datasource_id=$datasource_id"))
{
  $sth->execute();
}
else
{
  $dbh->disconnect;
  die "Can't prepare SQL statement: $DBI::errstr\n";
}

if( !open(REPEATMASKERFILE, $repeatmasker_file))
{
  $dbh->disconnect;
  die "Can't open file: $!\n";
}

my $db_count=0;
while(<REPEATMASKERFILE>)
{
  if ( !/^\s*[0-9]+/)
  {
    next;
  }

  chomp;

  my ($score, $substitute, $delete, $insert, $gene_id, $seq_start, $seq_end, $seq_left, $complement, $transposon_name, $transposon_type, $poson_start, $poson_end, $poson_left, $masker_id, $higher) = split;

  if( ! ($transposon_type=~/^SINE|^LINE|^LTR|^DNA/) )
  {
    next;
  }

  if( $poson_start=~/\(\d+\)/ )
  {
    $poson_start=$poson_end;
    $poson_end=$poson_left;
  }

  if ($sth = $dbh->prepare( "insert tb_transposon (transposon_name,transposon_type,seq_start,seq_end,poson_start,poson_end,gene_id,datasource_id) value (\"$transposon_name\",\"$transposon_type\",$seq_start,$seq_end,$poson_start,$poson_end,$gene_id,$datasource_id)"))
  {
    if( $sth->execute())
    {
      $db_count++;
    }
    else
    {
      print "\nError in:$_\n";
      print "Can't prepare SQL statement: $DBI::errstr\n";
      next;
    }
  }
  else
  {
    print "\nCan't prepare SQL statement: $DBI::errstr\n";
    next;
  }

  print "Writen db records: $db_count\r";
}
print "\n";

$sth->finish;

$dbh->disconnect;

close(REPEATMASKERFILE);

sub usage
{
  print "-h [Database Hostname], default is 'localhost'.\n";
  print "-d [Database name], default is 'lncrna'.\n";
  print "-u [Database username], default is 'master'.\n";
  print "-p [Database password], default is null.\n";
  print "-f [RepeatMasker output file].\n";
}

