#!/usr/bin/perl

use strict;
use Getopt::Std;
use Time::Local;
use DBI;

getopts('h:u:p:d:o:s:', \ my %opts);

my ($dbhost, $dbuser, $dbpasswd, $dbname, $output_path, $datasource_id);

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

if ( defined($opts{o}) )
{
  $output_path=$opts{o};
}
else
{
  $output_path='./';
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

if( !open(LNCRNAFASTA, ">$output_path/lncrna.fa"))
{
  die "Can't create fasta: $!\n";
}

$sth = $dbh->prepare("select gene_id,flank,sequence from tb_gene where datasource_id=$datasource_id")
  or die "Can't prepare SQL statement: $DBI::errstr\n";
$sth->execute()  or die "Can't execute SQL statement: $DBI::errstr\n";
my $gene_count=0;
while( my($gene_id, $flank, $sequence) = $sth->fetchrow_array() )
{
  $sequence=substr($sequence,$flank,length($sequence)-$flank*2);
  print LNCRNAFASTA ">$gene_id\n";
  print LNCRNAFASTA "$sequence\n";
  $gene_count++;
  print "Writen gene records: $gene_count\r";
}
print "\n";

$sth->finish;

$dbh->disconnect;

close(LNCRNAFASTA);

sub usage
{
  print "-h [Database Hostname], default is 'localhost'.\n";
  print "-d [Database name], default is 'lncrna'.\n";
  print "-u [Database username], default is 'master'.\n";
  print "-p [Database password], default is null.\n";
  print "-o [Output result path].\n";
}

